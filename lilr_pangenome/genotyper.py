from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .models import AlleleCall, HaplotypeSample, LocusAssignment
from .resources import LilrResources
from .utils import extract_fasta_region, get_logger, minimap2_to_paf, miniprot_to_gff3, write_fasta

_LOG = get_logger(__name__)

_FLANK = 1_500
_MIN_PAF_COVERAGE = 0.5

def run_genotyping(
    sample: HaplotypeSample,
    resources: LilrResources,
    *,
    threads: int = 4,
    force_rerun: bool = False,
) -> HaplotypeSample:
    if sample.failed:
        return sample

    if sample.lilr_contigs_fasta is None:
        sample.failed = True
        sample.failure_message = "Step 3: lilr_contigs_fasta not set - run extraction first."
        _LOG.error("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    if not sample.locus_assignments:
        sample.failed = True
        sample.failure_message = "Step 3: no locus assignments - run copy_number step first."
        _LOG.error("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    out_dir = sample.output_dir / "genotyping"
    out_dir.mkdir(parents=True, exist_ok=True)

    _LOG.info("[%s] Step 3: Extracting per-locus sequences...", sample.sample_hap_id)

    by_locus: dict[str, list[LocusAssignment]] = defaultdict(list)
    for asgn in sample.locus_assignments:
        by_locus[asgn.locus].append(asgn)

    for locus, assignments in sorted(by_locus.items()):
        calls_for_locus: list[AlleleCall] = []

        for copy_idx, asgn in enumerate(assignments):
            try:
                seq = extract_fasta_region(
                    sample.lilr_contigs_fasta,
                    asgn.contig_name,
                    max(0, asgn.query_start - _FLANK),
                    asgn.query_end + _FLANK,
                    strand=asgn.strand,
                )
            except Exception as exc:
                _LOG.warning(
                    "[%s] Could not extract %s copy %d from %s: %s - skipping.",
                    sample.sample_hap_id, locus, copy_idx + 1, asgn.contig_name, exc,
                )
                continue
            query_id = f"{sample.sample_hap_id}_{locus}_copy{copy_idx + 1}"
            query_rec = SeqRecord(seq, id=query_id, description="")

            query_fasta = out_dir / f"{query_id}.fasta"

            if not force_rerun and query_fasta.exists():
                _LOG.debug("[%s] Skipping %s - fasta exists.", sample.sample_hap_id, query_id)
            else:
                write_fasta([query_rec], query_fasta)

            paf_path = out_dir / f"{query_id}.paf"

            if not force_rerun and paf_path.exists():
                _LOG.debug("[%s] Skipping %s alignment - PAF exists.", sample.sample_hap_id, query_id)
            else:
                try:
                    minimap2_to_paf(
                        resources.combined_locus_reference,
                        query_fasta,
                        paf_path,
                        preset="asm5",
                        threads=threads,
                    )
                except Exception as exc:
                    _LOG.warning(
                        "[%s] minimap2 failed for %s: %s", sample.sample_hap_id, query_id, exc
                    )
                    paf_path.write_text("")

            confirmed_locus, identity, mismatches = _parse_best_paf_hit(
                paf_path, query_len=len(seq)
            )

            if confirmed_locus and confirmed_locus != locus:
                _LOG.warning(
                    "[%s] Locus reassignment %s -> %s (identity=%.4f) for copy %d",
                    sample.sample_hap_id, locus, confirmed_locus, identity, copy_idx + 1,
                )
                asgn.locus = confirmed_locus

            cdna_fasta    = out_dir / f"{query_id}_cdna.fasta"
            protein_fasta = out_dir / f"{query_id}_protein.fasta"

            if not force_rerun and cdna_fasta.exists() and protein_fasta.exists():
                _LOG.debug("[%s] Skipping cDNA/protein for %s - files exist.", sample.sample_hap_id, query_id)
            else:
                try:
                    cdna_seq, prot_seq = _extract_cdna_protein(
                        query_fasta,
                        asgn.locus,
                        resources,
                        query_id,
                        out_dir,
                    )
                    write_fasta([SeqRecord(cdna_seq, id=query_id, description="cDNA")], cdna_fasta)
                    write_fasta([SeqRecord(prot_seq, id=query_id, description="protein")], protein_fasta)
                except Exception as exc:
                    _LOG.warning("[%s] cDNA/protein extraction failed for %s: %s", sample.sample_hap_id, query_id, exc)
                    cdna_fasta    = None
                    protein_fasta = None

            call = AlleleCall(
                locus=asgn.locus,
                best_allele=f"{asgn.locus}_copy{copy_idx + 1}",
                mismatch_score=mismatches,
                identity_to_ref=identity,
                novel_allele=mismatches > 0,
                sequence_fasta=query_fasta,
                cdna_fasta=cdna_fasta,
                protein_fasta=protein_fasta,
                confirmed_locus=confirmed_locus or locus,
                locus_confirmed=(confirmed_locus == locus) if confirmed_locus else True,
            )
            calls_for_locus.append(call)
            _LOG.info(
                "[%s] %s copy%d: confirmed_locus=%s identity=%.4f mismatches=%d",
                sample.sample_hap_id, locus, copy_idx + 1,
                call.confirmed_locus, call.identity_to_ref, call.mismatch_score,
            )

        if not calls_for_locus:
            continue

        merged = calls_for_locus[0]
        for extra in calls_for_locus[1:]:
            merged.best_allele = f"{merged.best_allele}+{extra.best_allele}"
            merged.mismatch_score += extra.mismatch_score
            merged.identity_to_ref = (
                merged.identity_to_ref + extra.identity_to_ref
            ) / 2
            merged.novel_allele = merged.novel_allele or extra.novel_allele
            merged.locus_confirmed = merged.locus_confirmed and extra.locus_confirmed

        sample.allele_calls[locus] = merged

    return sample

def _extract_cdna_protein(
    gdna_fasta: Path,
    locus: str,
    resources: LilrResources,
    query_id: str,
    out_dir: Path,
) -> tuple[Seq, Seq]:
    from Bio import SeqIO as _SeqIO

    bait_fasta  = out_dir / f"{query_id}_bait.faa"
    gff3_path   = out_dir / f"{query_id}_cdna_miniprot.gff3"

    _write_locus_protein(resources.lilr_proteins, locus, bait_fasta)

    miniprot_to_gff3(
        target=gdna_fasta,
        proteins=bait_fasta,
        output_gff3=gff3_path,
        threads=1,
        max_hits=1,
        min_score_fraction=0.50,
    )

    cds_blocks = _parse_gff3_cds(gff3_path)
    if not cds_blocks:
        raise ValueError(f"miniprot found no CDS in {gff3_path}")

    records = list(_SeqIO.parse(str(gdna_fasta), "fasta"))
    if not records:
        raise ValueError(f"Empty gDNA fasta: {gdna_fasta}")
    gdna_seq = records[0].seq

    cdna = Seq("").join([gdna_seq[s:e] for s, e in sorted(cds_blocks)])
    protein = cdna.translate(to_stop=True)
    return cdna, protein

def _write_locus_protein(lilr_proteins: Path, locus: str, output: Path) -> None:
    from Bio import SeqIO as _SeqIO
    with open(lilr_proteins) as fh:
        records = [r for r in _SeqIO.parse(fh, "fasta") if r.id == locus]
    if not records:
        raise ValueError(f"No protein found for locus {locus} in {lilr_proteins}")
    write_fasta(records[:1], output)

def _parse_gff3_cds(gff3_path: Path) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            start = int(cols[3]) - 1
            end   = int(cols[4])
            blocks.append((start, end))
    return blocks

def _parse_best_paf_hit(
    paf_path: Path,
    query_len: int,
) -> tuple[str | None, float, int]:
    best_ref: str | None = None
    best_score: float = 0.0
    best_identity: float = 0.0
    best_mismatches: int = 0

    if not paf_path.exists() or paf_path.stat().st_size == 0:
        return best_ref, best_identity, best_mismatches

    if query_len == 0:
        return best_ref, best_identity, best_mismatches

    with open(paf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue

            ref_name   = cols[5]
            q_start    = int(cols[2])
            q_end      = int(cols[3])
            matches    = int(cols[9])
            aln_len    = int(cols[10])

            if aln_len == 0:
                continue

            q_cov = (q_end - q_start) / query_len
            if q_cov < _MIN_PAF_COVERAGE:
                continue

            score    = matches / query_len
            identity = matches / aln_len

            mismatches = aln_len - matches
            for tag in cols[12:]:
                if tag.startswith("NM:i:"):
                    try:
                        mismatches = int(tag[5:])
                    except ValueError:
                        pass
                    break

            if score > best_score:
                best_score      = score
                best_identity   = identity
                best_ref        = ref_name
                best_mismatches = mismatches

    return best_ref, best_identity, best_mismatches

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd

from .models import HaplotypeSample, LILR_LOCI, LocusAssignment, SamplePair
from .resources import LilrResources
from .utils import get_logger

_LOG = get_logger(__name__)

MIN_IDENTITY  = 0.75
MIN_COVERAGE  = 0.65
MAX_GENE_SPAN = 20_000

_LILRA6_CLUSTER_GAP        = 8_000
_LILRA6_MIN_COMBINED_SCORE = 0.90
_LILRA6_MIN_COPY_SPAN      = 3_000
_LILRA6_MIN_CLUSTER_BLOCKS = 3

_LOCUS_MAX_COPIES: dict[str, int] = {
    "LILRA3": 1,
    "LILRA6": 8,
}

_KNOWN_LILR_GENES: frozenset[str] = frozenset(LILR_LOCI)

def run_copy_number(
    sample: HaplotypeSample,
    resources: LilrResources,
    *,
    threads: int = 4,
    min_identity: float = MIN_IDENTITY,
    min_query_cov: float = MIN_COVERAGE,
    force_rerun: bool = False,
) -> HaplotypeSample:
    if sample.failed:
        return sample

    if sample.lilr_contigs_fasta is None:
        sample.failed = True
        sample.failure_message = "Step 2: lilr_contigs_fasta not set - run extraction first."
        _LOG.error("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    gff3_path = (
        sample.output_dir / "extractedFasta"
        / f"{sample.sample_hap_id}_lilr_miniprot.gff3"
    )

    if not gff3_path.exists():
        sample.failed = True
        sample.failure_message = (
            f"Step 2: miniprot GFF3 not found at {gff3_path}. "
            "Re-run extraction (force_rerun=True on run_extraction)."
        )
        _LOG.error("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    _LOG.info("[%s] Step 2: Parsing miniprot GFF3 for locus assignments...", sample.sample_hap_id)

    raw_hits = _parse_gff3_hits(gff3_path)
    _LOG.debug("[%s] %d raw mRNA hits from GFF3.", sample.sample_hap_id, len(raw_hits))

    accepted = _filter_hits(
        raw_hits,
        min_identity=min_identity,
        min_coverage=min_query_cov,
        sample_id=sample.sample_hap_id,
    )

    if not accepted:
        sample.failed = True
        sample.failure_message = (
            f"Step 2: No hits passed identity >= {min_identity:.0%}, "
            f"coverage >= {min_query_cov:.0%}, span <= {MAX_GENE_SPAN} bp, "
            "and fs=0, st=0 thresholds."
        )
        _LOG.warning("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    assignments: list[LocusAssignment] = []
    for h in accepted:
        assignments.append(LocusAssignment(
            contig_name=h["contig"],
            locus=h["gene"],
            identity=h["identity"],
            query_coverage=h["coverage"],
            ref_allele=h["gene"],
            query_start=h["start"],
            query_end=h["end"],
            strand=h["strand"],
            cds_coords=h.get("exons", []),
        ))

    sample.locus_assignments = assignments
    sample.locus_contig_map = _build_locus_contig_map(assignments)

    by_locus: dict[str, list[LocusAssignment]] = defaultdict(list)
    for asgn in assignments:
        by_locus[asgn.locus].append(asgn)

    sample.copy_numbers = {locus: len(lst) for locus, lst in by_locus.items()}
    for locus in LILR_LOCI:
        sample.copy_numbers.setdefault(locus, 0)

    _log_copy_numbers(sample)
    return sample

def write_copy_number_frames(
    sample_pairs: list[SamplePair],
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for pair in sample_pairs:
        row = {"sample": pair.name}
        row.update(pair.diploid_copy_numbers())
        rows.append(row)

    cn_df = pd.DataFrame(rows).set_index("sample")[LILR_LOCI]
    cn_df.to_csv(output_dir / "assemblyCopyNumberFrame.csv")
    _LOG.info("Wrote assemblyCopyNumberFrame.csv")

    asgn_rows: list[dict] = []
    for pair in sample_pairs:
        for hap in pair.haplotypes():
            for asgn in hap.locus_assignments:
                asgn_rows.append({
                    "sample": pair.name,
                    "haplotype": hap.hap_id,
                    "contig": asgn.contig_name,
                    "locus": asgn.locus,
                    "ref_allele": asgn.ref_allele,
                    "contig_start": asgn.query_start,
                    "contig_end": asgn.query_end,
                    "strand": asgn.strand,
                    "identity": round(asgn.identity, 4),
                    "coverage": round(asgn.query_coverage, 4),
                })

    asgn_df = pd.DataFrame(asgn_rows) if asgn_rows else pd.DataFrame()
    asgn_df.to_csv(output_dir / "assemblyContigAssignments.csv", index=False)
    _LOG.info("Wrote assemblyContigAssignments.csv")

    return cn_df, asgn_df

def _parse_paf_tags(paf_line: str) -> tuple[int, int]:
    fs = st = 0
    for field in paf_line.split("\t"):
        if field.startswith("fs:i:"):
            try:
                fs = int(field[5:])
            except ValueError:
                pass
        elif field.startswith("st:i:"):
            try:
                st = int(field[5:])
            except ValueError:
                pass
    return fs, st

def _parse_gff3_hits(gff3_path: Path) -> list[dict]:
    raw: list[dict] = []
    cur: dict | None = None
    exons: list[tuple[int, int]] = []
    pending_paf: str | None = None

    with open(gff3_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith("##PAF"):
                pending_paf = line
                continue
            if line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            ftype  = cols[2]
            contig = cols[0]
            start  = int(cols[3]) - 1
            end    = int(cols[4])
            strand = cols[6]
            attrs  = cols[8]

            if ftype == "mRNA":
                if cur is not None and exons:
                    raw.append({**cur, "exons": list(exons)})

                gene = None
                identity = 0.0
                coverage = 0.0
                for attr in attrs.split(";"):
                    attr = attr.strip()
                    if attr.startswith("Target="):
                        gene = attr[7:].split()[0]
                    elif attr.startswith("Identity="):
                        try:
                            identity = float(attr[9:])
                        except ValueError:
                            pass
                    elif attr.startswith("Positive="):
                        try:
                            coverage = float(attr[9:])
                        except ValueError:
                            pass

                fs, st = _parse_paf_tags(pending_paf) if pending_paf else (0, 0)

                cur = {
                    "gene":     gene,
                    "contig":   contig,
                    "strand":   strand,
                    "start":    start,
                    "end":      end,
                    "identity": identity,
                    "coverage": coverage,
                    "fs":       fs,
                    "st":       st,
                }
                exons = []
                pending_paf = None

            elif ftype == "CDS" and cur is not None:
                exons.append((start, end))

    if cur is not None and exons:
        raw.append({**cur, "exons": list(exons)})

    return raw

def _filter_hits(
    raw_hits: list[dict],
    *,
    min_identity: float,
    min_coverage: float,
    sample_id: str = "",
) -> list[dict]:
    stage1_other: list[dict] = []
    stage1_lilra6: list[dict] = []

    for h in raw_hits:
        gene = h.get("gene")
        if gene not in _KNOWN_LILR_GENES:
            continue
        if h["identity"] < min_identity:
            continue
        if h["coverage"] < min_coverage:
            continue

        if gene == "LILRA6":
            combined = h["identity"] * h["coverage"]
            if combined < _LILRA6_MIN_COMBINED_SCORE:
                _LOG.debug(
                    "[%s] Dropped LILRA6 on %s: identityxcoverage=%.3f (likely LILRB3 cross-hit)",
                    sample_id, h["contig"], combined,
                )
                continue
            stage1_lilra6.append(h)
        else:
            span = h["end"] - h["start"]
            if span > MAX_GENE_SPAN:
                _LOG.debug(
                    "[%s] Dropped %s on %s span=%d bp (cross-paralogue ghost)",
                    sample_id, gene, h["contig"], span,
                )
                continue
            if h["fs"] != 0:
                _LOG.debug(
                    "[%s] Dropped %s on %s: %d frameshift(s)",
                    sample_id, gene, h["contig"], h["fs"],
                )
                continue
            if h["st"] != 0:
                _LOG.debug(
                    "[%s] Dropped %s on %s: %d in-translation stop(s)",
                    sample_id, gene, h["contig"], h["st"],
                )
                continue
            stage1_other.append(h)

    _LOG.debug(
        "[%s] Stage 1: %d non-LILRA6, %d LILRA6 quality hits.",
        sample_id, len(stage1_other), len(stage1_lilra6),
    )

    stage1_other.sort(key=lambda h: -h["identity"])
    claimed: list[tuple[str, int, int]] = []

    def _overlaps(contig: str, s: int, e: int) -> bool:
        return any(ac == contig and s < ae and e > as_ for ac, as_, ae in claimed)

    stage2_other: list[dict] = []
    for h in stage1_other:
        c, s, e = h["contig"], h["start"], h["end"]
        if not _overlaps(c, s, e):
            claimed.append((c, s, e))
            stage2_other.append(h)
        else:
            _LOG.debug(
                "[%s] Dropped %s on %s [%d-%d]: overlaps higher-identity hit",
                sample_id, h["gene"], c, s, e,
            )

    stage2_lilra6 = _count_lilra6_copies(stage1_lilra6, sample_id=sample_id)

    all_accepted = stage2_other + stage2_lilra6

    _LOG.debug(
        "[%s] Stage 2: %d non-LILRA6 + %d LILRA6 copies.",
        sample_id, len(stage2_other), len(stage2_lilra6),
    )

    hits_per_locus: dict[str, int] = {}
    for h in all_accepted:
        hits_per_locus[h["gene"]] = hits_per_locus.get(h["gene"], 0) + 1

    accepted_loci = set(hits_per_locus)
    non_variable = accepted_loci - {"LILRA3"}
    if len(non_variable) < 3:
        _LOG.warning(
            "[%s] Only %d non-variable LILR loci detected (%s) - "
            "possible extraction failure or assembly gap.",
            sample_id, len(non_variable), ", ".join(sorted(non_variable)),
        )

    for locus, max_cn in _LOCUS_MAX_COPIES.items():
        n = hits_per_locus.get(locus, 0)
        if n > max_cn:
            _LOG.warning(
                "[%s] %s has %d copies - exceeds expected maximum of %d; "
                "check for cross-paralog ghost hits.",
                sample_id, locus, n, max_cn,
            )

    return all_accepted

def _count_lilra6_copies(
    hits: list[dict],
    *,
    cluster_gap: int = _LILRA6_CLUSTER_GAP,
    min_cluster_blocks: int = _LILRA6_MIN_CLUSTER_BLOCKS,
    min_copy_span: int = _LILRA6_MIN_COPY_SPAN,
    sample_id: str = "",
) -> list[dict]:
    if not hits:
        return []

    by_contig: dict[str, list[dict]] = defaultdict(list)
    for h in hits:
        for cds_s, cds_e in h.get("exons", []):
            by_contig[h["contig"]].append({"start": cds_s, "end": cds_e, "hit": h})

    copies: list[dict] = []

    for contig, blocks in by_contig.items():
        if not blocks:
            continue
        blocks.sort(key=lambda b: b["start"])

        clusters: list[list[dict]] = []
        cur: list[dict] = [blocks[0]]
        for blk in blocks[1:]:
            if blk["start"] - cur[-1]["end"] > cluster_gap:
                clusters.append(cur)
                cur = [blk]
            else:
                cur.append(blk)
        clusters.append(cur)

        for cluster in clusters:
            c_span = cluster[-1]["end"] - cluster[0]["start"]
            if len(cluster) < min_cluster_blocks or c_span < min_copy_span:
                _LOG.debug(
                    "[%s] LILRA6 cluster on %s [%d-%d] skipped: blocks=%d span=%d (spurious).",
                    sample_id, contig, cluster[0]["start"], cluster[-1]["end"],
                    len(cluster), c_span,
                )
                continue

            c_start = cluster[0]["start"]
            c_end   = cluster[-1]["end"]

            hit_block_count: dict[int, int] = defaultdict(int)
            for blk in cluster:
                hit_block_count[id(blk["hit"])] += 1
            best_id = max(hit_block_count, key=hit_block_count.__getitem__)
            rep = next(b["hit"] for b in cluster if id(b["hit"]) == best_id)

            copy_hit = dict(rep)
            copy_hit["start"] = c_start
            copy_hit["end"]   = c_end
            copy_hit["exons"] = [(b["start"], b["end"]) for b in cluster]
            copies.append(copy_hit)

            _LOG.debug(
                "[%s] LILRA6 copy on %s [%d-%d] span=%d bp, %d CDS blocks, identity=%.4f",
                sample_id, contig, c_start, c_end, c_end - c_start,
                len(cluster), rep["identity"],
            )

    _LOG.info(
        "[%s] LILRA6: %d cop%s detected from %d quality hit(s).",
        sample_id, len(copies), "y" if len(copies) == 1 else "ies", len(hits),
    )
    return copies

def _build_locus_contig_map(assignments: list[LocusAssignment]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = defaultdict(list)
    for asgn in assignments:
        if asgn.contig_name not in result[asgn.locus]:
            result[asgn.locus].append(asgn.contig_name)
    return dict(result)

def _log_copy_numbers(sample: HaplotypeSample) -> None:
    present = {k: v for k, v in sample.copy_numbers.items() if v > 0}
    _LOG.info(
        "[%s] Copy numbers: %s",
        sample.sample_hap_id,
        ", ".join(f"{k}={v}" for k, v in sorted(present.items())),
    )

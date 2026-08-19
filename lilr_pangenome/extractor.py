from __future__ import annotations

import logging
from pathlib import Path

from .models import HaplotypeSample
from .resources import LilrResources
from .utils import get_logger, miniprot_to_gff3, run_cmd

_LOG = get_logger(__name__)

_MIN_IDENTITY = 0.75
_MIN_COVERAGE = 0.60

def run_extraction(
    sample: HaplotypeSample,
    resources: LilrResources,
    *,
    threads: int = 4,
    force_rerun: bool = False,
) -> HaplotypeSample:
    out_dir = sample.output_dir / "extractedFasta"
    out_dir.mkdir(parents=True, exist_ok=True)

    gff3_path = out_dir / f"{sample.sample_hap_id}_lilr_miniprot.gff3"
    fai_path  = Path(str(sample.fasta_path) + ".fai")

    if not force_rerun and gff3_path.exists() and gff3_path.stat().st_size > 0 and fai_path.exists():
        _LOG.info("[%s] Skipping extraction - GFF3 and index already exist.", sample.sample_hap_id)
        sample.lilr_contigs_fasta = sample.fasta_path
        return sample

    _LOG.info(
        "[%s] Step 1: Running miniprot against %d LILR protein baits...",
        sample.sample_hap_id, 11,
    )

    try:
        miniprot_to_gff3(
            target=sample.fasta_path,
            proteins=resources.lilr_proteins,
            output_gff3=gff3_path,
            threads=threads,
            max_hits=12,
            min_score_fraction=0.80,
        )
    except Exception as exc:
        sample.failed = True
        sample.failure_message = f"Step 1 miniprot failed: {exc}"
        _LOG.error("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    lilr_contig_names = _get_contig_names_from_gff3(gff3_path)

    if not lilr_contig_names:
        sample.failed = True
        sample.failure_message = "Step 1: No contigs with LILR protein hits found."
        _LOG.warning("[%s] %s", sample.sample_hap_id, sample.failure_message)
        return sample

    _LOG.info(
        "[%s] %d contig(s) with LILR protein hits: %s",
        sample.sample_hap_id, len(lilr_contig_names),
        ", ".join(sorted(lilr_contig_names)),
    )

    if not fai_path.exists():
        _LOG.info("[%s] Indexing assembly with samtools faidx...", sample.sample_hap_id)
        try:
            run_cmd(["samtools", "faidx", str(sample.fasta_path)])
        except Exception as exc:
            sample.failed = True
            sample.failure_message = f"Step 1 samtools faidx failed: {exc}"
            _LOG.error("[%s] %s", sample.sample_hap_id, sample.failure_message)
            return sample

    sample.lilr_contigs_fasta = sample.fasta_path
    return sample

def _get_contig_names_from_gff3(
    gff3_path: Path,
    *,
    min_identity: float = _MIN_IDENTITY,
    min_coverage: float = _MIN_COVERAGE,
) -> set[str]:
    names: set[str] = set()

    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue

            contig = cols[0]
            attrs  = cols[8]

            identity = coverage = 0.0
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("Identity="):
                    try:
                        identity = float(attr[9:])
                    except ValueError:
                        pass
                elif attr.startswith("Positive="):
                    try:
                        coverage = float(attr[9:])
                    except ValueError:
                        pass

            if identity >= min_identity and coverage >= min_coverage:
                names.add(contig)
                _LOG.debug(
                    "  Accepted contig %s (identity=%.3f, coverage=%.3f)",
                    contig, identity, coverage,
                )

    return names

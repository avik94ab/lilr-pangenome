from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd

from .models import LILR_LOCI, HaplotypeSample, SamplePair
from .utils import get_logger

_LOG = get_logger(__name__)

_NOVEL_IDENTITY_THRESHOLD = 0.95

def write_outputs(
    sample_pairs: list[SamplePair],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    _write_copy_number_frame(sample_pairs, output_dir)
    _write_identity_frame(sample_pairs, output_dir)
    _write_sequence_metrics(sample_pairs, output_dir)
    _write_contig_assignments(sample_pairs, output_dir)
    _copy_novel_sequences(sample_pairs, output_dir)
    _write_failure_log(sample_pairs, output_dir)

def _write_copy_number_frame(
    sample_pairs: list[SamplePair], output_dir: Path
) -> None:
    rows = []
    for pair in sample_pairs:
        row = {"sample": pair.name}
        row.update(pair.diploid_copy_numbers())
        rows.append(row)

    df = pd.DataFrame(rows).set_index("sample")[LILR_LOCI]
    df.to_csv(output_dir / "assemblyCopyNumberFrame.csv")
    _LOG.info("Wrote assemblyCopyNumberFrame.csv (%d samples)", len(rows))

def _write_identity_frame(
    sample_pairs: list[SamplePair], output_dir: Path
) -> None:
    rows = []
    for pair in sample_pairs:
        row = {"sample": pair.name}
        row.update({
            locus: round(v, 6) if v == v else None
            for locus, v in pair.diploid_identity().items()
        })
        rows.append(row)

    df = pd.DataFrame(rows).set_index("sample")[LILR_LOCI]
    df.to_csv(output_dir / "assemblyIdentityFrame.csv")
    _LOG.info("Wrote assemblyIdentityFrame.csv")

def _write_sequence_metrics(
    sample_pairs: list[SamplePair], output_dir: Path
) -> None:
    rows = []
    for pair in sample_pairs:
        for hap in pair.haplotypes():
            if hap.failed:
                continue
            for locus, call in hap.allele_calls.items():
                copies = call.best_allele.split("+")
                for copy_label in copies:
                    rows.append({
                        "sample": pair.name,
                        "haplotype": hap.hap_id,
                        "locus": locus,
                        "copy_label": copy_label,
                        "confirmed_locus": call.confirmed_locus,
                        "locus_confirmed": call.locus_confirmed,
                        "identity_to_ref": round(call.identity_to_ref, 6),
                        "mismatches": call.mismatch_score,
                        "novel_allele": call.novel_allele,
                        "sequence_fasta": str(call.sequence_fasta) if call.sequence_fasta else "",
                    })

    df = pd.DataFrame(rows) if rows else pd.DataFrame()
    df.to_csv(output_dir / "assemblySequenceMetrics.csv", index=False)
    _LOG.info("Wrote assemblySequenceMetrics.csv (%d rows)", len(rows))

def _write_contig_assignments(
    sample_pairs: list[SamplePair], output_dir: Path
) -> None:
    rows = []
    for pair in sample_pairs:
        for hap in pair.haplotypes():
            for asgn in hap.locus_assignments:
                rows.append({
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

    df = pd.DataFrame(rows) if rows else pd.DataFrame()
    df.to_csv(output_dir / "assemblyContigAssignments.csv", index=False)
    _LOG.info("Wrote assemblyContigAssignments.csv")

def _copy_novel_sequences(
    sample_pairs: list[SamplePair], output_dir: Path
) -> None:
    novel_dir = output_dir / "novel_sequences"
    n_copied = 0

    for pair in sample_pairs:
        for hap in pair.haplotypes():
            if hap.failed:
                continue
            for call in hap.allele_calls.values():
                if call.sequence_fasta is None:
                    continue
                if call.identity_to_ref >= _NOVEL_IDENTITY_THRESHOLD:
                    continue
                dest_dir = novel_dir / pair.name
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest = dest_dir / call.sequence_fasta.name
                shutil.copy2(call.sequence_fasta, dest)
                n_copied += 1

    if n_copied:
        _LOG.info(
            "Copied %d novel sequence(s) (identity < %.0f%%) to %s/",
            n_copied, _NOVEL_IDENTITY_THRESHOLD * 100, novel_dir,
        )

def _write_failure_log(
    sample_pairs: list[SamplePair], output_dir: Path
) -> None:
    lines: list[str] = []
    for pair in sample_pairs:
        for hap in pair.haplotypes():
            if hap.failed:
                lines.append(f"{hap.sample_hap_id}\t{hap.failure_message}")

    log_path = output_dir / "failure.log"
    log_path.write_text("\n".join(lines) + ("\n" if lines else ""))

    if lines:
        _LOG.warning("%d haplotype(s) failed - see %s", len(lines), log_path)
    else:
        _LOG.info("No failures.")

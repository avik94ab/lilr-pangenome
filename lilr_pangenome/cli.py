from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .copy_number import run_copy_number
from .extractor import run_extraction
from .finalize import write_outputs
from .genotyper import run_genotyping
from .models import HaplotypeSample, SamplePair
from .resources import load_resources
from .utils import check_tools, configure_logging, get_logger

_LOG = get_logger(__name__)

def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    configure_logging(args.log_level)

    _LOG.info("lilr-pangenome starting.")

    try:
        check_tools()
    except RuntimeError as exc:
        _LOG.error("%s", exc)
        return 1

    try:
        resources = load_resources()
    except FileNotFoundError as exc:
        _LOG.error("Resource validation failed:\n%s", exc)
        return 1

    assembly_dir = Path(args.assembly_dir)
    output_dir   = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_pairs = _discover_samples(
        assembly_dir,
        hap1_suffix=args.hap1_suffix,
        hap2_suffix=args.hap2_suffix,
        sample_filter=set(args.samples) if args.samples else None,
        output_dir=output_dir,
    )

    if not sample_pairs:
        _LOG.error(
            "No sample pairs found in %s with suffixes '%s' / '%s'.",
            assembly_dir, args.hap1_suffix, args.hap2_suffix,
        )
        return 1

    _LOG.info("Processing %d sample(s).", len(sample_pairs))

    for pair in sample_pairs:
        _LOG.info("=== Sample: %s ===", pair.name)
        for hap in pair.haplotypes():
            hap = run_extraction(
                hap, resources,
                threads=args.threads,
                force_rerun=args.force_rerun,
            )
            hap = run_copy_number(
                hap, resources,
                threads=args.threads,
                min_identity=args.min_identity,
                min_query_cov=args.min_query_cov,
                force_rerun=args.force_rerun,
            )
            hap = run_genotyping(
                hap, resources,
                threads=args.threads,
                force_rerun=args.force_rerun,
            )

    _LOG.info("Writing outputs...")
    write_outputs(sample_pairs, output_dir)
    _LOG.info("Done. Results in %s", output_dir)
    return 0

def _discover_samples(
    assembly_dir: Path,
    *,
    hap1_suffix: str,
    hap2_suffix: str,
    sample_filter: set[str] | None,
    output_dir: Path,
) -> list[SamplePair]:
    hap1_files = {
        f.name[: -len(hap1_suffix)]: f
        for f in assembly_dir.iterdir()
        if f.name.endswith(hap1_suffix)
    }
    hap2_files = {
        f.name[: -len(hap2_suffix)]: f
        for f in assembly_dir.iterdir()
        if f.name.endswith(hap2_suffix)
    }

    common = sorted(set(hap1_files) & set(hap2_files))
    if sample_filter:
        common = [s for s in common if s in sample_filter]

    pairs = []
    for name in common:
        sample_out = output_dir / name
        hap1 = HaplotypeSample(
            name=name, hap_id="hap1",
            fasta_path=hap1_files[name],
            output_dir=sample_out / "hap1",
        )
        hap2 = HaplotypeSample(
            name=name, hap_id="hap2",
            fasta_path=hap2_files[name],
            output_dir=sample_out / "hap2",
        )
        pairs.append(SamplePair(name=name, hap1=hap1, hap2=hap2))

    return pairs

def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="lilr-pangenome",
        description="LILR copy number and sequence genotyping from haplotype-resolved assemblies.",
    )
    p.add_argument("--assembly-dir", required=True,
                   help="Directory containing paired haplotype FASTA files.")
    p.add_argument("--output-dir", required=True,
                   help="Directory for output files (created if absent).")
    p.add_argument("--hap1-suffix", default="_hap1.fa.gz",
                   help="Filename suffix identifying hap1 files. [%(default)s]")
    p.add_argument("--hap2-suffix", default="_hap2.fa.gz",
                   help="Filename suffix identifying hap2 files. [%(default)s]")
    p.add_argument("--samples", nargs="*", metavar="SAMPLE",
                   help="Process only these sample IDs (default: all).")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for minimap2/miniprot. [%(default)s]")
    p.add_argument("--min-identity", type=float, default=0.75,
                   help="Minimum protein identity for locus assignment. [%(default)s]")
    p.add_argument("--min-query-cov", type=float, default=0.65,
                   help="Minimum protein coverage for locus assignment. [%(default)s]")
    p.add_argument("--force-rerun", action="store_true",
                   help="Ignore cached intermediate files and rerun all steps.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                   help="Logging verbosity. [%(default)s]")
    return p.parse_args(argv)

if __name__ == "__main__":
    sys.exit(main())

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)

def configure_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"Required external tool '{name}' not found on PATH. "
            "Please install it and ensure it is accessible."
        )
    return path

def check_tools() -> None:
    require_tool("minimap2")
    require_tool("samtools")
    require_tool("miniprot")

def run_cmd(cmd: list[str], *, check: bool = True, capture: bool = False) -> subprocess.CompletedProcess:
    logger = get_logger(__name__)
    logger.debug("Running: %s", " ".join(str(c) for c in cmd))

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logger.debug("[%s stderr] %s", cmd[0], line)
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, result.stdout, result.stderr
        )
    return result

def run_piped(cmds: list[list[str]], output_path: Path) -> None:
    logger = get_logger(__name__)
    logger.debug(
        "Running pipeline: %s",
        " | ".join(" ".join(str(c) for c in cmd) for cmd in cmds),
    )

    procs: list[subprocess.Popen] = []
    prev_stdout = None
    for i, cmd in enumerate(cmds):
        is_last = i == len(cmds) - 1
        stdout = open(output_path, "wb") if is_last else subprocess.PIPE
        proc = subprocess.Popen(
            cmd,
            stdin=prev_stdout,
            stdout=stdout,
            stderr=subprocess.PIPE,
        )
        if prev_stdout is not None:
            prev_stdout.close()
        prev_stdout = proc.stdout
        procs.append(proc)

    for proc in procs:
        proc.wait()

    for proc in procs:
        if proc.returncode != 0:
            stderr = proc.stderr.read().decode() if proc.stderr else ""
            raise subprocess.CalledProcessError(proc.returncode, proc.args, stderr=stderr)

def extract_fasta_region(
    fasta_path: Path,
    contig: str,
    start: int,
    end: int,
    *,
    strand: str = "+",
):
    import io
    from Bio import SeqIO as _SeqIO

    region = f"{contig}:{start + 1}-{end}"
    result = run_cmd(["samtools", "faidx", str(fasta_path), region], capture=True)

    records = list(_SeqIO.parse(io.StringIO(result.stdout), "fasta"))
    if not records:
        raise ValueError(f"samtools faidx returned no sequence for {region}")

    seq = records[0].seq
    return seq.reverse_complement() if strand == "-" else seq

def minimap2_to_bam(
    reference: Path,
    query: Path,
    output_bam: Path,
    *,
    preset: str = "asm5",
    threads: int = 4,
    extra_args: list[str] | None = None,
) -> None:
    mm2_cmd = [
        "minimap2",
        "-ax", preset,
        "--cs",
        "--secondary=no",
        "-t", str(threads),
    ]
    if extra_args:
        mm2_cmd.extend(extra_args)
    mm2_cmd += [str(reference), str(query)]

    filter_cmd = ["samtools", "view", "-F", "4", "-u", "-"]

    sort_cmd = ["samtools", "sort", "-@", str(threads)]

    run_piped([mm2_cmd, filter_cmd, sort_cmd], output_bam)

    run_cmd(["samtools", "index", str(output_bam)])

def miniprot_to_gff3(
    target: Path,
    proteins: Path,
    output_gff3: Path,
    *,
    threads: int = 4,
    max_intron: int = 200_000,
    max_hits: int = 6,
    min_score_fraction: float = 0.80,
) -> None:
    cmd = [
        "miniprot",
        "--gff",
        "-G", str(max_intron),
        f"--outn={max_hits}",
        f"--outs={min_score_fraction:.2f}",
        "-t", str(threads),
        str(target),
        str(proteins),
    ]
    result = run_cmd(cmd, capture=True)
    output_gff3.write_text(result.stdout)

def minimap2_to_paf(
    reference: Path,
    query: Path,
    output_paf: Path,
    *,
    preset: str = "asm5",
    threads: int = 4,
    extra_args: list[str] | None = None,
) -> None:
    cmd = [
        "minimap2",
        "-x", preset,
        "--cs",
        "-c",
        "-t", str(threads),
    ]
    if extra_args:
        cmd.extend(extra_args)
    cmd += [str(reference), str(query)]

    result = run_cmd(cmd, capture=True)
    output_paf.write_text(result.stdout)

def _open_fasta(path: Path):
    if path.suffix in (".gz", ".bgz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def read_fasta(path: Path) -> dict[str, SeqRecord]:
    with _open_fasta(path) as fh:
        return SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

def iter_fasta(path: Path) -> Iterator[SeqRecord]:
    with _open_fasta(path) as fh:
        yield from SeqIO.parse(fh, "fasta")

def write_fasta(records: list[SeqRecord], output_path: Path, *, compress: bool = False) -> None:
    if compress:
        with gzip.open(output_path, "wt") as fh:
            SeqIO.write(records, fh, "fasta")
    else:
        with open(output_path, "w") as fh:
            SeqIO.write(records, fh, "fasta")

def filter_fasta_by_names(
    input_fasta: Path,
    names: set[str],
    output_fasta: Path,
    *,
    compress: bool = True,
) -> int:
    kept: list[SeqRecord] = [
        rec for rec in iter_fasta(input_fasta) if rec.id in names
    ]
    write_fasta(kept, output_fasta, compress=compress)
    return len(kept)

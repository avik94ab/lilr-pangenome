from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

LILR_LOCI: list[str] = [
    "LILRA1",
    "LILRA2",
    "LILRA3",
    "LILRA4",
    "LILRA5",
    "LILRA6",
    "LILRB1",
    "LILRB2",
    "LILRB3",
    "LILRB4",
    "LILRB5",
]

@dataclass
class AlleleCall:

    locus: str
    best_allele: str
    mismatch_score: int
    identity_to_ref: float
    novel_allele: bool
    sequence_fasta: Optional[Path] = None
    cdna_fasta:     Optional[Path] = None
    protein_fasta:  Optional[Path] = None

    confirmed_locus: str = ""
    locus_confirmed: bool = True

    ambiguous_candidates: list[str] = field(default_factory=list)

    def display(self) -> str:
        if self.ambiguous_candidates:
            return " ".join(sorted({self.best_allele} | set(self.ambiguous_candidates)))
        return self.best_allele

@dataclass
class LocusAssignment:

    contig_name: str
    locus: str
    identity: float
    query_coverage: float
    ref_allele: str
    query_start: int = 0
    query_end: int = 0
    strand: str = "+"
    cds_coords: list = field(default_factory=list)

@dataclass
class HaplotypeSample:

    name: str
    hap_id: str
    fasta_path: Path
    output_dir: Path

    failed: bool = False
    failure_message: str = ""

    lilr_contigs_fasta: Optional[Path] = None

    locus_contig_map: dict[str, list[str]] = field(default_factory=dict)
    copy_numbers: dict[str, int] = field(default_factory=dict)
    locus_assignments: list[LocusAssignment] = field(default_factory=list)
    allele_calls: dict[str, AlleleCall] = field(default_factory=dict)

    @property
    def sample_hap_id(self) -> str:
        return f"{self.name}_{self.hap_id}"

@dataclass
class SamplePair:

    name: str
    hap1: HaplotypeSample
    hap2: HaplotypeSample

    def haplotypes(self) -> list[HaplotypeSample]:
        return [self.hap1, self.hap2]

    def diploid_copy_numbers(self) -> dict[str, int]:
        result: dict[str, int] = {}
        for locus in LILR_LOCI:
            result[locus] = (
                self.hap1.copy_numbers.get(locus, 0)
                + self.hap2.copy_numbers.get(locus, 0)
            )
        return result

    def diploid_identity(self) -> dict[str, float]:
        result: dict[str, float] = {}
        for locus in LILR_LOCI:
            values: list[float] = []
            for hap in self.haplotypes():
                call = hap.allele_calls.get(locus)
                if call is not None:
                    values.append(call.identity_to_ref)
            result[locus] = sum(values) / len(values) if values else float("nan")
        return result

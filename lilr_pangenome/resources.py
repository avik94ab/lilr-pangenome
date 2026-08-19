from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .models import LILR_LOCI

_RESOURCES_DIR = Path(__file__).parent.parent / "resources"

@dataclass(frozen=True)
class LilrResources:

    @property
    def lilr_proteins(self) -> Path:
        p = _RESOURCES_DIR / "lilr_proteins.faa"
        if not p.exists():
            raise FileNotFoundError(f"LILR protein bait not found: {p}")
        return p

    @property
    def hg38_refs_dir(self) -> Path:
        return _RESOURCES_DIR / "lilr_hg38_refs"

    def locus_reference(self, locus: str) -> Path:
        p = self.hg38_refs_dir / f"{locus}.fasta"
        if not p.exists():
            raise FileNotFoundError(f"No hg38 reference for locus {locus}: {p}")
        return p

    @property
    def combined_locus_reference(self) -> Path:
        p = self.hg38_refs_dir / "all_loci.fasta"
        if not p.exists():
            raise FileNotFoundError(
                f"Combined locus reference not found: {p}\n"
                "Recreate with: cat resources/lilr_hg38_refs/*.fasta > "
                "resources/lilr_hg38_refs/all_loci.fasta"
            )
        return p

def load_resources() -> LilrResources:
    res = LilrResources()
    _validate(res)
    return res

def _validate(res: LilrResources) -> None:
    missing: list[str] = []
    try:
        _ = res.lilr_proteins
    except FileNotFoundError as exc:
        missing.append(str(exc))
    try:
        _ = res.combined_locus_reference
    except FileNotFoundError as exc:
        missing.append(str(exc))
    for locus in LILR_LOCI:
        try:
            _ = res.locus_reference(locus)
        except FileNotFoundError as exc:
            missing.append(str(exc))
    if missing:
        raise FileNotFoundError(
            "Missing required resource files:\n  " + "\n  ".join(missing)
        )

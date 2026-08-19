# lilr-pangenome

Copy-number and sequence genotyping of the eleven human LILR genes
(*LILRA1-6*, *LILRB1-5*; chr19q13.4, leukocyte receptor complex) from
haplotype-resolved pangenome assemblies.

Given a pair of haplotype FASTA files per sample, the tool extracts the
LILR-containing contigs, assigns each hit to one of the 11 loci, counts
per-gene copy number, and scores each copy's sequence identity to a bundled
per-locus reference (flagging novel sequences). All reference data required to
run is bundled in [`resources/`](resources/); no external database download is
needed.

## How it works

The pipeline runs four steps per haplotype:

1. **Extract** LILR-containing contigs with `miniprot`, using the 11 canonical
   LILR proteins ([`resources/lilr_proteins.faa`](resources/lilr_proteins.faa)) as bait.
2. **Copy number**: parse the `miniprot` GFF3, assign each hit to a locus, and
   count copies per gene.
3. **Genotype**: extract each per-locus copy and align it to the bundled hg38
   reference with `minimap2` to score identity and flag novel sequences.
4. **Finalize**: write cohort-level CSVs and novel-sequence FASTAs.

Ten of the eleven locus references derive from GRCh38 (NC_000019.10). *LILRA3*
is absent from GRCh38, so its reference comes from an alternate assembly;
identity scores for *LILRA3* are relative to that alternate.

## Requirements

- Python >= 3.10 with `biopython` >= 1.81 and `pandas` >= 2.0 (installed below).
- Three external tools on `PATH`: [`miniprot`](https://github.com/lh3/miniprot),
  [`minimap2`](https://github.com/lh3/minimap2), and
  [`samtools`](https://github.com/samtools/samtools).

## Install

The bundled reference data in `resources/` is resolved relative to the package
location, so install in editable mode from the cloned repository (this keeps the
package next to `resources/`):

```bash
git clone https://github.com/avik94ab/lilr-pangenome.git
cd lilr-pangenome
pip install -e .
```

Equivalently, without installing, run the module directly from the repo root:
`python -m lilr_pangenome.cli ...`.

## Usage

Input is a directory of paired haplotype assemblies named
`{SAMPLE}{hap1_suffix}` and `{SAMPLE}{hap2_suffix}`
(default suffixes `_hap1.fa.gz` / `_hap2.fa.gz`):

```
assemblies/
  NA12878_hap1.fa.gz
  NA12878_hap2.fa.gz
```

Run:

```bash
lilr-pangenome --assembly-dir assemblies/ --output-dir results/
```

Common options:

| Option | Default | Meaning |
|---|---|---|
| `--hap1-suffix` / `--hap2-suffix` | `_hap1.fa.gz` / `_hap2.fa.gz` | Haplotype filename suffixes |
| `--samples S1 S2 ...` | all | Restrict to specific sample IDs |
| `--threads` | 4 | Threads for `minimap2` / `miniprot` |
| `--min-identity` | 0.75 | Minimum protein identity for locus assignment |
| `--min-query-cov` | 0.65 | Minimum protein coverage for locus assignment |
| `--force-rerun` | off | Ignore cached intermediates and rerun all steps |

Run `lilr-pangenome --help` for the full list.

## Smoke test

The bundled combined reference doubles as a self-contained functional test
input (no external assembly needed). It contains the 11 locus references, so
using it for both haplotypes should call all 11 loci at diploid copy number 2
(one copy per haplotype). Input must be `bgzip`-compressed, as `samtools faidx`
requires, which is how real assembly FASTAs are distributed:

```bash
mkdir -p smoke/in
bgzip -c resources/lilr_hg38_refs/all_loci.fasta > smoke/in/TEST_hap1.fa.gz
bgzip -c resources/lilr_hg38_refs/all_loci.fasta > smoke/in/TEST_hap2.fa.gz
lilr-pangenome --assembly-dir smoke/in --output-dir smoke/out --threads 4
cat smoke/out/assemblyCopyNumberFrame.csv
# sample,LILRA1,LILRA2,...,LILRB5
# TEST,2,2,2,2,2,2,2,2,2,2,2
```

## Repository layout

```
lilr_pangenome/        the genotyping pipeline (Python package)
  cli.py               command-line entry point
  extractor.py         step 1: pull LILR contigs from each assembly (miniprot)
  copy_number.py       step 2: locus assignment and per-gene copy number
  genotyper.py         step 3: per-copy sequence extraction and identity scoring
  finalize.py          step 4: cohort-level output writers
  resources.py         reference-path resolution and validation
  models.py            pipeline data classes
  utils.py             subprocess, FASTA I/O, logging helpers
resources/             bundled reference data (the "database")
  lilr_proteins.faa            11 canonical LILR protein baits
  lilr_hg38_refs/*.fasta       per-locus references + combined all_loci.fasta
  lilr_uniprot_domains.json    UniProt domain annotations per locus
pyproject.toml         pip-installable package definition
```

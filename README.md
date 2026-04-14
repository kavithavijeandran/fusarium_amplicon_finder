# fusarium-amplicon-finder

**In silico amplicon screening tool for evaluating PCR primer specificity across Fusarium genomes**

## Overview

This repository contains the custom Python scripts developed for the in silico evaluation of published PCR primers for *Fusarium* detection in wheat, as described in:

> Vijeandran K, Larionov A, Cervini C, Campbell K, Walkowiak S, Pasquali M, Richard-Forget F, Brown N, Rose LJ, Verheecke-Vaessen C. *Revisiting the PCR-based molecular approaches for mycotoxigenic Fusarium species detection and quantification in wheat.* [Journal name] [Year]. DOI: [to be added upon publication]

The tool performs genome-wide in silico PCR screening — identifying primer binding sites across genome assemblies (FASTA format), evaluating all possible primer orientations including reverse complements, and detecting putative amplicons within a user-defined size range. It supports degenerate IUPAC primer codes and flexible mismatch tolerance, making it suitable for evaluating legacy primers against modern genome assemblies.

---

## Features

- Screens multiple genome FASTA files in a single run
- Supports all IUPAC degenerate base codes (R, Y, S, W, K, M, B, D, H, V, N)
- Configurable mismatch tolerance (0–5 mismatches per primer)
- Evaluates all four primer orientation combinations (plus/minus strand, forward/reverse)
- Filters amplicons by size range
- Memory-efficient deduplication of identical amplicons across contigs
- Outputs structured CSV with amplicon coordinates, sequences, mismatch counts, orientation, and panel metadata

---

## Repository structure

```
fusarium-amplicon-finder/
├── amplicon_finder.py          # Main Python script
├── requirements.txt            # Python dependencies
├── LICENSE
└── README.md
```

---

## Requirements

- Python 3.8 or higher
- [Biopython](https://biopython.org/)

Install dependencies with:

```bash
pip install -r requirements.txt
```

---

## Usage

```bash
python amplicon_finder.py \
  -d /path/to/genomes/ \
  -f FORWARD_PRIMER_SEQUENCE \
  -r REVERSE_PRIMER_SEQUENCE \
  --forward_name PRIMER_NAME_F \
  --reverse_name PRIMER_NAME_R \
  --min 50 \
  --max 300 \
  --mismatches 5 \
  --outdir /path/to/output/ \
  --progress
```

### Arguments

| Argument | Description | Default |
|---|---|---|
| `-d`, `--directory` | Directory containing genome FASTA files (searches `core/` and `extended/` subdirectories) | Required |
| `-f`, `--forward` | Forward primer sequence (IUPAC codes supported) | Required |
| `-r`, `--reverse` | Reverse primer sequence (IUPAC codes supported) | Required |
| `--forward_name` | Forward primer label (e.g. ITS1) | — |
| `--reverse_name` | Reverse primer label (e.g. ITS4) | — |
| `--min` | Minimum amplicon size in bp | 100 |
| `--max` | Maximum amplicon size in bp | 5000 |
| `--mismatches` | Maximum mismatches allowed per primer | 5 |
| `--no_dedup` | Disable deduplication of identical amplicons | Off |
| `--outdir` | Output directory | `./amplicon_results` |
| `--progress` | Print progress for each genome processed | Off |

### Example

```bash
python amplicon_finder.py \
  -d genomes/ \
  -f TGGCCTGAATGAAGGATTTCTAG \
  -r CATCGTTGTTAACTTATTGGAGATG \
  --forward_name COB1 \
  --reverse_name COB2 \
  --min 50 \
  --max 300 \
  --mismatches 5 \
  --outdir results/ \
  --progress
```

---

## Genome directory structure

The script expects genomes to be organised into `core/` and `extended/` subdirectories, which are used to assign a panel label to each result:

```
genomes/
├── core/
│   ├── GCA_000240135.fna
│   └── ...
└── extended/
    ├── GCA_003012315.fna
    └── ...
```

Accepted file extensions: `.fna`, `.fasta`

---

## Output

Results are saved to a CSV file at `<outdir>/amplicon_finder/amplicon_results_<FORWARD_NAME>-<REVERSE_NAME>.csv`

The CSV includes a metadata header block followed by one row per detected amplicon:

| Column | Description |
|---|---|
| Species Name | Extracted from the FASTA header |
| Accession | Derived from the genome filename |
| Contig | Contig/chromosome ID |
| Contig Length | Total length of the contig (bp) |
| Start | Amplicon start position (1-based) |
| End | Amplicon end position |
| Amplicon Sequence | Full amplicon sequence |
| Amplicon Size | Amplicon length in bp |
| Forward Primer Match Site | Genomic sequence at forward primer binding site |
| Reverse Primer Match Site | Genomic sequence at reverse primer binding site |
| Orientation | Strand orientation of the amplicon |
| Forward Primer Mismatches | Number of mismatches at forward primer site |
| Reverse Primer Mismatches | Number of mismatches at reverse primer site |
| Mismatch Tier | Maximum of forward/reverse mismatch counts |
| Panel | `core` or `extended` based on genome directory |

---

## HPC usage

An example PBS job script (`amplicon_finder.pbs`) is provided for running the tool on an HPC cluster with a PBS/Torque scheduler. Edit the primer sequences, size range, and directory paths at the top of the script before submitting.

```bash
qsub amplicon_finder.pbs
```

---

## Citation

If you use this tool, please cite the associated publication:

> Vijeandran K, et al. *Revisiting the PCR-based molecular approaches for mycotoxigenic Fusarium species detection and quantification in wheat.* [Journal] [Year]. DOI: [to be added]

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Kavitha Vijeandran — Magan Centre of Applied Mycology, Cranfield University
Kavitha.Vijeandran@cranfield.ac.uk

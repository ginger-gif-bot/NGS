# NGS Short-Read Processing Pipeline

A reproducible, end-to-end Bash and Python pipeline for processing paired-end Illumina short-read sequencing data — from raw SRA accessions through to extracted coding sequences (CDS) ready for downstream analysis.

Built as part of an M.Sc. Bioinformatics thesis project. Currently applied to a comparative codon usage bias study of drug-sensitive versus multidrug-resistant *Mycobacterium tuberculosis* isolates using the CRyPTIC consortium dataset.

---

## Pipeline Overview

```
SRA Accessions
      ↓
Raw FASTQ (prefetch + fasterq-dump)
      ↓
Quality Control (FastQC + MultiQC)
      ↓
Adapter Trimming (Trimmomatic PE)
      ↓
Post-trim QC (FastQC + MultiQC)
      ↓
Reference Indexing (BWA index)
      ↓
Alignment (BWA MEM)
      ↓
Sorting + Indexing (samtools sort)
      ↓
Variant Calling (bcftools mpileup + call)
      ↓
Consensus Sequence (bcftools consensus)
      ↓
CDS Extraction (Biopython + pyliftover)
```

---

## Repository Structure

```
NGS_pipeline/
├── scripts/
│   ├── config.sh               # Shared configuration — paths, threads, reference
│   ├── getting_data.sh         # SRA download and FASTQ conversion
│   ├── fastqc_raw.sh           # Quality control on raw reads
│   ├── trimming.sh             # Adapter and quality trimming
│   ├── fastqc_trim.sh          # Quality control on trimmed reads
│   ├── multiqc_raw.sh          # Aggregate raw QC report
│   ├── multiqc_trim.sh         # Aggregate trimmed QC report
│   ├── reference.sh            # Reference genome indexing
│   ├── bwa_alignment.sh        # Read alignment to reference
│   ├── samtools_sort.sh        # SAM to sorted, indexed BAM
│   ├── variant_calling.sh      # Variant calling per sample
│   ├── consensus.sh            # Consensus FASTA + chain file generation
│   └── cleanup_sams.sh         # Safe SAM deletion after BAM verification
│
├── python_scripts/
│   └── cds_extraction.py       # CDS extraction using annotation + coordinate liftover
│
├── metadata/                   # Sample accession lists and CRyPTIC metadata
├── data/
│   ├── raw_reads/              # Raw SRA and FASTQ files
│   └── reference/              # Reference genome and annotation files
├── results/                    # All pipeline outputs (gitignored)
├── logs/                       # Per-stage and per-sample log files
└── .gitignore                  # Excludes large data files from version control
```

---

## Key Design Principles

### Idempotency
Every stage checks for existing output before processing. Re-running any script safely skips completed samples and resumes from where it left off. An interrupted run at sample 47 of 100 resumes at sample 48 — no reprocessing of completed work.

### Atomic Outputs
All tools that write output use a `.tmp` staging pattern: output is written to a temporary file and renamed to its final name only after the tool exits successfully. An interrupted or failed run never leaves a file that looks complete but isn't.

### Failure Handling
Every script uses `set -e` (and `set -o pipefail` where pipes are used), so failures halt execution immediately rather than silently producing downstream garbage. Cleanup guards (`||  { rm -f ...; exit 1; }`) remove partial outputs on failure, preventing future runs from treating them as complete.

### Per-Sample Logging
Alignment logs are written per sample to `logs/bwa/`, allowing per-sample mapping rates and errors to be inspected independently. Stage-level logs capture full stdout and stderr for every other stage.

### Reproducibility
- All shared paths and parameters live in `scripts/config.sh` — one file to update when adapting to new data or a new reference genome
- Reference genome provenance is documented in `scripts/reference.sh` with exact NCBI accession and download steps
- Sample selection is scripted from metadata (no manual curation)
- All scripts are version-controlled; all large data files are gitignored

### Timing
Every processing stage reports elapsed time per sample, making it straightforward to estimate runtime for new datasets and identify unexpectedly slow samples.

---

## Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| SRA Toolkit (prefetch, fasterq-dump) | 3.x | SRA download and FASTQ conversion |
| FastQC | 0.12.x | Per-sample read quality assessment |
| MultiQC | 1.x | Aggregate quality reports |
| Trimmomatic | 0.40 | Adapter trimming and quality filtering |
| BWA | 0.7.19 | Short-read alignment |
| samtools | 1.x | SAM/BAM processing and indexing |
| bcftools | 1.24 | Variant calling and consensus generation |
| Biopython | 1.8x | GBFF parsing and sequence extraction |
| pyliftover | 0.4 | Coordinate liftover using chain files |

---

## Configuration

All shared variables are defined in `scripts/config.sh`:

```bash
# Reference genome path
REFERENCE="data/reference/.../genomic.fna"

# Input/output directories
TRIMMED_READS="results/trimmed"
SAM_FILES="results/alignment"
SORTED_BAM="results/sorted"
VCF_DIR="results/vcf"
CONSENSUS="results/consensus"
CHAIN="results/chain"

# Compute resources
THREADS=8
```

Adapting the pipeline to a new organism or dataset requires updating this file and the reference genome — all downstream scripts source it automatically.

---

## CDS Extraction

The final stage extracts target coding sequences from each consensus genome using:

1. **GBFF annotation parsing** — Biopython reads the reference annotation to locate CDS features by gene name, capturing coordinates and strand
2. **Coordinate liftover** — pyliftover uses the chain files generated during consensus building to translate reference coordinates to consensus coordinates, correctly accounting for insertions and deletions
3. **Sequence extraction** — the CDS region is sliced from the consensus sequence and reverse-complemented if on the minus strand
4. **Output** — one FASTA file per sample per gene, organized into gene-named subdirectories

The `gene_name_extraction()` helper function allows searching the annotation by substring, making it straightforward to identify correct gene names for any new target before running the full extraction.

---

## Current Application

This pipeline is being applied to 100 *Mycobacterium tuberculosis* isolates from the CRyPTIC consortium dataset (PLOS Biology, 2022) — 50 drug-sensitive and 50 multidrug-resistant — for a comparative codon usage bias analysis of resistance-associated genes (*rpoB*, *katG*, *gyrA*).

Reference genome: H37Rv (NC_000962.3)

---

## Environment Setup

Due to conda dependency conflicts encountered during development, this pipeline uses four separate environments. If you are lucky, fewer environments may work — attempt to install everything into one environment first, and only create separate ones if conflicts arise.

**Attempt a single environment first:**

```bash
conda create -n ngs_env -c bioconda -c conda-forge \
    sra-tools fastqc trimmomatic multiqc \
    bwa samtools bcftools biopython
pip install pyliftover
```

If conflicts arise (which they likely will), use the four-environment setup below.

**ngs_env** — alignment, sorting, SRA download, CDS extraction:

```bash
conda create -n ngs_env -c bioconda -c conda-forge \
    sra-tools bwa samtools biopython
pip install pyliftover
```

**qc_env** — FastQC and MultiQC:

```bash
conda create -n qc_env -c bioconda -c conda-forge \
    fastqc multiqc
```

**trim_env** — Trimmomatic:

```bash
conda create -n trim_env -c bioconda -c conda-forge \
    trimmomatic
```

**bcftools_env** — variant calling and consensus generation:

```bash
conda create -n bcftools_env -c bioconda -c conda-forge \
    bcftools
```

**Which environment to activate per stage:**

| Stage | Script | Environment |
|-------|--------|-------------|
| Download + FASTQ conversion | getting_data.sh | ngs_env |
| Raw QC | fastqc_raw.sh | qc_env |
| Trimming | trimming.sh | trim_env |
| Trimmed QC | fastqc_trim.sh, multiqc_*.sh | qc_env |
| Reference indexing + Alignment | reference.sh, bwa_alignment.sh | ngs_env |
| Sorting | samtools_sort.sh | ngs_env |
| Variant calling + Consensus | variant_calling.sh, consensus.sh | bcftools_env |
| CDS extraction | cds_extraction.py | ngs_env |

All scripts were developed and tested on Ubuntu 24.04 (WSL2 on Windows 11) with Python 3.13.

---

## Running the Pipeline

```bash
# From the project root directory

conda activate ngs_env

bash scripts/getting_data.sh      # Download + convert SRA
bash scripts/fastqc_raw.sh        # Raw QC
bash scripts/trimming.sh          # Trim reads
bash scripts/fastqc_trim.sh       # Post-trim QC
bash scripts/multiqc_raw.sh       # Aggregate raw QC
bash scripts/multiqc_trim.sh      # Aggregate trimmed QC
bash scripts/reference.sh         # Index reference
bash scripts/bwa_alignment.sh     # Align reads
bash scripts/samtools_sort.sh     # Sort and index BAMs
bash scripts/variant_calling.sh   # Call variants
bash scripts/consensus.sh         # Generate consensus + chain files
python python_scripts/cds_extraction.py  # Extract CDS
```

Or chain stages with automatic handoff:

```bash
bash scripts/getting_data.sh && \
bash scripts/fastqc_raw.sh && \
bash scripts/trimming.sh && \
bash scripts/fastqc_trim.sh
```

---

## Notes

- All paths are relative to the project root — always run scripts from `NGS_pipeline/`
- The pipeline was developed and tested on Ubuntu 24 (WSL2 on Windows 11)
- Large data files (FASTQ, BAM, VCF, SRA) are excluded from version control via `.gitignore`; only scripts, metadata, and empty directory placeholders (`.gitkeep`) are tracked
- SAM files are deleted after BAM verification using `cleanup_sams.sh` to manage disk usage

---

## Author

M.Sc. Biotechnology — Integrated Programme  
Pipeline developed as part of thesis preparation, 2026

# NGS Pipeline — Influenza A H1N1 Codon Usage Bias Study

A reproducible bash pipeline for processing Influenza A (H1N1) raw sequencing reads from SRA through quality control and trimming, in preparation for reference-based alignment, variant calling, and downstream codon usage bias analysis comparing avian (duck) and human host adaptation.

---

## Project Structure

```
NGS_pipeline/
├── data/
│   └── raw_reads/
│       ├── raw_sra/          # Downloaded .sra files from NCBI SRA (prefetch output)
│       └── raw_fastq/        # Converted paired-end FASTQ files (fasterq-dump output)
├── metadata/
│   ├── SraRunTable.csv       # Full metadata exported from SRA Run Selector
│   └── selected_sra_run_ids  # Filtered accession list (H1N1, duck + human)
├── results/
│   ├── fastqc_raw/           # FastQC HTML + ZIP reports for raw reads
│   ├── multiqc_raw/          # MultiQC combined report for raw reads
│   ├── trimmed/              # Trimmomatic output — trimmed paired/unpaired FASTQ files
│   ├── fastqc_trimmed/       # FastQC HTML + ZIP reports for trimmed reads
│   ├── multiqc_trimmed/      # MultiQC combined report for trimmed reads
│   ├── alignment/            # BWA-MEM SAM alignment files (next stage)
│   ├── bam/                  # Sorted and indexed BAM files (next stage)
│   └── failed_samples.txt    # Samples that failed trimming — excluded from analysis
├── logs/                     # Pipeline run logs
├── scripts/                  # All pipeline bash scripts (see below)
└── .gitignore                # Excludes large data and result files from GitHub
```

---

## Scripts

Each stage of the pipeline is a separate script. Run them in order from the project root.

| Script | Stage | Tool | Input | Output |
|--------|-------|------|-------|--------|
| `getting_data.sh` | Data acquisition | prefetch, fasterq-dump | SraRunTable.csv + accession list | raw_sra/, raw_fastq/ |
| `fastqc_raw.sh` | Quality check (raw) | FastQC | raw_fastq/*.fastq | fastqc_raw/ |
| `multiqc_raw.sh` | Aggregate QC (raw) | MultiQC | fastqc_raw/ | multiqc_raw/ |
| `trimming.sh` | Adapter trimming | Trimmomatic PE | raw_fastq/ | trimmed/ |
| `fastqc_trim.sh` | Quality check (trimmed) | FastQC | trimmed/*.fastq | fastqc_trimmed/ |
| `multiqc_trim.sh` | Aggregate QC (trimmed) | MultiQC | fastqc_trimmed/ | multiqc_trimmed/ |

---

## How to Run

### 1. Set up the environment

```bash
conda activate trim_env
```

### 2. Download data

```bash
bash scripts/getting_data.sh metadata/SraRunTable.csv metadata/selected_sra_run_ids
```

### 3. QC on raw reads

```bash
bash scripts/fastqc_raw.sh
bash scripts/multiqc_raw.sh
```

### 4. Trim

```bash
bash scripts/trimming.sh
```

### 5. QC on trimmed reads

```bash
bash scripts/fastqc_trim.sh
bash scripts/multiqc_trim.sh
```

---

## Study Design

| Parameter | Decision |
|-----------|----------|
| Research question | Do avian (duck) and human Influenza A H1N1 viruses differ in codon usage bias? |
| Host groups | Human (Homo sapiens) vs Duck (Anas platyrhynchos) |
| Subtype | H1N1 — present in both hosts, avoids host-subtype confounding |
| Genes | HA (segment 4), NA (segment 6), PB2 (segment 1) |
| Read type | Paired-end, Illumina |
| Mapping approach | Reference-based (not de novo assembly) |
| Sample size | ~100 total (50 duck + 50 human) from SRA raw reads |

---

## Pipeline Features

- **Skip logic** — every script checks if output already exists before re-running. Safe to re-run after a crash without redoing completed samples.
- **Error handling** — trimming failures are logged to `results/failed_samples.txt` and skipped; the pipeline continues.
- **Resumable** — all scripts read from accession lists, so interrupted runs pick up where they left off.

---

## Data Sources

- Raw reads: [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) — filtered to Influenza A virus, H1N1, paired-end, Illumina
- Host filter: strain name token search (`duck` / H1N1 strain names for human)
- Reference sequences: NCBI RefSeq, one per gene per subtype

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| SRA-tools (prefetch, fasterq-dump) | 2.10.0 | Download raw reads |
| FastQC | — | Per-sample quality reports |
| MultiQC | — | Aggregate quality reports |
| Trimmomatic | 0.40 | Adapter trimming |
| BWA | — | Reference alignment (next stage) |
| samtools | — | BAM processing (next stage) |
| bcftools | — | Variant calling + consensus (next stage) |

---

## Next Steps (after trimming)

1. BWA alignment to H1N1 gene references (HA, NA, PB2)
2. samtools — SAM to sorted indexed BAM
3. bcftools — variant calling + consensus sequence
4. Biopython — CDS extraction, codon usage analysis (RSCU, ENC, CAI, GC content)
5. Statistical comparison: avian vs human codon usage patterns

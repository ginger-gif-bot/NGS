# NGS Path
```
*raw reads (FASTQ) → FastQC/QC → trim → assemble → predict genes → extract CDS → analyze codons*
```
 

## File formats

1. *FASTQ*
- Raw sequence reads (raw output from the machine)
- Each read has four lines
    - Line1 : Identifier `(@Read_001)`
    - Line2 : Sequence `(ATGCTAGCTAGA)`
    - Line3: Separator `(+)`
    - Line4 : Quality Scores (ASCII) `(IIIHHHGGFF)`
- **Alignment softwares (like BWA, Bowtie2, HISAT2) take as input**

1. *FASTA*
- used to store reference genome in NGS pipeline
- Two components 
```
     >chr1
     ATGCTGACTAGCATAGACA
```

1. *SAM/BAM*
- **After aligning the FASTQ reads to the ref genome, the aligner produces a SAM/BAM file**
- SAM = Sequence Alignment/Map -> human readable text
- BAM = Binary Alignment/Map   -> compressed binary format (smaller size and much faster)
- SAM and BAM contain the same info
- Components of SAM file

```
     @HD VN:1.6 SO:coordinate
     @SQ SN:chr1 LN:248956422

     Read_001  0  chr1  1050  60  12M  *  0  0  ATCGTGACCTGA  IIIIHHHHGGFF
```
- Meaning 
     - Field         |          Meaning
     - Read_001      |          Read Name
     - 0             |          FLAG(Alignment information)
     - chr1          |          Reference genome
     - 1050          |          Position where the read aligned
     - 12M           |          CIGAR string (12 bases matched)
     - 12M           |          CIGAR string (12 bases matched)
     - 60            |          Mapping Quality
     - 12M           |          CIGAR string (12 bases matched)
     - "*"           |          Mate reference
     - 0             |          Mate position
     - 0             |          Tempalte length
     - ATCGTGACCTGA  |          Sequence
     - IIIIHHHHGGFF  |          Quality Scores

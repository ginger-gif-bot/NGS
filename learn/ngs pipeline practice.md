# NGS Pipeline: Raw Reads → Clean CDS (Practice Worksheet)

The capstone. This ties together everything — bash plumbing, file handling, and
biology — into the actual pipeline your thesis runs on: **take raw sequencing
reads and turn them into clean coding sequences (CDS) ready for codon-usage
analysis.**

Structured in three layers, per your request:
1. **The map** — the route from reads to CDS, and *why* each step exists.
2. **Each step** — concept first, then the real commands (FastQC, trimming,
   assembly, gene prediction, samtools).
3. **Wiring it** — assembling the steps into one runnable bash pipeline script.

Not every step needs heavy drilling — the *map* is the thing to internalize.

> **Read this first — the whole journey in one line:**
> ```
> *raw reads (FASTQ) → QC → trim → assemble → predict genes → extract CDS → analyze codons*
> ```
> Everything below is just expanding each arrow. Once you hold this shape in
> your head, no individual command is scary — you always know where you are.

> **A note on your two data types:** the *route* is nearly identical for
> Illumina (short reads) and Nanopore/PacBio (long reads) — QC, trim, assemble,
> predict, extract. Only the *tools* at a couple of steps differ (noted where it
> matters). So learn the map once; it serves both.

---

## PART 1 — The Map (understand the route before any command)

No commands here. Just make sure you can answer these in plain English. This is
the part that actually matters — the rest is looking up syntax.

1. **The file-format journey.** Your data changes format as it moves down the
   pipeline. Put these in order and say what each holds:
   `FASTA (.fasta)`, `FASTQ (.fastq)`, `SAM/BAM (.sam/.bam)`. Which one carries
   *quality scores*? Which is the raw sequencer output? Which is aligned reads?
2. **Why QC first?** Why do you *always* run FastQC on raw reads before doing
   anything else? What goes wrong downstream if you skip it?
3. **Why trim?** What two things does trimming remove, and why does leaving them
   in wreck an assembly? (Hint: adapters + low-quality ends.)
4. **Two routes to CDS — the fork in the road.** From clean reads there are two
   ways to get coding sequences:
   - **(A) Assembly route:** assemble reads into contigs/genome → predict genes
     → extract CDS. Used when you have *no reference*.
   - **(B) Reference route:** align reads to a known reference genome → call the
     coding regions. Used when a reference *exists*.
   For a codon-usage study on an organism *without* a good reference, which
   route do you take? Why?
5. **Why is FastQC run twice?** (Once on raw, once after trimming.) What are you
   checking the second time?
6. **Where does your Biopython codon work plug in?** After you have CDS FASTA
   files — that's the handoff point. The NGS pipeline's *job* is to deliver
   clean CDS; your codon-usage analysis takes over from there. Draw the boundary
   in your head: where does "NGS pipeline" end and "codon analysis" begin?
7. **Paired-end reads.** Illumina often gives you *two* files per sample
   (`_R1.fastq`, `_R2.fastq`) — forward and reverse reads of the same fragment.
   Why does the pipeline need to keep these two together at every step?

---

## PART 2 — The Steps (concept, then commands)

For each step: understand *what it does and why*, then the command. Don't
memorize flags — understand the shape. You'll always look up exact flags.

### Step 1 — Quality control (FastQC)

8. **Concept:** FastQC reads your FASTQ and produces a report of quality
   metrics — per-base quality, adapter contamination, GC content, duplication.
   It doesn't *change* your data; it *tells you what's wrong* so you know what to
   trim. Why is "look before you cut" the right order?
9. **Command:** run FastQC on a reads file:
   ```bash
   fastqc reads.fastq -o qc_reports/
   ```
   What does `-o` do? What files does FastQC produce?
10. **MultiQC:** with many samples, you run FastQC on each, then `multiqc` to
    combine all reports into one. Why is one combined report better than reading
    30 separate ones?

### Step 2 — Trimming (fastp / Trimmomatic)

11. **Concept:** trimming removes adapter sequences (leftover lab artifacts, not
    biology) and low-quality bases (usually at read ends, where the sequencer
    gets unreliable). The output is a cleaned FASTQ. Why must this happen
    *before* assembly, not after?
12. **Command (fastp — recommended, one tool does trim + report):**
    ```bash
    fastp -i reads.fastq -o trimmed.fastq
    ```
    For **paired** reads:
    ```bash
    fastp -i R1.fastq -I R2.fastq -o R1_trimmed.fastq -O R2_trimmed.fastq
    ```
    Why do the paired flags come in pairs (`-i/-I`, `-o/-O`)?
13. **The QC loop:** after trimming, run FastQC *again*. What are you confirming?
    (That the bad stuff is gone — *and* that you didn't over-trim into real
    sequence.)

### Step 3 — Assembly (SPAdes for Illumina / Flye for long reads)

14. **Concept:** individual reads are short fragments. Assembly stitches
    overlapping reads into longer continuous sequences (**contigs**), rebuilding
    the original genome/transcriptome as far as the data allows. Output is a
    FASTA of contigs. Why can't you predict genes directly on raw reads?
15. **Command (Illumina, SPAdes):**
    ```bash
    spades.py -1 R1_trimmed.fastq -2 R2_trimmed.fastq -o assembly/
    ```
    (Long reads use **Flye** instead: `flye --nano-raw reads.fastq -o
    assembly/`. Same *role*, different tool — this is the one real fork between
    your two data types.)
16. **Assembly QC:** how would you check the assembly is any good? (Think back to
    your bash skills: count contigs with `grep -c ">"`, check lengths.) What does
    a *good* assembly look like vs a fragmented one?

### Step 4 — Gene prediction / annotation (Prodigal / Augustus)

17. **Concept:** a contig is just a long DNA string. Gene prediction *finds the
    coding regions* within it — where genes start (ATG), stop, and the reading
    frame. Output includes the **CDS** (coding sequences) as FASTA. Why is this
    the step that actually gives you what codon analysis needs?
18. **Command (Prodigal — for prokaryotes/simple genomes):**
    ```bash
    prodigal -i assembly/contigs.fasta -d cds.fasta -a proteins.fasta
    ```
    What's the difference between the `-d` output (nucleotide CDS) and `-a`
    output (protein)? Which one feeds your codon-usage work, and why?
19. **In frame:** why is it *critical* that the CDS from this step are in-frame
    (length divisible by 3, starting at a real start codon)? Tie this back to the
    codon-column work — what breaks if the frame is off?

### Step 5 — Handling BAM/SAM (samtools) — the reference route

*(Only if you go the reference-alignment route instead of assembly. Skim if
you're doing pure assembly, but know it exists.)*

20. **Concept:** if a reference genome exists, you *align* reads to it (with
    `bwa` or `minimap2`) instead of assembling. The alignment is a SAM file
    (text) or BAM (compressed binary). `samtools` is the swiss-army knife for
    these. Why is BAM preferred over SAM for storage?
21. **Core samtools commands** — match each to its job:
    ```bash
    samtools view      # look at / filter alignments
    samtools sort      # sort by position (needed before indexing)
    samtools index     # build an index for fast random access
    samtools flagstat  # summary stats: how many reads aligned?
    ```
    Which one tells you "what fraction of my reads actually mapped"?
22. **The align→sort→index chain:** why must you *sort* a BAM before you can
    *index* it? (Same logic as `sort` before `uniq` — order enables the next
    step.)

---

## PART 3 — Wiring It Into One Pipeline Script

Now the payoff: everything you learned in the bash sheet — variables, `$1`,
loops, `basename`, `mkdir -p`, redirects, `set -e` — assembles into one runnable
pipeline. This is what "an NGS pipeline" literally *is*: a bash script that runs
the steps in order.

23. **The skeleton.** Write the *structure* (comments only, no real tools yet)
    of a script that takes a sample name as `$1` and walks it through: make
    output folders → QC → trim → QC again → assemble → predict genes → report.
    Just the comments and folder setup. This is the backbone.
24. **`set -e` — the safety net.** Put `set -e` at the top. Why is this
    non-negotiable in a pipeline? (If step 3 fails, you do NOT want step 4
    running on garbage — halt immediately.)
25. **Organized outputs.** Use `mkdir -p` to make a clean folder structure per
    sample (`results/$SAMPLE/qc`, `.../trimmed`, `.../assembly`, `.../cds`). Why
    does organizing outputs by sample matter when you have 30 samples?
26. **Progress messages.** Between steps, `echo` what's happening ("Trimming
    $SAMPLE...", "Assembly done"). Why does a long-running pipeline need loud
    progress output?
27. **The batch wrapper.** Write an *outer* loop that runs your whole
    single-sample pipeline over every sample in a folder. This is the "process
    all samples" move — the entire reason pipelines exist.
28. **Argument guard.** Add the `if [ -z "$1" ]` usage check from the bash sheet
    — refuse to run with no sample name. Why does every real pipeline start with
    input validation?

---

## PART 4 — Mini-projects (the real shape)

29. **Dry-run pipeline:** write the complete single-sample script with all the
    structure (folders, `set -e`, echoes, guard) but with the actual tool
    commands as `echo "would run: fastqc..."` placeholders. Run it — confirm the
    *flow* works before you ever touch a real tool. (This is genuinely how people
    debug pipelines — get the plumbing right first.)
30. **CDS delivery check:** after the (real or simulated) pipeline, write a bash
    check that confirms you got usable CDS: count sequences (`grep -c ">"`),
    verify they're in frame (length ÷ 3), and report. This is the handoff
    validation before Biopython takes over.
31. **The full picture:** in a comment block at the top of your script, draw the
    ASCII flow of your whole thesis: `raw reads → [NGS pipeline] → clean CDS →
    [Biopython codon analysis] → codon usage bias results`. Mark exactly where
    this script's job ends.

---

---

# ✅ Answer Key & Explanations

### Part 1 — The Map

**Q1 — format journey:**
- **FASTQ** = raw sequencer output. Sequence **+ quality scores** per base.
- **FASTA** = sequences only, no quality. (Assembled contigs, CDS, proteins.)
- **SAM/BAM** = reads *aligned to a reference*, with position info.
Order: FASTQ (raw) → [align] → SAM/BAM, or FASTQ → [assemble] → FASTA.
FASTQ carries quality scores; that's the whole reason it exists.

**Q2:** QC first because every downstream step assumes decent data. Bad reads →
bad assembly → wrong genes → wrong codons. FastQC tells you what to fix *before*
errors propagate and silently corrupt everything.

**Q3:** Trimming removes (1) **adapters** (lab artifacts, not biology) and (2)
**low-quality bases** (sequencer gets unreliable at read ends). Left in, they
create false overlaps and errors that break or mislead assembly.

**Q4:** For an organism **without a reference**, take the **assembly route (A)** —
you have nothing to align to, so you build the sequence from the reads
themselves, then predict genes. The reference route needs a reference to exist.

**Q5:** Second FastQC confirms trimming *worked* (adapters gone, quality up)
*and* that you didn't **over-trim** — cut so aggressively you removed real
biological sequence. It's a before/after check.

**Q6:** The NGS pipeline ends when you have **clean CDS FASTA files**. That's the
handoff. Your Biopython codon-column/usage work begins there. Everything in this
sheet exists to deliver that one deliverable: in-frame CDS.

**Q7:** Paired reads (R1/R2) are two ends of the *same* DNA fragment — together
they give more information (the gap between them, orientation). Tools use both;
if they get separated or misordered, the pairing information is lost and
alignment/assembly quality drops.

### Part 2 — Steps

**Q8/Q9/Q10:**
```bash
fastqc reads.fastq -o qc_reports/     # -o = output folder for the report
```
FastQC produces an `.html` report + a `.zip` of the raw metrics. **MultiQC**
(`multiqc qc_reports/`) merges many reports into one — essential at 30 samples,
because you spot *systematic* problems (all samples share an issue) that you'd
miss reading reports one at a time.

**Q11/Q12/Q13:**
```bash
fastp -i R1.fastq -I R2.fastq -o R1_trim.fastq -O R2_trim.fastq
```
Paired flags come in pairs because the two files must be trimmed *together* to
stay synchronized. Trim *before* assembly so the assembler never sees adapters
or junk. Re-run FastQC after to confirm clean-but-not-over-trimmed.

**Q14/Q15/Q16:**
```bash
spades.py -1 R1_trim.fastq -2 R2_trim.fastq -o assembly/   # Illumina
flye --nano-raw reads.fastq -o assembly/                    # long reads
```
You can't predict genes on raw reads because reads are too short to contain
whole genes — you need the assembled, continuous contigs first. Check assembly
with `grep -c ">" assembly/contigs.fasta` (contig count) — fewer, longer contigs
= better assembly; thousands of tiny ones = fragmented.

**Q17/Q18/Q19:**
```bash
prodigal -i contigs.fasta -d cds.fasta -a proteins.fasta
```
`-d` = **nucleotide CDS** (what codon analysis needs — codons are nucleotide
triplets). `-a` = protein translation. You feed the **`-d` nucleotide CDS** to
your codon work. It's critical they're in-frame because your codon-column
slicing (`[k*3:k*3+3]`) assumes correct frame — off by one nucleotide and every
codon is wrong. This is the exact frame lesson from the codon sheet, now at the
pipeline level.

**Q20/Q21/Q22:**
BAM is compressed binary — far smaller and faster than text SAM.
- `samtools flagstat` → **"what fraction of reads mapped"** (the mapping-rate check)
- `sort` before `index` because indexing needs position-ordered data to build
  its fast-lookup table — same "order enables the next step" logic as `sort |
  uniq`.

### Part 3 — Wiring (skeleton answer)

```bash
#!/bin/bash
set -e                                    # Q24: halt on any failure

if [ -z "$1" ]; then                      # Q28: input guard
    echo "Usage: $0 <sample_name>"
    exit 1
fi

SAMPLE=$1
OUT=results/$SAMPLE
mkdir -p "$OUT"/{qc,trimmed,assembly,cds} # Q25: organized per-sample folders

echo "=== [$SAMPLE] QC on raw reads ==="  # Q26: progress
fastqc data/${SAMPLE}.fastq -o "$OUT"/qc/

echo "=== [$SAMPLE] Trimming ==="
fastp -i data/${SAMPLE}.fastq -o "$OUT"/trimmed/${SAMPLE}.fastq

echo "=== [$SAMPLE] QC on trimmed reads ==="
fastqc "$OUT"/trimmed/${SAMPLE}.fastq -o "$OUT"/qc/

echo "=== [$SAMPLE] Assembly ==="
spades.py -s "$OUT"/trimmed/${SAMPLE}.fastq -o "$OUT"/assembly/

echo "=== [$SAMPLE] Gene prediction ==="
prodigal -i "$OUT"/assembly/contigs.fasta -d "$OUT"/cds/${SAMPLE}_cds.fasta

echo "=== [$SAMPLE] Done. CDS: $OUT/cds/${SAMPLE}_cds.fasta ==="
```

**Q27 — batch wrapper** (runs the above for every sample):
```bash
for f in data/*.fastq; do
    name=$(basename "$f" .fastq)          # basename again — same move as bash sheet
    bash pipeline.sh "$name"
done
```

### Part 4

**Q29 — dry run:** replace each tool line with `echo "would run: <tool> ..."` and
run it. If the folders get made and the echoes fire in order, your plumbing is
correct — *then* swap in real tools one at a time.

**Q30 — CDS check:**
```bash
echo "sequences: $(grep -c '>' cds.fasta)"
# frame check: every CDS length divisible by 3 (loop or awk)
```

**Q31 — the full picture:**
```
raw reads (FASTQ)
      │  ┌─────────── THIS SCRIPT'S JOB ───────────┐
      ▼  │ QC → trim → QC → assemble → predict genes │
   clean CDS (FASTA) ◄───────────────────────────────┘
      │
      ▼  [Biopython: codon columns, usage counting]
   codon usage bias results
```

---

## The one thing to actually carry
**An NGS pipeline is not magic — it's a bash script that runs known tools in a
fixed order, each turning one file format into the next, until you have clean
CDS.** You already have every skill it needs:
- **bash** (variables, loops, `$1`, `basename`, `mkdir -p`, `set -e`) → the wiring
- **file-format sense** (FASTQ → FASTA → CDS) → knowing where you are
- **Biopython** → what happens *after* the pipeline delivers CDS

The tools (fastqc, fastp, spades, prodigal, samtools) are just commands you look
up and slot into the bash skeleton you already know how to write. The map is the
hard part, and now you have it. Start with the **dry-run** (Q29) — get the flow
right with echoes, then drop in real tools one step at a time, checking output at
each stage. That's how every pipeline gets built.

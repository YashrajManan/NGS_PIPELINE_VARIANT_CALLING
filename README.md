# NGS Pipeline + Variant Calling — NA12878, Validated Against GIAB

A from-scratch NGS pipeline — raw reads → alignment → duplicate marking → variant calling →
filtering — run on **real human sequencing data**, then validated against an independent,
expert-curated truth set. Real data at every step: reads and reference are fetched **live** from
public genomic servers (not packaged files), and the final calls are checked against NIST's Genome in a
Bottle (GIAB) benchmark for the same sample.

## Aim

Can a standard statistical NGS pipeline (alignment plus genotype-likelihood variant calling) recover
real human genetic variants from raw reads — and how accurate are the calls when checked against an
independent, multi-technology expert truth set, not just trusted on faith?

## Objective

Build the complete real NGS variant-calling pipeline end to end (read QC, alignment, duplicate marking,
statistical variant calling, hard filtering) on real human sequencing reads, then quantitatively measure
precision, recall, and F1 against NIST's Genome in a Bottle high-confidence truth set for the identical
sample and genomic region.

## Data fetch

All real data, fetched live, nothing packaged. **Reads**: NA12878 (the single most validated human
genome in existence — used to benchmark essentially every sequencing method and variant caller ever
published), 1000 Genomes phase-3 low-coverage CRAM (~4-8x genome-wide), remote-sliced directly to one
100kb window of chromosome 20 (`chr20:10,000,000-10,100,000`) over HTTP — no whole-genome download.
**Reference**: GRCh38 chromosome 20 (UCSC mirror of the primary assembly), downloaded directly. **Truth
set**: NIST/Genome in a Bottle's high-confidence benchmark VCF and confident-regions BED for NA12878 — an
answer key built by reconciling multiple independent sequencing technologies (Illumina, PacBio, Nanopore,
and more), not one machine's output trusted alone.

## Data describe

Real Illumina paired-end short reads (~150bp each), remote-sliced from the original low-coverage 1000
Genomes CRAM alignment for NA12878 — meaning this window's average depth is genuinely low (~4-8x), not a
higher-coverage 30x resequencing, which directly shapes the pipeline's results (see below). GIAB's truth
set additionally ships a confident-regions BED: the subset of the genome where even GIAB's own reconciled
truth set is considered reliable, since some repetitive/duplicated regions can't be confidently resolved
by any current technology — comparing calls to truth only within that BED is essential to avoid
manufacturing false positives from regions where the truth set itself is unreliable.

## Methods / Workflow — what we did

```
Real NA12878 CRAM (remote-sliced) ── real GRCh38 chr20 reference (downloaded)
   │
   ▼ fastp (read QC/trim)
   ▼ BWA-MEM (align) → samtools sort/index
   ▼ samtools fixmate/markdup (duplicate marking)
   ▼ bcftools mpileup | bcftools call -mv (statistical variant calling)
   ▼ bcftools filter (QUAL/DP hard filters)
   ▼ restrict to GIAB confident regions → bcftools isec vs. GIAB truth
   → precision / recall / F1  [validated against real ground truth]
```

1. Slice the real NA12878 CRAM to a 100kb chr20 window and fetch the matching GRCh38 reference.
2. QC/trim reads with `fastp`.
3. Align with **BWA-MEM**, sort and index with `samtools`.
4. Mark PCR duplicates (`samtools fixmate`/`markdup`) — two reads from the same original DNA molecule are
   one piece of evidence, not two, and failing to mark them inflates confidence in variants supported by
   fewer independent observations than the raw read count suggests.
5. Call variants statistically with `bcftools mpileup | bcftools call -mv` — computing the likelihood of
   each possible genotype (homozygous-reference, heterozygous, homozygous-alternate) from the pileup's
   base calls and quality scores, not a raw vote count.
6. Apply hard filters on QUAL (statistical confidence) and depth (DP) to separate believable calls from
   noise.
7. Restrict both the calls and the truth set to GIAB's confident regions, then compare with
   `bcftools isec` to compute true positives, false positives, and false negatives.

## Results

| Metric | Value |
|---|---|
| True Positives | 99 |
| False Positives | 18 |
| False Negatives | 110 |
| **Precision** | **0.846** |
| **Recall** | **0.474** |
| **F1** | **0.607** |

## Biology interpretation of results

Precision is high — when this pipeline calls a variant, it's real 85% of the time, meaning the
statistical genotype-calling and hard-filtering steps are doing their job: they are not flooding the
output with spurious calls. Recall is the real, honestly reported weak point: at only ~5x average depth
in this specific window (the original low-coverage 1000 Genomes alignment, not a deeper 30x
resequencing), many true heterozygous variants simply don't have enough overlapping reads for the caller
to call them with confidence, even after loosening the depth filter to match the lower coverage. This is
the textbook, well-established coverage/recall relationship in variant calling, not a defect in the
pipeline's logic — more sequencing depth would recover more of the missed variants, and the gap between
precision (0.846) and recall (0.474) directly diagnoses *which* failure mode dominates here: the
pipeline is too conservative (under-calling from shallow coverage), not too permissive. This
demonstrates, on real data with a real external answer key, that a fully statistical pipeline — no manual
curation of any call — can recover real human genetic variation with quantifiable, bounded,
externally-validated error, which is exactly the kind of pipeline underlying real clinical genetic
diagnosis, population genetics, and cancer genomics variant calling at scale.

## Learning through project

Precision and recall diagnose two structurally different failure modes, and reporting only one (or only
F1) would hide which one actually applies — here, the gap between a strong 0.846 precision and a weak
0.474 recall specifically points to insufficient sequencing depth, not to the caller or filters being
poorly tuned, and that diagnosis only becomes visible by looking at both numbers together, not the
combined F1 alone. Ground-truth validation is only meaningful within the region the truth set itself
considers reliable — comparing calls against GIAB's truth set genome-wide, rather than restricted to its
own published confident regions, would manufacture false positives out of areas where even the expert
truth set can't confidently resolve variants, unfairly penalizing the pipeline for a limitation that
isn't its own. More generally: validating a real bioinformatics pipeline against an independent,
multi-technology expert truth set — rather than trusting the pipeline's own internal confidence scores —
is what actually earns confidence in a result, the same principle applied differently across this whole
portfolio (cross-language validation, cross-cohort external validation, and here, cross-technology
ground-truth validation).

## Limitations

A single 100kb window, not genome-wide — this precision/recall should not be extrapolated to whole-genome
performance. Hard filters (QUAL/DP thresholds) are cruder than model-based filtering (e.g. GATK's VQSR).
Results depend on the specific aligner/caller/parameter choices made, not on the underlying biology
alone. The low-coverage source data specifically depresses recall here, as directly observed and
explained above, rather than reflecting a universal ceiling on this pipeline's achievable recall at
higher depth.

## Reproduce

Requires WSL/Linux with `sudo` access (installs `bwa`, `samtools`, `bcftools`, `fastp`, `tabix` via `apt`
on first run). Needs internet access throughout (live data fetch from EBI/NCBI/UCSC).

```bash
bash ngs_variant_calling.sh
```

## Tech

`BWA-MEM` · `samtools` · `bcftools` · `fastp` · 1000 Genomes CRAM streaming · GIAB benchmark validation

## Files

```
ngs_variant_calling.sh   # the full pipeline (fetch → align → call → filter → validate)
results/                 # fastp report, flagstat, VCFs, precision_recall.txt, isec/
```

## License

All rights reserved — see `LICENSE`. This repository is public for portfolio/demonstration purposes
only; no permission is granted to copy, modify, or reuse any part of it.

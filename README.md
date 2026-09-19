# NGS Pipeline + Variant Calling — NA12878, validated against GIAB

A from-scratch NGS pipeline — raw reads → alignment → duplicate marking → variant calling →
filtering — run on **real human sequencing data**, then validated against an independent,
expert-curated truth set. Real data at every step: reads and reference are fetched **live** from
public genomic servers (not packaged files), and the final calls are checked against NIST's
Genome in a Bottle (GIAB) benchmark for the same sample.

## Research question

Can a standard statistical NGS pipeline (alignment + genotype-likelihood variant calling)
recover real human genetic variants from raw reads — and how do the calls compare to an
independent, multi-technology expert truth set?

## Data (all real, fetched live)

- **Reads:** NA12878 (the most-validated human genome in existence), 1000 Genomes phase-3
  low-coverage CRAM (~4-8x genome-wide), remote-sliced to one **100kb window of chr20**
  (`chr20:10,000,000-10,100,000`) — no whole-genome download.
- **Reference:** GRCh38 chr20 (UCSC mirror of the primary assembly), downloaded directly.
- **Truth set:** NIST/Genome in a Bottle's high-confidence benchmark VCF + confident-regions
  BED for NA12878 — an answer key built from multiple independent sequencing technologies.

## Pipeline

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

## Key result

| Metric | Value |
|---|---|
| True Positives | 99 |
| False Positives | 18 |
| False Negatives | 110 |
| **Precision** | **0.846** |
| **Recall** | **0.474** |
| **F1** | **0.607** |

Precision is high — when the pipeline calls a variant, it's real 85% of the time. Recall is the
weak point: at only ~5x average depth in this window (this CRAM is the *original* ~4-8x
low-coverage 1000 Genomes alignment, not a 30x resequencing), many true heterozygous variants
don't have enough overlapping reads to be called confidently — even after loosening the depth
filter (DP<4) to match the lower coverage. This is the textbook coverage/recall relationship in
variant calling, not a pipeline defect: **more depth would recover more of the missed variants**.

## Limitations

Single 100kb window, not genome-wide — don't extrapolate this precision/recall to whole-genome
performance. Hard filters (QUAL/DP thresholds) are cruder than model-based filtering (GATK's
VQSR). Results depend on the specific aligner/caller/parameters chosen, not biology alone.
Low-coverage source data specifically depresses recall, as observed here.

## Files

```
ngs_variant_calling.sh   # the full pipeline (fetch → align → call → filter → validate)
results/                 # fastp report, flagstat, VCFs, precision_recall.txt, isec/
```

## Run

Requires WSL/Linux with `sudo` access (installs `bwa`, `samtools`, `bcftools`, `fastp`, `tabix`
via `apt` on first run). Needs internet access throughout (live data fetch from EBI/NCBI/UCSC).

```bash
bash ngs_variant_calling.sh
```

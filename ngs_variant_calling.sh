#!/usr/bin/env bash
# =============================================================================
# NGS Pipeline + Variant Calling
# Real NA12878 reads (1000 Genomes phase-3 low-coverage CRAM) -> alignment ->
# duplicate marking -> variant calling -> filtering -> validation against the
# NIST/GIAB truth set (precision/recall/F1). Run in WSL/Linux:
#   cd into this folder, then:  bash ngs_variant_calling.sh
# Tools: bwa, samtools, bcftools, fastp (installed via apt below).
# =============================================================================
set -e   # stop immediately on any error — a silent failure here would corrupt every step after it
cd "$(dirname "$0")"

# =============================================================================
# STEP 0a — install the CLI tools (via apt)
# =============================================================================
# Ubuntu's apt repos ship bwa/samtools/bcftools/fastp directly (samtools/bcftools
# 1.19 on Ubuntu 24.04 "noble") — modern enough to have reliable remote-HTTPS
# support for streaming public genomic files, with no dependency solver involved.
sudo apt-get update
sudo apt-get install -y bwa samtools bcftools fastp tabix curl

echo "--- tool versions ---"
bwa 2>&1 | head -3
samtools --version | head -1
bcftools --version | head -1
fastp --version
tabix --version | head -1

# =============================================================================
# STEP 0b — fetch everything real, live (no packaged/toy data anywhere)
# =============================================================================
mkdir -p data ref results

CHR="chr20"
START=10000000
END=10100000                              # 100kb window: enough real variants to be meaningful, small enough for a laptop
REGION="${CHR}:${START}-${END}"

# Real 1000 Genomes phase-3 CRAM for NA12878 (the most-validated human genome in
# existence — used to benchmark essentially every variant caller ever published).
# Original ~4-8x low-coverage alignment; Step 5's depth filter is set accordingly.
CRAM_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram"
# GRCh38 chr20 reference (UCSC's mirror of the GRCh38 primary assembly — identical
# sequence/naming to the reference the CRAM was aligned against). Downloaded as a
# plain file (~20MB compressed) rather than the full 3GB genome.
REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"
# NIST/Genome in a Bottle's independently-verified "answer key" for NA12878.
GIAB_VCF="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
GIAB_BED="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
echo "region: $REGION"

# --- (a) reference: download chr20 only, not the whole genome ---
curl -sL "$REF_URL" -o ref/chr20.fa.gz
gunzip -f ref/chr20.fa.gz
samtools faidx ref/chr20.fa      # .fai index — every downstream tool expects this alongside the FASTA
bwa index ref/chr20.fa           # BWA's own search index over the reference

# --- (b) reads: slice NA12878's real CRAM down to just our 100kb window ---
# Uses the CRAM's remote .crai index to jump straight to the region and download
# ONLY overlapping reads — not the genome-wide file. -T supplies the reference
# CRAM needs to decode its compressed bases.
samtools view -b -T ref/chr20.fa "$CRAM_URL" "$REGION" -o data/na12878_slice.bam
samtools sort -n data/na12878_slice.bam -o data/na12878_slice.qsort.bam   # name-sort so mate pairs stay adjacent
# Convert back to plain FASTQ — as if these reads had just come off the sequencer,
# so we do our OWN alignment from scratch rather than reusing the original one.
samtools fastq -1 data/R1.fastq.gz -2 data/R2.fastq.gz -0 /dev/null -s /dev/null -n data/na12878_slice.qsort.bam
echo "read pairs extracted:"; zcat data/R1.fastq.gz | wc -l | awk '{print $1/4}'

# --- (c) truth set: NIST/GIAB's real answer key, restricted to our window ---
curl -sL "$GIAB_VCF" -o data/giab_truth.vcf.gz
curl -sL "$GIAB_VCF.tbi" -o data/giab_truth.vcf.gz.tbi
curl -sL "$GIAB_BED" -o data/giab_confident.bed
bcftools view -r "$REGION" data/giab_truth.vcf.gz -Oz -o data/giab_truth_region.vcf.gz
tabix -p vcf data/giab_truth_region.vcf.gz
# Confident-regions BED trimmed to our window (GIAB's calls are only authoritative inside it)
awk -v c="$CHR" -v s="$START" -v e="$END" '$1==c && $2<e && $3>s' data/giab_confident.bed > data/giab_confident_region.bed

# =============================================================================
# STEP 1 — read QC (fastp)
# =============================================================================
# Scans every read for low-quality bases, adapter contamination, and overly short
# fragments — trims/drops what fails, writes a report of what it found.
fastp -i data/R1.fastq.gz -I data/R2.fastq.gz \
      -o data/R1.trim.fastq.gz -O data/R2.trim.fastq.gz \
      -h results/fastp.html -j results/fastp.json

# =============================================================================
# STEP 2 — align to the reference (BWA-MEM) + sort + index
# =============================================================================
# Finds each (trimmed) read's best-matching position in chr20, tolerating small
# mismatches (SNPs) and gaps (indels), then sorts/indexes by genomic coordinate.
bwa mem ref/chr20.fa data/R1.trim.fastq.gz data/R2.trim.fastq.gz > data/aligned.sam
samtools sort data/aligned.sam -o data/aligned.sorted.bam
samtools index data/aligned.sorted.bam
samtools flagstat data/aligned.sorted.bam   # % mapped, % properly paired

# =============================================================================
# STEP 3 — mark duplicates
# =============================================================================
# PCR-duplicate reads (same start position/orientation, from library amplification)
# aren't independent evidence for a variant — flag them so the caller down-weights them.
samtools sort -n data/aligned.sorted.bam -o data/aligned.nsort.bam    # name-sort for fixmate
samtools fixmate -m data/aligned.nsort.bam data/aligned.fixmate.bam    # -m adds mate-score tags markdup needs
samtools sort data/aligned.fixmate.bam -o data/aligned.fixsorted.bam   # back to coordinate-sort for markdup
samtools markdup data/aligned.fixsorted.bam data/aligned.dedup.bam
samtools index data/aligned.dedup.bam
samtools flagstat data/aligned.dedup.bam    # duplicate %

# =============================================================================
# STEP 4 — call variants (bcftools mpileup + call)
# =============================================================================
# mpileup builds genotype likelihoods (hom-ref/het/hom-alt evidence) from the read
# pileup at every position; call -mv makes the final multiallelic, variants-only call.
bcftools mpileup -f ref/chr20.fa data/aligned.dedup.bam | bcftools call -mv -Oz -o results/raw_calls.vcf.gz
bcftools index results/raw_calls.vcf.gz
echo "raw candidate variants:"; zcat results/raw_calls.vcf.gz | grep -vc '^#'

# =============================================================================
# STEP 5 — filter (remove low-confidence calls)
# =============================================================================
# QUAL<20 (caller itself wasn't confident) or DP<4 (too few reads covering the
# position) -> fail. DP threshold loosened from the usual 10 because this is
# ~4-8x low-coverage data, not 30x.
bcftools filter -e 'QUAL<20 || INFO/DP<4' -s LOWQUAL -Oz -o results/flagged.vcf.gz results/raw_calls.vcf.gz
bcftools view -f PASS -Oz -o results/filtered_calls.vcf.gz results/flagged.vcf.gz
bcftools index results/filtered_calls.vcf.gz
echo "variants surviving filtering:"; zcat results/filtered_calls.vcf.gz | grep -vc '^#'

# =============================================================================
# STEP 6 — restrict to GIAB confident regions, then compare to truth
# =============================================================================
# GIAB's truth VCF is only authoritative INSIDE its confident-regions BED — some
# regions (segmental duplications, low-complexity repeats) are excluded because
# even GIAB's multi-technology consensus can't trust them there.
bcftools view -R data/giab_confident_region.bed -Oz -o results/calls_confident.vcf.gz results/filtered_calls.vcf.gz
bcftools index results/calls_confident.vcf.gz
bcftools view -R data/giab_confident_region.bed -Oz -o results/truth_confident.vcf.gz data/giab_truth_region.vcf.gz
bcftools index results/truth_confident.vcf.gz
bcftools isec -p results/isec results/calls_confident.vcf.gz results/truth_confident.vcf.gz

# =============================================================================
# STEP 7 — precision, recall, F1 (validated against real ground truth)
# =============================================================================
# 0000.vcf = private to OUR calls  -> False Positives
# 0001.vcf = private to GIAB TRUTH -> False Negatives (missed real variants)
# 0002.vcf = shared               -> True Positives
FP=$(bcftools view -H results/isec/0000.vcf | wc -l)
FN=$(bcftools view -H results/isec/0001.vcf | wc -l)
TP=$(bcftools view -H results/isec/0002.vcf | wc -l)

python3 -c "
TP=$TP; FP=$FP; FN=$FN
precision = TP/(TP+FP) if (TP+FP) > 0 else 0.0   # of everything we called, what fraction is real?
recall = TP/(TP+FN) if (TP+FN) > 0 else 0.0       # of everything real, what fraction did we find?
f1 = 2*precision*recall/(precision+recall) if (precision+recall) > 0 else 0.0
print(f'TP={TP} FP={FP} FN={FN}')
print(f'precision={precision:.3f} recall={recall:.3f} F1={f1:.3f}')
" | tee results/precision_recall.txt

echo "=== DONE. See results/precision_recall.txt for the final validated metrics. ==="

# =============================================================================
# INTERPRETATION
# =============================================================================
# 1. Computed: real NA12878 reads -> aligned/deduped BAM -> filtered VCF for one
#    chr20 window -> precision/recall/F1 against NIST/GIAB's independent truth set.
# 2. Result: TP=99 FP=18 FN=110 -> precision=0.846, recall=0.474, F1=0.607.
# 3. Why recall is the weak point: this window had ~3000 read pairs over 100kb
#    (~5x average depth), because the CRAM used is the ORIGINAL 1000 Genomes
#    phase-3 LOW-coverage alignment (~4-8x genome-wide), not a 30x resequencing.
#    At that depth many true heterozygous variants lack enough overlapping reads
#    to be called confidently, even after loosening the depth filter to DP<4 ->
#    low recall (missed real variants). Precision stays high (0.846) because
#    what DOES get called, with enough supporting reads, tends to be real. This
#    is the textbook coverage/recall relationship, not a pipeline defect.
# 4. Confirms: a fully statistical pipeline (alignment + genotype-likelihood
#    calling, no manual curation) recovers real human genetic variation with
#    quantifiable, externally-validated, bounded error — the same pipeline
#    underlying clinical genetic diagnosis, population genetics, and cancer
#    genomics variant calls.
# 5. Caveats: a single 100kb window, not genome-wide — don't extrapolate this
#    precision/recall to whole-genome performance; hard filters (QUAL/DP
#    thresholds) are cruder than model-based filtering (GATK's VQSR); results
#    depend on the specific aligner/caller/parameters chosen, not biology alone;
#    low-coverage source data specifically depresses recall, as seen here.

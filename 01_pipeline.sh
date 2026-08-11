#!/usr/bin/env bash
# ==============================================================================
# 01_pipeline.sh — Whole Exome Sequencing: FASTQ -> ClinVar-annotated VCF
#
# Usage:
#   ./01_pipeline.sh -r <reference.fa> -1 <R1.fastq[.gz]> -2 <R2.fastq[.gz]> \
#                     -c <clinvar.vcf.gz> -s <sample_name> [-t <threads>] \
#                     [-b <target_regions.bed>] [-o <output_root>]
#
# -b is optional: a BED file of your exome capture kit's target regions
#    (Agilent SureSelect / IDT xGen / Twist Exome, etc). If given, variant
#    calling and reporting are restricted to on-target regions.
# ==============================================================================
set -euo pipefail

THREADS=4
SORT_MEM="512M"
OUTROOT="results"
TARGET_BED=""

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while getopts "r:1:2:c:s:t:b:o:h" opt; do
  case "${opt}" in
    r) REF="${OPTARG}" ;;
    1) R1="${OPTARG}" ;;
    2) R2="${OPTARG}" ;;
    c) CLINVAR="${OPTARG}" ;;
    s) SAMPLE="${OPTARG}" ;;
    t) THREADS="${OPTARG}" ;;
    b) TARGET_BED="${OPTARG}" ;;
    o) OUTROOT="${OPTARG}" ;;
    h|*) usage ;;
  esac
done

: "${REF:?-r <reference.fa> required}"
: "${R1:?-1 <R1.fastq> required}"
: "${R2:?-2 <R2.fastq> required}"
: "${CLINVAR:?-c <clinvar.vcf.gz> required}"
: "${SAMPLE:?-s <sample_name> required}"

for f in "${REF}" "${R1}" "${R2}" "${CLINVAR}"; do
  [ -f "${f}" ] || { echo "ERROR: missing file: ${f}"; exit 1; }
done

OUTDIR="${OUTROOT}/${SAMPLE}"
mkdir -p "${OUTDIR}"/{qc,align,vcf,report}

echo "=== Sample: ${SAMPLE} | Threads: ${THREADS} | Target BED: ${TARGET_BED:-none (whole reference)} ==="

echo "[1/8] Reference indices"
[ -f "${REF}.bwt" ] || bwa index "${REF}"
[ -f "${REF}.fai" ] || samtools faidx "${REF}"
DICT="${REF%.*}.dict"
[ -f "${DICT}" ] || gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"

echo "[2/8] QC + adapter trimming (fastp)"
# --detect_adapter_for_pe: auto-detect Illumina adapters from PE overlap
# -q 20: min base quality kept as high-quality        -l 50: drop reads <50bp post-trim
# -5 -3: sliding-window quality trim from both ends
fastp \
  -i "${R1}" -I "${R2}" \
  -o "${OUTDIR}/qc/${SAMPLE}_R1.trim.fastq.gz" \
  -O "${OUTDIR}/qc/${SAMPLE}_R2.trim.fastq.gz" \
  --detect_adapter_for_pe \
  -q 20 -l 50 -5 -3 \
  --thread "${THREADS}" \
  -j "${OUTDIR}/qc/${SAMPLE}_fastp.json" \
  -h "${OUTDIR}/qc/${SAMPLE}_fastp.html"

echo "[3/8] Alignment (BWA-MEM) + coordinate sort"
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib1"
bwa mem -t "${THREADS}" -R "${RG}" "${REF}" \
  "${OUTDIR}/qc/${SAMPLE}_R1.trim.fastq.gz" \
  "${OUTDIR}/qc/${SAMPLE}_R2.trim.fastq.gz" \
  | samtools sort -@ "${THREADS}" -m "${SORT_MEM}" \
    -o "${OUTDIR}/align/${SAMPLE}.sorted.bam" -
samtools index "${OUTDIR}/align/${SAMPLE}.sorted.bam"

echo "[3b/8] Alignment QC"
samtools flagstat "${OUTDIR}/align/${SAMPLE}.sorted.bam" | tee "${OUTDIR}/align/${SAMPLE}.flagstat.txt"
if [ -n "${TARGET_BED}" ]; then
  samtools depth -a -b "${TARGET_BED}" "${OUTDIR}/align/${SAMPLE}.sorted.bam" \
    | awk '{s+=$3;n++} END{if(n>0) printf "  On-target mean depth: %.2fx\n", s/n}'
else
  samtools depth -a "${OUTDIR}/align/${SAMPLE}.sorted.bam" \
    | awk '{s+=$3;n++} END{if(n>0) printf "  Mean depth: %.2fx\n", s/n}'
fi

echo "[4/8] Mark duplicates (PCR/optical dups from library prep)"
gatk MarkDuplicates \
  -I "${OUTDIR}/align/${SAMPLE}.sorted.bam" \
  -O "${OUTDIR}/align/${SAMPLE}.dedup.bam" \
  -M "${OUTDIR}/align/${SAMPLE}.dup_metrics.txt"
samtools index "${OUTDIR}/align/${SAMPLE}.dedup.bam"

echo "[5/8] Variant calling (GATK HaplotypeCaller, diploid germline)"
HC_ARGS=(-R "${REF}" -I "${OUTDIR}/align/${SAMPLE}.dedup.bam" \
  -O "${OUTDIR}/vcf/${SAMPLE}.raw.vcf.gz" --sample-ploidy 2)
[ -n "${TARGET_BED}" ] && HC_ARGS+=(-L "${TARGET_BED}" -ip 20)  # -ip 20: pad 20bp around targets
gatk HaplotypeCaller "${HC_ARGS[@]}"

echo "[6/8] Hard-filter + normalize"
# QD<2.0: quality normalized by depth, low = likely artifact
# FS>60.0: strand bias (Fisher strand), high = likely artifact
# MQ<40.0: mapping quality, low = ambiguous/multi-mapped region
# DP<10:   read depth, below this genotype confidence is unreliable
gatk VariantFiltration \
  -R "${REF}" \
  -V "${OUTDIR}/vcf/${SAMPLE}.raw.vcf.gz" \
  -O "${OUTDIR}/vcf/${SAMPLE}.filtered.vcf.gz" \
  --filter-expression "QD < 2.0"  --filter-name "lowQD" \
  --filter-expression "FS > 60.0" --filter-name "highFS" \
  --filter-expression "MQ < 40.0" --filter-name "lowMQ" \
  --filter-expression "DP < 10"   --filter-name "lowDP"

bcftools view -f PASS "${OUTDIR}/vcf/${SAMPLE}.filtered.vcf.gz" -Oz \
  -o "${OUTDIR}/vcf/${SAMPLE}.pass.vcf.gz"
bcftools norm -f "${REF}" -m -both "${OUTDIR}/vcf/${SAMPLE}.pass.vcf.gz" -Oz \
  -o "${OUTDIR}/vcf/${SAMPLE}.norm.vcf.gz"
tabix -p vcf "${OUTDIR}/vcf/${SAMPLE}.norm.vcf.gz"

echo "[7/8] Annotate against ClinVar"
SnpSift annotate \
  -info CLNSIG,CLNDN,CLNREVSTAT,CLNHGVS,CLNDISDB,CLNVC,GENEINFO,MC,ORIGIN \
  "${CLINVAR}" \
  "${OUTDIR}/vcf/${SAMPLE}.norm.vcf.gz" \
  > "${OUTDIR}/vcf/${SAMPLE}.clinvar_annotated.vcf"

echo "[8/8] Summary"
N_TOTAL=$(grep -vc "^#" "${OUTDIR}/vcf/${SAMPLE}.clinvar_annotated.vcf" || true)
N_ANNOT=$(grep -v "^#" "${OUTDIR}/vcf/${SAMPLE}.clinvar_annotated.vcf" | grep -c "CLNSIG=" || true)
echo "  Total PASS variants called : ${N_TOTAL}"
echo "  Matched a ClinVar record   : ${N_ANNOT}"
echo ""
echo "Annotated VCF: ${OUTDIR}/vcf/${SAMPLE}.clinvar_annotated.vcf"
echo "Next: python 02_generate_report.py ${OUTDIR}/vcf/${SAMPLE}.clinvar_annotated.vcf ${SAMPLE} ${OUTDIR}/report"

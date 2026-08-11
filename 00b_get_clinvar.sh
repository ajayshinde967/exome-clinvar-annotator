#!/usr/bin/env bash
# ==============================================================================
# 00b_get_clinvar.sh
# Downloads ClinVar and produces a VCF whose chromosome naming is auto-matched
# to your reference FASTA. Works for any exome (whole reference, not just chr21).
#
# Usage:
#   ./00b_get_clinvar.sh <reference.fa> <output_dir> [genome_build]
#   genome_build: GRCh38 (default) or GRCh37
# ==============================================================================
set -euo pipefail

REF="${1:?Usage: $0 <reference.fa> <output_dir> [GRCh38|GRCh37]}"
OUTDIR="${2:?Usage: $0 <reference.fa> <output_dir> [GRCh38|GRCh37]}"
BUILD="${3:-GRCh38}"

mkdir -p "${OUTDIR}"

if [[ "${BUILD}" == "GRCh38" ]]; then
  URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
elif [[ "${BUILD}" == "GRCh37" ]]; then
  URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
else
  echo "ERROR: genome_build must be GRCh38 or GRCh37"; exit 1
fi

echo "[1/4] Downloading ClinVar (${BUILD})"
curl -s "${URL}" -o "${OUTDIR}/clinvar_full.vcf.gz"
curl -s "${URL}.tbi" -o "${OUTDIR}/clinvar_full.vcf.gz.tbi"

echo "[2/4] Detecting chromosome naming style in your reference"
REF_CONTIGS=$(grep "^>" "${REF}" | sed 's/^>//' | awk '{print $1}')
FIRST_CONTIG=$(echo "${REF_CONTIGS}" | head -1)
if [[ "${FIRST_CONTIG}" == chr* ]]; then
  STYLE="chr"
else
  STYLE="plain"
fi
echo "  Reference naming style: ${STYLE} (e.g. '${FIRST_CONTIG}')"

echo "[3/4] Renaming ClinVar contigs to match, if needed (NCBI ClinVar uses plain '1','2'...'X','Y','MT')"
if [[ "${STYLE}" == "chr" ]]; then
  {
    for c in $(seq 1 22) X Y MT; do echo "${c} chr${c}"; done
  } > "${OUTDIR}/rename_map.txt"
  # MT in UCSC-style refs is usually chrM, not chrMT - handle common case
  sed -i 's/MT chrMT/MT chrM/' "${OUTDIR}/rename_map.txt"
  bcftools annotate --rename-chrs "${OUTDIR}/rename_map.txt" \
    "${OUTDIR}/clinvar_full.vcf.gz" -Oz -o "${OUTDIR}/clinvar_matched.vcf.gz"
else
  cp "${OUTDIR}/clinvar_full.vcf.gz" "${OUTDIR}/clinvar_matched.vcf.gz"
fi
tabix -p vcf "${OUTDIR}/clinvar_matched.vcf.gz"

echo "[4/4] Restricting ClinVar to contigs actually present in your reference"
echo "${REF_CONTIGS}" > "${OUTDIR}/ref_contigs.txt"
bcftools view -R <(awk '{print $1"\t0\t500000000"}' "${OUTDIR}/ref_contigs.txt") \
  "${OUTDIR}/clinvar_matched.vcf.gz" -Oz -o "${OUTDIR}/clinvar_exome.vcf.gz" || \
  cp "${OUTDIR}/clinvar_matched.vcf.gz" "${OUTDIR}/clinvar_exome.vcf.gz"
tabix -p vcf "${OUTDIR}/clinvar_exome.vcf.gz"

rm -f "${OUTDIR}/clinvar_full.vcf.gz"* "${OUTDIR}/clinvar_matched.vcf.gz"*

echo ""
echo "Done: ${OUTDIR}/clinvar_exome.vcf.gz"

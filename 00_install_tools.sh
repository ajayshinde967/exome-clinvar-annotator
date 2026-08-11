#!/usr/bin/env bash
# ==============================================================================
# 00_install_tools.sh
# Installs every tool used by the WES → ClinVar pipeline into one conda env.
# Versions pinned to known-compatible releases (Aug 2026).
# ==============================================================================
set -euo pipefail

ENV_NAME="${1:-wes-pipeline}"

if ! command -v mamba &>/dev/null && ! command -v conda &>/dev/null; then
  echo "conda/mamba not found. Installing Miniforge (conda + mamba)..."
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    -O /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p "${HOME}/miniforge3"
  export PATH="${HOME}/miniforge3/bin:${PATH}"
  echo "Add this to your ~/.bashrc: export PATH=\"${HOME}/miniforge3/bin:\$PATH\""
fi

CONDA_BIN="mamba"
command -v mamba &>/dev/null || CONDA_BIN="conda"

echo "Creating environment '${ENV_NAME}' with pinned tool versions..."
${CONDA_BIN} create -n "${ENV_NAME}" -y \
  -c bioconda -c conda-forge \
  fastp=0.23.4 \
  bwa=0.7.17 \
  samtools=1.19 \
  bcftools=1.19 \
  tabix \
  gatk4=4.5.0.0 \
  snpsift=5.2 \
  snpeff=5.2 \
  seqtk \
  bedtools=2.31.1

echo ""
echo "Done. Activate with:"
echo "  ${CONDA_BIN} activate ${ENV_NAME}"
echo ""
echo "Then verify:"
echo "  fastp --version"
echo "  bwa 2>&1 | head -3"
echo "  samtools --version | head -1"
echo "  bcftools --version | head -1"
echo "  gatk --version"
echo "  SnpSift 2>&1 | head -5"
echo "  bedtools --version"

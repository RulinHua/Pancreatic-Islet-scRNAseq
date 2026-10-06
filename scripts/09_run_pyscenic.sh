#!/usr/bin/env bash
set -euo pipefail

# Activate the conda environment containing pySCENIC before running:
#   conda activate SCENIC
#
# Required:
#   export MOTIF_DIR=/path/to/cisTarget_databases_hg38
#
# Optional:
#   export SCENIC_DIR=/path/to/output
#   export NWORKERS=40

SCENIC_DIR="${SCENIC_DIR:-$(pwd)/results/SCENIC}"
MOTIF_DIR="${MOTIF_DIR:?Please set MOTIF_DIR to the hg38 cisTarget database directory}"
NWORKERS="${NWORKERS:-40}"

mkdir -p "${SCENIC_DIR}"

TF_LIST="${MOTIF_DIR}/hs_hgnc_curated_tfs.txt"
ANNOT_TBL="${MOTIF_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
FEATHERS=( "${MOTIF_DIR}"/*mc9nr*.feather )

# GRNBoost2 co-expression network
pyscenic grn \
  "${SCENIC_DIR}/HPAP_V3_GRN.loom" \
  "${TF_LIST}" \
  --method grnboost2 \
  --seed 123 \
  --num_workers "${NWORKERS}" \
  --sparse \
  --cell_id_attribute CellID \
  --gene_attribute Gene \
  --output "${SCENIC_DIR}/HPAP_V3_adjacencies.tsv" \
  2>&1 | tee "${SCENIC_DIR}/01_HPAP_V3_GRN.log"

# cisTarget motif enrichment/pruning
pyscenic ctx \
  "${SCENIC_DIR}/HPAP_V3_adjacencies.tsv" \
  "${FEATHERS[@]}" \
  --annotations_fname "${ANNOT_TBL}" \
  --expression_mtx_fname "${SCENIC_DIR}/HPAP_V3_GRN.loom" \
  --cell_id_attribute CellID \
  --gene_attribute Gene \
  --mode custom_multiprocessing \
  --num_workers "${NWORKERS}" \
  --output "${SCENIC_DIR}/HPAP_V3_regulons.csv" \
  2>&1 | tee "${SCENIC_DIR}/02_HPAP_V3_CTX.log"

# AUCell in all HPAP V3 beta cells
pyscenic aucell \
  "${SCENIC_DIR}/HPAP_V3_ALL.loom" \
  "${SCENIC_DIR}/HPAP_V3_regulons.csv" \
  --cell_id_attribute CellID \
  --gene_attribute Gene \
  --seed 123 \
  --num_workers "${NWORKERS}" \
  --output "${SCENIC_DIR}/HPAP_V3_AUC.csv" \
  2>&1 | tee "${SCENIC_DIR}/03_HPAP_V3_AUC.log"

# Apply the same regulons to IIDP V3 beta cells
pyscenic aucell \
  "${SCENIC_DIR}/IIDP_V3_ALL.loom" \
  "${SCENIC_DIR}/HPAP_V3_regulons.csv" \
  --cell_id_attribute CellID \
  --gene_attribute Gene \
  --seed 123 \
  --num_workers "${NWORKERS}" \
  --output "${SCENIC_DIR}/IIDP_V3_AUC.csv" \
  2>&1 | tee "${SCENIC_DIR}/04_IIDP_V3_AUC.log"

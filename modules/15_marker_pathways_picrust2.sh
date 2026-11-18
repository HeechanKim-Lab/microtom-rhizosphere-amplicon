#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 15_marker_pathways_picrust2.sh
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1-3, GJ1-2; N=110 Locked Cohort)
# Platform: macOS (conda environment picrust2)
# Description: Infers MetaCyc pathways separately for each RF marker ASV
#              (marker-only tables from module 14) with PICRUSt2
#              pathway_pipeline.py in per-sequence mode, and adds pathway
#              names. Run after module 14 and before module 16.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

THREADS=4
OUT_DIR="rf_picrust2"

echo "=== [$(date +'%T')] Initializing PICRUSt2 Pathway Inference for RF Marker ASVs ==="

# Initialize Conda for subshell scripting (as in pipeline/run_pipeline.sh)
eval "$(conda shell.bash hook)"
conda activate picrust2

# Verify core inputs
for f in picrust2_out_16S/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
         "${OUT_DIR}/seqtab_norm_markers.tsv.gz" \
         "${OUT_DIR}/EC_predicted_markers.tsv.gz"; do
  if [[ ! -f "$f" ]]; then
    echo "[-] ERROR: Required file '$f' not found (run module 14 first)." >&2
    exit 1
  fi
done

# ==============================================================================
# 1. Pathways per marker ASV (MinPath run on each marker ASV's predicted genome)
# ==============================================================================
echo "--- Step 1: pathway_pipeline.py, per-sequence mode ---"

pathway_pipeline.py \
  -i picrust2_out_16S/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
  -o "${OUT_DIR}/pathways_markers" \
  --per_sequence_contrib \
  --per_sequence_abun "${OUT_DIR}/seqtab_norm_markers.tsv.gz" \
  --per_sequence_function "${OUT_DIR}/EC_predicted_markers.tsv.gz" \
  --coverage \
  -p "${THREADS}"

# Main output: path_abun_predictions.tsv.gz
#   one row per marker ASV, one column per pathway; > 0 = pathway predicted present

# ==============================================================================
# 2. MetaCyc pathway names
# ==============================================================================
echo "--- Step 2: add_descriptions.py ---"

add_descriptions.py \
  -i "${OUT_DIR}/pathways_markers/path_abun_unstrat_per_seq.tsv.gz" \
  -m METACYC \
  -o "${OUT_DIR}/pathways_markers/path_abun_unstrat_per_seq_descrip.tsv.gz"

echo "=== [$(date +'%T')] PICRUSt2 Pathway Inference for RF Marker ASVs Complete ==="

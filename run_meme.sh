#!/usr/bin/env bash

set -u
set -o pipefail

# =============================================================================
# MEME ANALYSIS FOR SELECTED GENES
# Linux bash version
# - selected genes only
# - max 8 parallel files
# - exactly CPU=1 per HyPhy process
# - analysis is run on Background branches
# =============================================================================

# =========================
# SETTINGS
# =========================

HYPHY_BIN="/home/poroshina/.conda/envs/based/bin/hyphy"

MAX_JOBS=8
HYPHY_CPU_PER_JOB=1

# 0 = skip existing non-empty json
# 1 = rerun even if json already exists
FORCE_RERUN_MEME_BACKGROUND=0

PAL2NAL_DIR="pipeline_results/05_pal2nal_normalized_no_duplicates"
TREE_DIR="pipeline_results/07_labeled_trees"
MEME_DIR="pipeline_results/09_meme_background"

ALIGNMENT_SUFFIX="_bird_only_nucleotide_codon_aligned.fasta"
TREE_SUFFIX="_bird_only_nucleotide_labeled.nwk"

GENES=(
  "GALNT7"
  "OR10C1"
  "ZNF532"
  "SELL"
  "CA14"
  "HLA-DMA"
  "EPHA6"
  "RAB18"
  "MAGI3"
)

mkdir -p "$MEME_DIR"

# =========================
# THREAD SAFETY LIMITS
# =========================

# Дополнительная защита: не даем OpenMP / BLAS распараллеливаться сверх CPU=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# =========================
# CHECK HYPHY
# =========================

if [[ ! -x "$HYPHY_BIN" ]]; then
  echo "[ERROR] HyPhy not found or not executable: $HYPHY_BIN"
  echo "Check with:"
  echo "  which hyphy"
  echo "  realpath \"\$(which hyphy)\""
  exit 1
fi

echo "Using HyPhy: $HYPHY_BIN"
echo "MAX_JOBS: $MAX_JOBS"
echo "HYPHY_CPU_PER_JOB: $HYPHY_CPU_PER_JOB"
echo "TOTAL CPU LIMIT: $((MAX_JOBS * HYPHY_CPU_PER_JOB))"
echo "MEME branch set: Background"
echo "Output dir: $MEME_DIR"
echo

# =========================
# FUNCTIONS
# =========================

find_alignment() {
  local gene="$1"
  local exact="${PAL2NAL_DIR}/${gene}${ALIGNMENT_SUFFIX}"

  if [[ -s "$exact" ]]; then
    echo "$exact"
    return 0
  fi

  local candidate
  candidate=$(find "$PAL2NAL_DIR" -maxdepth 1 -type f -name "${gene}*${ALIGNMENT_SUFFIX}" | sort | head -n 1)

  if [[ -n "${candidate:-}" && -s "$candidate" ]]; then
    echo "$candidate"
    return 0
  fi

  return 1
}

find_tree() {
  local gene="$1"
  local exact="${TREE_DIR}/${gene}${TREE_SUFFIX}"

  if [[ -s "$exact" ]]; then
    echo "$exact"
    return 0
  fi

  local candidate
  candidate=$(find "$TREE_DIR" -maxdepth 1 -type f -name "${gene}*${TREE_SUFFIX}" | sort | head -n 1)

  if [[ -n "${candidate:-}" && -s "$candidate" ]]; then
    echo "$candidate"
    return 0
  fi

  return 1
}

run_meme_background_one_gene() {
  local gene="$1"
  local alignment="$2"
  local tree="$3"
  local output_json="$4"
  local output_log="$5"

  echo "[START BACKGROUND] $gene"

  {
    echo "GENE: $gene"
    echo "START: $(date)"
    echo "ALIGNMENT: $alignment"
    echo "TREE: $tree"
    echo "BRANCH SET: Background"
    echo "OUTPUT: $output_json"
    echo
    echo "COMMAND:"
    echo "$HYPHY_BIN CPU=${HYPHY_CPU_PER_JOB} ENV=TOLERATE_NUMERICAL_ERRORS=1 meme --alignment \"$alignment\" --tree \"$tree\" --branches Background --output \"$output_json\""
    echo
  } > "$output_log"

  "$HYPHY_BIN" \
    "CPU=${HYPHY_CPU_PER_JOB}" \
    "ENV=TOLERATE_NUMERICAL_ERRORS=1" \
    meme \
    --alignment "$alignment" \
    --tree "$tree" \
    --branches Background \
    --output "$output_json" \
    >> "$output_log" 2>&1

  local rc=$?

  {
    echo
    echo "END: $(date)"
    echo "RETURNCODE: $rc"
  } >> "$output_log"

  if [[ "$rc" -eq 0 && -s "$output_json" ]]; then
    echo "[OK BACKGROUND] $gene -> $output_json"
  else
    echo "[FAILED BACKGROUND] $gene | returncode=$rc | log=$output_log"
  fi

  return "$rc"
}

wait_for_slot() {
  while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do
    wait -n
  done
}

# =========================
# MAIN LOOP
# =========================

tasks_started=0
tasks_skipped=0
tasks_missing=0

for gene in "${GENES[@]}"; do
  alignment=""
  tree=""

  if ! alignment=$(find_alignment "$gene"); then
    echo "[MISSING ALIGNMENT] $gene"
    tasks_missing=$((tasks_missing + 1))
    continue
  fi

  if ! tree=$(find_tree "$gene"); then
    echo "[MISSING TREE] $gene"
    echo "  alignment found: $alignment"
    tasks_missing=$((tasks_missing + 1))
    continue
  fi

  output_json="${MEME_DIR}/${gene}_meme_background.json"
  output_log="${MEME_DIR}/${gene}_meme_background.log"

  if [[ "$FORCE_RERUN_MEME_BACKGROUND" -eq 0 && -s "$output_json" ]]; then
    echo "[SKIP READY BACKGROUND] $gene -> $output_json"
    tasks_skipped=$((tasks_skipped + 1))
    continue
  fi

  wait_for_slot

  run_meme_background_one_gene "$gene" "$alignment" "$tree" "$output_json" "$output_log" &
  tasks_started=$((tasks_started + 1))
done

wait

echo
echo "MEME background batch finished."
echo "Started: $tasks_started"
echo "Skipped ready: $tasks_skipped"
echo "Missing input: $tasks_missing"
echo "Output dir: $MEME_DIR"

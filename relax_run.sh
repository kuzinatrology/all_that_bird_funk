#!/usr/bin/env bash

set -u
set -o pipefail

RELAX_WORKERS=8

HYPHY_BIN="${HYPHY_BIN:-/home/poroshina/.conda/envs/based/bin/hyphy}"

FORCE_RERUN_RELAX=0

PAL2NAL_DIR="pipeline_results/05_pal2nal_normalized_no_duplicates"
LABELED_TREES_DIR="pipeline_results/07_labeled_trees"
RELAX_DIR="pipeline_results/08b_relax"
SUMMARY_DIR="pipeline_results/11_summary"

TEST_LABEL="Foreground"
REFERENCE_LABEL="Background"

RUN_LIST_DIR="${RELAX_DIR}/_run_lists"

MASTER_LOG="${RELAX_DIR}/run_relax.master.log"

COMPLETED_LIST="${RUN_LIST_DIR}/completed.tsv"
SKIPPED_LIST="${RUN_LIST_DIR}/skipped.tsv"
FAILED_LIST="${RUN_LIST_DIR}/failed.tsv"
MISSING_TREES_LIST="${RUN_LIST_DIR}/missing_trees.tsv"
BAD_INPUT_LIST="${RUN_LIST_DIR}/bad_inputs.tsv"

SUMMARY_CSV="${SUMMARY_DIR}/relax_results.csv"

mkdir -p "$RELAX_DIR" "$RUN_LIST_DIR" "$SUMMARY_DIR"

: > "$COMPLETED_LIST"
: > "$SKIPPED_LIST"
: > "$FAILED_LIST"
: > "$MISSING_TREES_LIST"
: > "$BAD_INPUT_LIST"
: > "$MASTER_LOG"

log_msg() {
    local msg="$1"
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$msg" | tee -a "$MASTER_LOG"
}

run_relax_one() {
    local codon_path="$1"
    local fname stem tree_fname tree_path output_json output_log status

    fname="$(basename "$codon_path")"
    stem="${fname%.*}"

    tree_fname="${fname/_codon_aligned.fasta/_labeled.nwk}"
    tree_path="${LABELED_TREES_DIR}/${tree_fname}"

    output_json="${RELAX_DIR}/${stem}_relax.json"
    output_log="${RELAX_DIR}/${stem}_relax.log"

    if [[ ! -s "$codon_path" ]]; then
        printf '%s\t%s\n' "$fname" "empty_alignment" >> "$BAD_INPUT_LIST"
        return 0
    fi

    if [[ "$FORCE_RERUN_RELAX" != "1" && -s "$output_json" ]]; then
        printf '%s\t%s\n' "$fname" "$output_json" >> "$SKIPPED_LIST"
        return 0
    fi

    if [[ ! -s "$tree_path" ]]; then
        printf '%s\t%s\n' "$fname" "tree_not_found:${tree_path}" >> "$MISSING_TREES_LIST"
        return 0
    fi

    {
        echo "===================================================================="
        echo "RELAX job"
        echo "alignment: $codon_path"
        echo "tree:      $tree_path"
        echo "test:      $TEST_LABEL"
        echo "reference: $REFERENCE_LABEL"
        echo "output:    $output_json"
        echo "started:   $(date '+%Y-%m-%d %H:%M:%S')"
        echo "===================================================================="
        echo

        if ! grep -q "$TEST_LABEL" "$tree_path"; then
            echo "WARNING: test label '${TEST_LABEL}' was not found literally in the tree file."
        fi

        if ! grep -q "$REFERENCE_LABEL" "$tree_path"; then
            echo "WARNING: reference label '${REFERENCE_LABEL}' was not found literally in the tree file."
            echo "RELAX may fail if this branch set is required but absent."
        fi

        echo
    } > "$output_log"

    "$HYPHY_BIN" relax \
        --alignment "$codon_path" \
        --tree "$tree_path" \
        --test "$TEST_LABEL" \
        --reference "$REFERENCE_LABEL" \
        --output "$output_json" \
        >> "$output_log" 2>&1

    status=$?

    {
        echo
        echo "finished: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "exit_code: $status"
    } >> "$output_log"

    if [[ "$status" -eq 0 && -s "$output_json" ]]; then
        printf '%s\t%s\t%s\n' "$fname" "$tree_path" "$output_json" >> "$COMPLETED_LIST"
    else
        printf '%s\t%s\t%s\t%s\n' "$fname" "$tree_path" "$output_json" "exit_code_${status}" >> "$FAILED_LIST"
    fi

    return 0
}

if [[ ! -x "$HYPHY_BIN" ]]; then
    echo "ERROR: HyPhy binary is not executable: $HYPHY_BIN" >&2
    echo "HYPHY_BIN=/path/to/hyphy bash run_relax.sh" >&2
    exit 1
fi

if [[ ! -d "$PAL2NAL_DIR" ]]; then
    echo "ERROR: PAL2NAL_DIR does not exist: $PAL2NAL_DIR" >&2
    exit 1
fi

if [[ ! -d "$LABELED_TREES_DIR" ]]; then
    echo "ERROR: LABELED_TREES_DIR does not exist: $LABELED_TREES_DIR" >&2
    exit 1
fi

log_msg "Starting RELAX"
log_msg "HyPhy: $HYPHY_BIN"
log_msg "Alignments: $PAL2NAL_DIR"
log_msg "Trees: $LABELED_TREES_DIR"
log_msg "Output: $RELAX_DIR"
log_msg "Parallel file jobs: $RELAX_WORKERS"
log_msg "Labels: test=${TEST_LABEL}, reference=${REFERENCE_LABEL}"
log_msg "Force rerun: $FORCE_RERUN_RELAX"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

job_count=0

while IFS= read -r -d '' codon_path; do
    while [[ "$(jobs -rp | wc -l)" -ge "$RELAX_WORKERS" ]]; do
        sleep 2
    done

    log_msg "Queued $(basename "$codon_path")"

    run_relax_one "$codon_path" &

    job_count=$((job_count + 1))

done < <(
    find "$PAL2NAL_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*_codon_aligned.fasta' \
        -print0 \
    | sort -z
)

wait || true

log_msg "All RELAX jobs finished"
log_msg "Queued files: $job_count"
log_msg "Completed: $(wc -l < "$COMPLETED_LIST")"
log_msg "Skipped existing: $(wc -l < "$SKIPPED_LIST")"
log_msg "Missing trees: $(wc -l < "$MISSING_TREES_LIST")"
log_msg "Bad inputs: $(wc -l < "$BAD_INPUT_LIST")"
log_msg "Failed: $(wc -l < "$FAILED_LIST")"

python3 - <<'PY'
import csv
import json
import math
from pathlib import Path

relax_dir = Path("pipeline_results/08b_relax")
summary_dir = Path("pipeline_results/11_summary")
summary_dir.mkdir(parents=True, exist_ok=True)

out_csv = summary_dir / "relax_results.csv"

rows = []

for json_path in sorted(relax_dir.glob("*_relax.json")):
    if json_path.stat().st_size == 0:
        continue

    stem = json_path.name.replace("_relax.json", "")

    gene = stem
    for suffix in [
        "_bird_only_nucleotide_codon_aligned",
        "_nucleotide_codon_aligned",
        "_codon_aligned",
    ]:
        if gene.endswith(suffix):
            gene = gene[: -len(suffix)]
            break

    try:
        with json_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        rows.append({
            "gene": gene,
            "json_file": str(json_path),
            "K": "",
            "p_value": "",
            "interpretation": "parse_error",
            "status": f"parse_error:{type(e).__name__}",
        })
        continue

    test_results = data.get("test results", {}) or {}

    k_value = test_results.get("relaxation or intensification parameter", None)
    p_value = test_results.get("p-value", None)

    interpretation = "unknown"

    try:
        k_float = float(k_value)
        if math.isfinite(k_float):
            if k_float < 1:
                interpretation = "relaxation"
            elif k_float > 1:
                interpretation = "intensification"
            else:
                interpretation = "neutral"
    except Exception:
        pass

    rows.append({
        "gene": gene,
        "json_file": str(json_path),
        "K": k_value if k_value is not None else "",
        "p_value": p_value if p_value is not None else "",
        "interpretation": interpretation,
        "status": "ok" if k_value is not None and p_value is not None else "missing_test_results",
    })

with out_csv.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "gene",
            "json_file",
            "K",
            "p_value",
            "interpretation",
            "status",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Saved summary: {out_csv}")
print(f"Parsed RELAX JSON files: {len(rows)}")
PY

log_msg "Summary CSV: $SUMMARY_CSV"
log_msg "Done"

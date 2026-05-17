#!/usr/bin/env bash

set -u
set -o pipefail

HYPHY_BIN="/home/poroshina/.conda/envs/based/bin/hyphy"
FITMG94_BF="/home/poroshina/dn_ds_pipeline_easy_search/hyphy-analyses/FitMG94/FitMG94.bf"

FORCE_RERUN_FITMG94=0
WORKERS=4

GENES_FILE="fitmg94_genes_2.csv"

PAL2NAL_DIR="pipeline_results/05_pal2nal_normalized_no_duplicates"
TREE_DIR="pipeline_results/07_labeled_trees"
FITMG94_DIR="pipeline_results/09_fitmg94"

ALIGNMENT_SUFFIX="_bird_only_nucleotide_codon_aligned.fasta"
TREE_SUFFIX="_bird_only_nucleotide_labeled.nwk"

TASKS_FILE="${FITMG94_DIR}/fitmg94_tasks.tsv"

mkdir -p "$FITMG94_DIR"

> "$TASKS_FILE"

tail -n +2 "$GENES_FILE" | while IFS= read -r gene_name; do
    gene_name="${gene_name//$'\r'/}"
    gene_name="${gene_name//\"/}"

    if [[ -z "$gene_name" ]]; then
        continue
    fi

    codon_path="${PAL2NAL_DIR}/${gene_name}${ALIGNMENT_SUFFIX}"
    tree_path="${TREE_DIR}/${gene_name}${TREE_SUFFIX}"

    output_json="${FITMG94_DIR}/${gene_name}_bird_only_nucleotide_codon_aligned_fitmg94.json"
    output_log="${FITMG94_DIR}/${gene_name}_bird_only_nucleotide_codon_aligned_fitmg94.log"

    if [[ "$FORCE_RERUN_FITMG94" == "0" && -s "$output_json" ]]; then
        continue
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$gene_name" \
        "$codon_path" \
        "$tree_path" \
        "$output_json" \
        "$output_log" >> "$TASKS_FILE"
done

echo "Prepared FitMG94 tasks: $(wc -l < "$TASKS_FILE")"
echo "Running FitMG94 jobs"
echo "Workers: $WORKERS"
echo "FORCE_RERUN_FITMG94=$FORCE_RERUN_FITMG94"

export HYPHY_BIN
export FITMG94_BF
export OMP_NUM_THREADS=1

run_one() {
    gene_name="$1"
    codon_path="$2"
    tree_path="$3"
    output_json="$4"
    output_log="$5"

    echo "[RUN] $gene_name"

    "$HYPHY_BIN" "$FITMG94_BF" \
        --alignment "$codon_path" \
        --tree "$tree_path" \
        --type lineage \
        --kill-zero-lengths No \
        --output "$output_json" \
        > "$output_log" 2>&1

    if [[ -s "$output_json" ]]; then
        echo "[OK] $gene_name"
    else
        echo "[FAIL] $gene_name | log: $output_log"
    fi
}

export -f run_one

cat "$TASKS_FILE" | xargs -r -P "$WORKERS" -n 5 bash -c 'run_one "$@"' _

echo "FitMG94 run finished"

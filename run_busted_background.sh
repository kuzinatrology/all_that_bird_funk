#!/usr/bin/env bash

set -u
set -o pipefail

HYPHY_BIN="/home/poroshina/.conda/envs/based/bin/hyphy"

BUSTED_WORKERS="${BUSTED_WORKERS:-4}"
FORCE_RERUN="${FORCE_RERUN:-0}"

BUSTED_CSV="pipeline_results/11_summary/busted_results.csv"

PAL2NAL_DIR="pipeline_results/05_pal2nal_normalized_no_duplicates"
LABELED_TREES_DIR="pipeline_results/07_labeled_trees"
BUSTED_BACKGROUND_DIR="pipeline_results/08c_busted_background"

FILE_SUFFIX="_bird_only_nucleotide_codon_aligned.fasta"
TREE_SUFFIX="_bird_only_nucleotide_labeled.nwk"

mkdir -p "$BUSTED_BACKGROUND_DIR"

MAIN_LOG="$BUSTED_BACKGROUND_DIR/_busted_background_run.log"
TASKS_TSV="$BUSTED_BACKGROUND_DIR/_busted_background_tasks.tsv"
RESULTS_TSV="$BUSTED_BACKGROUND_DIR/_busted_background_results.tsv"
GENE_LIST="$BUSTED_BACKGROUND_DIR/_genes_p001_from_busted_results.txt"

: > "$MAIN_LOG"
printf "gene\tcodon_file\tcodon_path\ttree_path\toutput_json\toutput_log\n" > "$TASKS_TSV"
printf "gene\tcodon_file\tstatus\treturncode\toutput_json\toutput_log\n" > "$RESULTS_TSV"

log() {
    printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$MAIN_LOG"
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ ! -x "$HYPHY_BIN" ]]; then
    log "ERROR: HyPhy binary not found or not executable: $HYPHY_BIN"
    exit 1
fi

if [[ ! -f "$BUSTED_CSV" ]]; then
    log "ERROR: CSV not found: $BUSTED_CSV"
    exit 1
fi

if [[ ! -d "$PAL2NAL_DIR" ]]; then
    log "ERROR: alignment directory not found: $PAL2NAL_DIR"
    exit 1
fi

if [[ ! -d "$LABELED_TREES_DIR" ]]; then
    log "ERROR: labeled trees directory not found: $LABELED_TREES_DIR"
    exit 1
fi

SELECTED_N=$(
python3 - "$BUSTED_CSV" "$GENE_LIST" <<'PY'
import csv
import sys

csv_path = sys.argv[1]
out_path = sys.argv[2]

genes = []
seen = set()

with open(csv_path, newline="") as f:
    reader = csv.DictReader(f)

    if "gene" not in reader.fieldnames or "p_value" not in reader.fieldnames:
        raise SystemExit("CSV must contain columns: gene, p_value")

    for row in reader:
        gene = (row.get("gene") or "").strip()
        p_raw = (row.get("p_value") or "").strip()

        if not gene or not p_raw:
            continue

        try:
            p_value = float(p_raw)
        except ValueError:
            continue

        if p_value < 0.01 and gene not in seen:
            seen.add(gene)
            genes.append(gene)

with open(out_path, "w") as out:
    for gene in genes:
        out.write(gene + "\n")

print(len(genes))
PY
)

log "CSV: $BUSTED_CSV"
log "Genes with p_value < 0.01 in CSV: $SELECTED_N"
log "OUT_DIR: $BUSTED_BACKGROUND_DIR"
log "BUSTED_WORKERS: $BUSTED_WORKERS"
log "FORCE_RERUN: $FORCE_RERUN"

declare -A GENE_SET

while IFS= read -r gene; do
    [[ -n "$gene" ]] && GENE_SET["$gene"]=1
done < "$GENE_LIST"

scanned_alignments=0
tasks_n=0
skipped_existing=0
missing_tree=0
empty_alignment=0

for codon_path in "$PAL2NAL_DIR"/*"$FILE_SUFFIX"; do
    [[ -e "$codon_path" ]] || continue

    codon_file="$(basename "$codon_path")"
    gene_name="${codon_file%"$FILE_SUFFIX"}"

    ((scanned_alignments++))

    if [[ -z "${GENE_SET[$gene_name]+x}" ]]; then
        continue
    fi

    if [[ ! -s "$codon_path" ]]; then
        ((empty_alignment++))
        log "SKIP empty alignment: $codon_file"
        continue
    fi

    tree_file="${gene_name}${TREE_SUFFIX}"
    tree_path="$LABELED_TREES_DIR/$tree_file"

    if [[ ! -f "$tree_path" ]]; then
        ((missing_tree++))
        log "SKIP tree not found for $codon_file: $tree_path"
        continue
    fi

    output_json="$BUSTED_BACKGROUND_DIR/${codon_file}_busted_background.json"
    output_log="$BUSTED_BACKGROUND_DIR/${codon_file}_busted_background.log"

    if [[ "$FORCE_RERUN" == "0" && -s "$output_json" ]]; then
        ((skipped_existing++))
        log "SKIP existing result: $gene_name | $output_json"
        continue
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$gene_name" \
        "$codon_file" \
        "$codon_path" \
        "$tree_path" \
        "$output_json" \
        "$output_log" >> "$TASKS_TSV"

    ((tasks_n++))
done

log "Scanned alignments: $scanned_alignments"
log "Skipped existing results: $skipped_existing"
log "Skipped missing trees: $missing_tree"
log "Skipped empty alignments: $empty_alignment"
log "Found $tasks_n genes for background BUSTED-E analysis"

run_background_busted() {
    local gene="$1"
    local codon_file="$2"
    local codon_path="$3"
    local tree_path="$4"
    local output_json="$5"
    local output_log="$6"

    log "START: $gene"

    "$HYPHY_BIN" CPU=1 busted \
        --alignment "$codon_path" \
        --tree "$tree_path" \
        --branches Background \
        --error-sink Yes \
        --output "$output_json" \
        > "$output_log" 2>&1

    local rc=$?

    if [[ "$rc" -eq 0 && -s "$output_json" ]]; then
        log "DONE: $gene"
        printf "%s\t%s\tsuccess\t%s\t%s\t%s\n" \
            "$gene" "$codon_file" "$rc" "$output_json" "$output_log" >> "$RESULTS_TSV"
    else
        log "FAILED: $gene | returncode=$rc | log=$output_log"
        printf "%s\t%s\tfailed\t%s\t%s\t%s\n" \
            "$gene" "$codon_file" "$rc" "$output_json" "$output_log" >> "$RESULTS_TSV"
    fi
}

if (( tasks_n > 0 )); then
    first_line=1

    while IFS=$'\t' read -r gene codon_file codon_path tree_path output_json output_log; do
        if (( first_line )); then
            first_line=0
            continue
        fi

        while (( $(jobs -rp | wc -l) >= BUSTED_WORKERS )); do
            wait -n || true
        done

        run_background_busted \
            "$gene" \
            "$codon_file" \
            "$codon_path" \
            "$tree_path" \
            "$output_json" \
            "$output_log" &
    done < "$TASKS_TSV"

    wait
fi

success_n=$(awk -F'\t' 'NR > 1 && $3 == "success" {c++} END {print c+0}' "$RESULTS_TSV")
failed_n=$(awk -F'\t' 'NR > 1 && $3 == "failed" {c++} END {print c+0}' "$RESULTS_TSV")

log "Completed $success_n / $tasks_n background BUSTED-E analyses"
log "Failed: $failed_n"
log "Main log: $MAIN_LOG"
log "Results table: $RESULTS_TSV"

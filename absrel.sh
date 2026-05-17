#!/usr/bin/env bash

set -u
set -o pipefail

HYPHY_BIN="/home/poroshina/.conda/envs/based/bin/hyphy"
ABSREL_WORKERS=4
FORCE_RERUN=0

BUSTED_CSV="pipeline_results/11_summary/busted_results.csv"

PAL2NAL_DIR="pipeline_results/05_pal2nal_normalized_no_duplicates"
LABELED_TREES_DIR="pipeline_results/07_labeled_trees"
ABSREL_DIR="pipeline_results/09_absrel"

FILE_SUFFIX="_bird_only_nucleotide_codon_aligned.fasta"
TREE_SUFFIX="_bird_only_nucleotide_labeled.nwk"

mkdir -p "$ABSREL_DIR"

MAIN_LOG="$ABSREL_DIR/_absrel_run.log"
RESULTS_TSV="$ABSREL_DIR/_absrel_results.tsv"
GENE_LIST="$ABSREL_DIR/_absrel_genes.txt"

exec > "$MAIN_LOG" 2>&1

echo "Started aBSREL run: $(date)"
echo "HYPHY_BIN=$HYPHY_BIN"
echo "ABSREL_WORKERS=$ABSREL_WORKERS"
echo "FORCE_RERUN=$FORCE_RERUN"
echo "BUSTED_CSV=$BUSTED_CSV"
echo "PAL2NAL_DIR=$PAL2NAL_DIR"
echo "LABELED_TREES_DIR=$LABELED_TREES_DIR"
echo "ABSREL_DIR=$ABSREL_DIR"
echo

if [[ ! -x "$HYPHY_BIN" ]]; then
    echo "ERROR: HYPHY_BIN is not executable: $HYPHY_BIN"
    exit 1
fi

if [[ ! -f "$BUSTED_CSV" ]]; then
    echo "ERROR: BUSTED CSV not found: $BUSTED_CSV"
    exit 1
fi

if [[ ! -d "$PAL2NAL_DIR" ]]; then
    echo "ERROR: PAL2NAL_DIR not found: $PAL2NAL_DIR"
    exit 1
fi

if [[ ! -d "$LABELED_TREES_DIR" ]]; then
    echo "ERROR: LABELED_TREES_DIR not found: $LABELED_TREES_DIR"
    exit 1
fi

python3 - "$BUSTED_CSV" "$GENE_LIST" <<'PY'
import csv
import sys

csv_path = sys.argv[1]
out_path = sys.argv[2]

with open(csv_path, newline="", encoding="utf-8-sig") as f:
    sample = f.read(4096)
    f.seek(0)

    try:
        dialect = csv.Sniffer().sniff(sample)
    except csv.Error:
        dialect = csv.excel

    reader = csv.DictReader(f, dialect=dialect)

    if reader.fieldnames is None:
        raise SystemExit("ERROR: empty CSV or no header found")

    fieldnames = [name.strip() for name in reader.fieldnames]

    if "gene" not in fieldnames:
        raise SystemExit(f"ERROR: column 'gene' not found in CSV. Found columns: {fieldnames}")

    if "p_value" not in fieldnames:
        raise SystemExit(f"ERROR: column 'p_value' not found in CSV. Found columns: {fieldnames}")

    genes = []
    seen = set()

    for row in reader:
        row = {str(k).strip(): v for k, v in row.items()}

        gene = str(row.get("gene", "")).strip()
        p_raw = str(row.get("p_value", "")).strip()

        if not gene or not p_raw:
            continue

        try:
            p_value = float(p_raw)
        except ValueError:
            continue

        if p_value < 0.01 and gene not in seen:
            genes.append(gene)
            seen.add(gene)

with open(out_path, "w") as out:
    for gene in genes:
        out.write(gene + "\n")

print(f"Selected genes for aBSREL with p_value < 0.01: {len(genes)}")
PY

GENE_COUNT=$(wc -l < "$GENE_LIST" | tr -d ' ')

echo "Running aBSREL on $GENE_COUNT genes using $ABSREL_WORKERS workers..."
echo

printf "gene\tstatus\tpath\telapsed_min\n" > "$RESULTS_TSV"

run_absrel() {
    local gene_name="$1"
    local start_ts end_ts elapsed_sec elapsed_min exit_code
    local codon_file tree_file codon_path tree_path output_json output_log

    start_ts=$(date +%s)

    codon_file="${gene_name}${FILE_SUFFIX}"
    tree_file="${gene_name}${TREE_SUFFIX}"

    codon_path="${PAL2NAL_DIR}/${codon_file}"
    tree_path="${LABELED_TREES_DIR}/${tree_file}"

    output_json="${ABSREL_DIR}/${gene_name}_absrel.json"
    output_log="${ABSREL_DIR}/${gene_name}_absrel.log"

    if [[ "$FORCE_RERUN" -eq 0 && -s "$output_json" ]]; then
        end_ts=$(date +%s)
        elapsed_sec=$((end_ts - start_ts))
        elapsed_min=$(awk -v s="$elapsed_sec" 'BEGIN { printf "%.2f", s / 60 }')
        printf "%s\t%s\t%s\t%s\n" "$gene_name" "skipped_existing" "$output_json" "$elapsed_min" >> "$RESULTS_TSV"
        echo "$gene_name: skipped_existing: $output_json"
        return 0
    fi

    if [[ ! -f "$codon_path" ]]; then
        end_ts=$(date +%s)
        elapsed_sec=$((end_ts - start_ts))
        elapsed_min=$(awk -v s="$elapsed_sec" 'BEGIN { printf "%.2f", s / 60 }')
        printf "%s\t%s\t%s\t%s\n" "$gene_name" "missing_alignment" "$codon_path" "$elapsed_min" >> "$RESULTS_TSV"
        echo "$gene_name: missing_alignment: $codon_path after ${elapsed_min} min"
        return 0
    fi

    if [[ ! -f "$tree_path" ]]; then
        end_ts=$(date +%s)
        elapsed_sec=$((end_ts - start_ts))
        elapsed_min=$(awk -v s="$elapsed_sec" 'BEGIN { printf "%.2f", s / 60 }')
        printf "%s\t%s\t%s\t%s\n" "$gene_name" "missing_tree" "$tree_path" "$elapsed_min" >> "$RESULTS_TSV"
        echo "$gene_name: missing_tree: $tree_path after ${elapsed_min} min"
        return 0
    fi

    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    "$HYPHY_BIN" CPU=1 absrel \
        --alignment "$codon_path" \
        --tree "$tree_path" \
        --output "$output_json" \
        > "$output_log" 2>&1

    exit_code=$?

    end_ts=$(date +%s)
    elapsed_sec=$((end_ts - start_ts))
    elapsed_min=$(awk -v s="$elapsed_sec" 'BEGIN { printf "%.2f", s / 60 }')

    if [[ "$exit_code" -ne 0 ]]; then
        printf "%s\t%s\t%s\t%s\n" "$gene_name" "failed" "$output_log" "$elapsed_min" >> "$RESULTS_TSV"
        echo "$gene_name: failed: $output_log after ${elapsed_min} min"
        return 0
    fi

    if [[ ! -s "$output_json" ]]; then
        printf "%s\t%s\t%s\t%s\n" "$gene_name" "no_json" "$output_log" "$elapsed_min" >> "$RESULTS_TSV"
        echo "$gene_name: no_json: $output_log after ${elapsed_min} min"
        return 0
    fi

    printf "%s\t%s\t%s\t%s\n" "$gene_name" "completed" "$output_json" "$elapsed_min" >> "$RESULTS_TSV"
    echo "Completed $gene_name in ${elapsed_min} min"
}

report_progress() {
    python3 - "$RESULTS_TSV" "$GENE_COUNT" "$ABSREL_WORKERS" <<'PY'
import sys
import statistics

path = sys.argv[1]
total = int(sys.argv[2])
workers = int(sys.argv[3])

rows = []

with open(path) as f:
    next(f, None)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 4:
            rows.append(parts)

done = len(rows)
times = []

for gene, status, path, elapsed in rows:
    if status == "completed":
        try:
            value = float(elapsed)
        except ValueError:
            continue
        if value > 0:
            times.append(value)

if times:
    median_time = statistics.median(times)
    remaining = max(total - done, 0)
    approx_remaining = (remaining / workers) * median_time
    print(f"Progress: {done}/{total}; median completed time: {median_time:.2f} min; approx remaining: {approx_remaining:.2f} min")
else:
    print(f"Progress: {done}/{total}")
PY
}

active_jobs=0

while IFS= read -r gene_name; do
    [[ -z "$gene_name" ]] && continue

    run_absrel "$gene_name" &

    active_jobs=$((active_jobs + 1))

    if (( active_jobs >= ABSREL_WORKERS )); then
        wait -n
        active_jobs=$((active_jobs - 1))
        report_progress
    fi
done < "$GENE_LIST"

while (( active_jobs > 0 )); do
    wait -n
    active_jobs=$((active_jobs - 1))
    report_progress
done

analyses_completed=$(awk -F'\t' 'NR > 1 && $2 == "completed" { count++ } END { print count + 0 }' "$RESULTS_TSV")
analyses_skipped=$(awk -F'\t' 'NR > 1 && $2 == "skipped_existing" { count++ } END { print count + 0 }' "$RESULTS_TSV")

echo
echo "Completed $analyses_completed new aBSREL analyses"
echo "Skipped $analyses_skipped existing aBSREL analyses"
echo
cat "$RESULTS_TSV"
echo
echo "Finished aBSREL run: $(date)"

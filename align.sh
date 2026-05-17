#!/usr/bin/env bash
set -u

NUC_DIR="pipeline_results/03_nucleotides/gene_specific_cds"
ALIGN_DIR="pipeline_results/04_alignmentss"
BUSTED_READY_DIR="pipeline_results/08_busted/ready_and_okay"

MACSE_BIN="/home/poroshina/.conda/envs/based/bin/macse"
MAX_WORKERS=12
TIMEOUT_HOURS=4

RUN_DIR="$ALIGN_DIR/_macse_run"
LOG_DIR="$RUN_DIR/logs"
mkdir -p "$ALIGN_DIR" "$RUN_DIR" "$LOG_DIR"

SUCCESS_FILE="$RUN_DIR/success.txt"
FAIL_FILE="$RUN_DIR/fail.txt"
SKIP_FILE="$RUN_DIR/skipped_has_busted.txt"
TOOLONG_GENES_FILE="$RUN_DIR/too_long_genes.txt"
MAIN_LOG="$RUN_DIR/macse_run.log"
SUMMARY_FILE="$RUN_DIR/summary.txt"

: > "$SUCCESS_FILE"
: > "$FAIL_FILE"
: > "$SKIP_FILE"
: > "$TOOLONG_GENES_FILE"
: > "$MAIN_LOG"
: > "$SUMMARY_FILE"

run_macse_one() {
    local nuc_path="$1"
    local nuc_file basename gene
    local output_nt output_aa
    local stdout_file stderr_file
    local rc

    nuc_file=$(basename "$nuc_path")
    basename="${nuc_file%.fasta}"
    gene="${basename%%_bird_only_nucleotide*}"

    output_nt="$ALIGN_DIR/${basename}_aligned_NT.fasta"
    output_aa="$ALIGN_DIR/${basename}_aligned_AA.fasta"

    stdout_file="$LOG_DIR/${basename}.stdout.log"
    stderr_file="$LOG_DIR/${basename}.stderr.log"

    if [[ -f "$output_nt" && -f "$output_aa" ]] && \
       find "$BUSTED_READY_DIR" -maxdepth 1 -type f -name "${gene}_*" | grep -q .; then
        printf '%s\n' "$gene" >> "$SKIP_FILE"
        printf '[SKIP] gene=%s reason=already_has_busted_result\n' "$gene" >> "$MAIN_LOG"
        return 0
    fi

    printf '[START] gene=%s file=%s\n' "$gene" "$nuc_file" >> "$MAIN_LOG"

    timeout "${TIMEOUT_HOURS}h" "$MACSE_BIN" \
        -prog alignSequences \
        -seq "$nuc_path" \
        -out_NT "$output_nt" \
        -out_AA "$output_aa" \
        >"$stdout_file" 2>"$stderr_file"

    rc=$?

    if [[ $rc -eq 0 && -s "$output_nt" && -s "$output_aa" ]]; then
        printf '%s\n' "$gene" >> "$SUCCESS_FILE"
        printf '[OK] gene=%s\n' "$gene" >> "$MAIN_LOG"

    elif [[ $rc -eq 124 ]]; then
        printf '%s\n' "$gene" >> "$TOOLONG_GENES_FILE"
        printf '[TOO_LONG] gene=%s timeout=%sh\n' "$gene" "$TIMEOUT_HOURS" >> "$MAIN_LOG"

    else
        printf '%s\n' "$gene" >> "$FAIL_FILE"
        {
            printf '[FAIL] gene=%s file=%s returncode=%s\n' "$gene" "$nuc_file" "$rc"
            printf 'stdout_log=%s\n' "$stdout_file"
            printf 'stderr_log=%s\n' "$stderr_file"
            printf '\n'
        } >> "$MAIN_LOG"
    fi
}

export ALIGN_DIR BUSTED_READY_DIR MACSE_BIN LOG_DIR
export SUCCESS_FILE FAIL_FILE SKIP_FILE TOOLONG_GENES_FILE MAIN_LOG
export TIMEOUT_HOURS
export -f run_macse_one

mapfile -d '' -t nuc_files < <(
    find "$NUC_DIR" -maxdepth 1 -type f -name "*.fasta" -print0 | sort -z
)

printf '%s\0' "${nuc_files[@]}" | \
    xargs -0 -n 1 -P "$MAX_WORKERS" bash -c 'run_macse_one "$1"' _

files_aligned=$(wc -l < "$SUCCESS_FILE")
files_failed=$(wc -l < "$FAIL_FILE")
files_toolong=$(wc -l < "$TOOLONG_GENES_FILE")
files_skipped=$(wc -l < "$SKIP_FILE")

{
    printf 'aligned=%s\n' "$files_aligned"
    printf 'failed=%s\n' "$files_failed"
    printf 'too_long=%s\n' "$files_toolong"
    printf 'skipped_has_busted=%s\n' "$files_skipped"
} > "$SUMMARY_FILE"

cat "$SUMMARY_FILE" >> "$MAIN_LOG"

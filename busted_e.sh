#!/usr/bin/env bash

set -u
set -o pipefail

# =========================
# НАСТРОЙКИ
# =========================
BUSTED_WORKERS=10
HYPHY_BIN="/home/poroshina/.conda/envs/based/bin/hyphy"

# 0 = не перезапускать json, которые уже есть в ready_and_okay
# 1 = перезапускать вообще всё, даже если json есть в ready_and_okay
FORCE_RERUN_BUSTED=0

# Максимальное время на один BUSTED job
# 4h = 4 часа
# Если хочешь отключить лимит, поставь BUSTED_TIMEOUT="0"
BUSTED_TIMEOUT="4h"

# HyPhy-аргумент: numerical errors будут warnings, а не fatal errors
HYPHY_ENV_ARG='ENV=TOLERATE_NUMERICAL_ERRORS=1;'

PAL2NAL_DIR="pipeline_results/05_pal2nal_normalized_no_duplicates"
LABELED_TREES_DIR="pipeline_results/07_labeled_trees"
BUSTED_DIR="pipeline_results/08_busted"

# ВАЖНО:
# готовыми считаются только json, которые лежат здесь
READY_AND_OKAY_DIR="${BUSTED_DIR}/ready_and_okay"

# =========================
# ПОДГОТОВКА
# =========================
mkdir -p "$BUSTED_DIR"
mkdir -p "$READY_AND_OKAY_DIR"

RUN_TS="$(date '+%Y%m%d_%H%M%S')"

MASTER_LOG="${BUSTED_DIR}/_busted_master_${RUN_TS}.log"
TASKS_FILE="${BUSTED_DIR}/_busted_tasks_${RUN_TS}.tsv"
STATUS_FILE="${BUSTED_DIR}/_busted_status_${RUN_TS}.tsv"
MISSING_TREES_FILE="${BUSTED_DIR}/_missing_trees_${RUN_TS}.txt"
EMPTY_FILES_FILE="${BUSTED_DIR}/_empty_codon_files_${RUN_TS}.txt"
SKIPPED_FILE="${BUSTED_DIR}/_skipped_existing_in_ready_and_okay_${RUN_TS}.txt"

: > "$MASTER_LOG"
: > "$TASKS_FILE"
: > "$STATUS_FILE"
: > "$MISSING_TREES_FILE"
: > "$EMPTY_FILES_FILE"
: > "$SKIPPED_FILE"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$MASTER_LOG"
}

find_tree_for_codon() {
    local codon_file="$1"
    local base tree_path

    base="${codon_file%_codon_aligned.fasta}"
    tree_path="${LABELED_TREES_DIR}/${base}_labeled.nwk"

    if [[ -f "$tree_path" ]]; then
        printf '%s\n' "$tree_path"
        return 0
    fi

    return 1
}

run_busted_job() {
    local codon_file="$1"
    local codon_path="$2"
    local tree_path="$3"
    local output_json="$4"
    local output_log="$5"

    local rc=0
    local status="FAIL"

    # ВАЖНО:
    # если файл есть только в BUSTED_DIR, но не в ready_and_okay,
    # он считается НЕготовым и будет перезаписан
    rm -f "$output_json"

    {
        echo "============================================================"
        echo "START: $(date '+%F %T')"
        echo "FILE:  $codon_file"
        echo "ALN:   $codon_path"
        echo "TREE:  $tree_path"
        echo "JSON:  $output_json"
        echo "LOG:   $output_log"
        echo "TIMEOUT: $BUSTED_TIMEOUT"
        echo "HYPHY_ENV_ARG: $HYPHY_ENV_ARG"
        echo "============================================================"

        if [[ "${BUSTED_TIMEOUT:-0}" == "0" ]]; then
            OMP_NUM_THREADS=1 \
            OPENBLAS_NUM_THREADS=1 \
            MKL_NUM_THREADS=1 \
            NUMEXPR_NUM_THREADS=1 \
            "$HYPHY_BIN" CPU=1 "$HYPHY_ENV_ARG" busted \
                --alignment "$codon_path" \
                --tree "$tree_path" \
                --branches Foreground \
                --error-sink Yes \
                --output "$output_json"
        else
            OMP_NUM_THREADS=1 \
            OPENBLAS_NUM_THREADS=1 \
            MKL_NUM_THREADS=1 \
            NUMEXPR_NUM_THREADS=1 \
            timeout --kill-after=60s "$BUSTED_TIMEOUT" \
            "$HYPHY_BIN" CPU=1 "$HYPHY_ENV_ARG" busted \
                --alignment "$codon_path" \
                --tree "$tree_path" \
                --branches Foreground \
                --error-sink Yes \
                --output "$output_json"
        fi

        rc=$?

        echo
        echo "END: $(date '+%F %T')"
        echo "RETURNCODE: $rc"

        if [[ $rc -eq 124 || $rc -eq 137 ]]; then
            echo "STATUS: TIMEOUT"
            echo "Reason: job exceeded BUSTED_TIMEOUT=$BUSTED_TIMEOUT"
        fi
    } > "$output_log" 2>&1

    if [[ $rc -eq 0 && -s "$output_json" ]]; then
        status="OK"
    elif [[ $rc -eq 124 || $rc -eq 137 ]]; then
        status="TIMEOUT"
        rm -f "$output_json"
    else
        status="FAIL"
    fi

    printf '%s\t%s\t%s\t%s\n' "$status" "$rc" "$codon_file" "$output_log" >> "$STATUS_FILE"
}

# =========================
# ПРОВЕРКИ
# =========================
if [[ ! -x "$HYPHY_BIN" ]]; then
    log "ERROR: hyphy не найден или не исполняемый: $HYPHY_BIN"
    exit 1
fi

if [[ ! -d "$PAL2NAL_DIR" ]]; then
    log "ERROR: не найдена папка с codon alignment: $PAL2NAL_DIR"
    exit 1
fi

if [[ ! -d "$LABELED_TREES_DIR" ]]; then
    log "ERROR: не найдена папка с labeled trees: $LABELED_TREES_DIR"
    exit 1
fi

if [[ "${BUSTED_TIMEOUT:-0}" != "0" ]] && ! command -v timeout >/dev/null 2>&1; then
    log "ERROR: не найдена команда timeout. На Linux она обычно есть в coreutils."
    exit 1
fi

log "Running BUSTED-E analysis..."
log "BUSTED_WORKERS=$BUSTED_WORKERS"
log "HYPHY_BIN=$HYPHY_BIN"
log "FORCE_RERUN_BUSTED=$FORCE_RERUN_BUSTED"
log "BUSTED_TIMEOUT=$BUSTED_TIMEOUT"
log "HYPHY_ENV_ARG=$HYPHY_ENV_ARG"
log "BUSTED_DIR=$BUSTED_DIR"
log "READY_AND_OKAY_DIR=$READY_AND_OKAY_DIR"
log "Only json files inside READY_AND_OKAY_DIR are treated as already completed."

# =========================
# СБОР ЗАДАЧ
# =========================
task_count=0
missing_count=0
empty_count=0
skipped_count=0

while IFS= read -r -d '' codon_path; do
    codon_file="$(basename "$codon_path")"

    [[ -f "$codon_path" ]] || continue

    if [[ ! -s "$codon_path" ]]; then
        echo "$codon_file" >> "$EMPTY_FILES_FILE"
        ((empty_count++))
        continue
    fi

    if ! tree_path="$(find_tree_for_codon "$codon_file")"; then
        echo "$codon_file" >> "$MISSING_TREES_FILE"
        ((missing_count++))
        continue
    fi

    output_json="${BUSTED_DIR}/${codon_file}_busted.json"
    output_log="${BUSTED_DIR}/${codon_file}_busted.log"

    ready_json="${READY_AND_OKAY_DIR}/${codon_file}_busted.json"

    # ВАЖНО:
    # скипаем только если json есть именно в ready_and_okay
    if [[ "$FORCE_RERUN_BUSTED" -eq 0 && -s "$ready_json" ]]; then
        echo "$codon_file" >> "$SKIPPED_FILE"
        ((skipped_count++))
        continue
    fi

    # если дошли сюда — значит файла нет в ready_and_okay
    # даже если файл есть в BUSTED_DIR, он будет перезаписан внутри run_busted_job

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$codon_file" \
        "$codon_path" \
        "$tree_path" \
        "$output_json" \
        "$output_log" >> "$TASKS_FILE"

    ((task_count++))

done < <(find "$PAL2NAL_DIR" -maxdepth 1 -type f -name '*.fasta' -print0 | sort -z)

log "Found $task_count BUSTED jobs to run"
log "Skipped existing in ready_and_okay: $skipped_count"
log "Missing trees: $missing_count"
log "Empty codon files: $empty_count"

# =========================
# ЗАПУСК
# =========================
running=0

if [[ -s "$TASKS_FILE" ]]; then
    while IFS=$'\t' read -r codon_file codon_path tree_path output_json output_log; do
        run_busted_job "$codon_file" "$codon_path" "$tree_path" "$output_json" "$output_log" &
        ((running++))

        if (( running >= BUSTED_WORKERS )); then
            wait -n
            ((running--))
        fi
    done < "$TASKS_FILE"

    while (( running > 0 )); do
        wait -n
        ((running--))
    done
fi

# =========================
# СВОДКА
# =========================
analyses_completed="$(awk -F'\t' '$1=="OK"{c++} END{print c+0}' "$STATUS_FILE")"
failed_count="$(awk -F'\t' '$1=="FAIL"{c++} END{print c+0}' "$STATUS_FILE")"
timeout_count="$(awk -F'\t' '$1=="TIMEOUT"{c++} END{print c+0}' "$STATUS_FILE")"

log "======================================================================"
log "BUSTED-E finished"
log "======================================================================"
log "Completed successfully: $analyses_completed"
log "Failed: $failed_count"
log "Timed out: $timeout_count"
log "Skipped existing in ready_and_okay: $skipped_count"
log "Missing trees: $missing_count"
log "Empty codon files: $empty_count"
log "Master log: $MASTER_LOG"
log "Status file: $STATUS_FILE"
log "Missing trees list: $MISSING_TREES_FILE"
log "Empty files list: $EMPTY_FILES_FILE"
log "Skipped existing in ready_and_okay list: $SKIPPED_FILE"

if (( failed_count > 0 || timeout_count > 0 )); then
    log "Failed/timeout jobs first 20:"
    shown=0

    while IFS=$'\t' read -r status rc codon_file output_log; do
        [[ "$status" == "FAIL" || "$status" == "TIMEOUT" ]] || continue

        printf '  %s | status=%s | returncode=%s\n    log: %s\n' \
            "$codon_file" "$status" "$rc" "$output_log" | tee -a "$MASTER_LOG"

        ((shown++))
        (( shown >= 20 )) && break
    done < "$STATUS_FILE"
fi
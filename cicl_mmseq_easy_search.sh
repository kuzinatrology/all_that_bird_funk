#!/bin/bash

set -uo pipefail

# Папка с входными fasta-файлами
QUERY_DIR="category_fastas_selection_criteria"

# База, по которой ищем
TARGET="all_birds_proteins_tagged_merged.faa"

# Общая папка для результатов
RESULTS_DIR="category_easy_search"
TMP_ROOT="${RESULTS_DIR}/tmp"
LOG_DIR="${RESULTS_DIR}/logs"

# Параметры easy-search
THREADS=8
PARAMS="-s 7 --min-seq-id 0.3 -e 1e-3 --alt-ali 30 --threads ${THREADS}"

mkdir -p "$RESULTS_DIR" "$TMP_ROOT" "$LOG_DIR"
shopt -s nullglob

# Цикл по всем fasta-файлам без исключений
for query in "$QUERY_DIR"/*.fasta; do
    [ -f "$query" ] || continue

    base=$(basename "$query")
    base_name="${base%.fasta}"

    output="${RESULTS_DIR}/homologs_${base_name}.m8"
    tmp_dir="${TMP_ROOT}/${base_name}"
    log="${LOG_DIR}/${base_name}.log"

    echo "Обработка: $query -> $output"
    mkdir -p "$tmp_dir"

    if mmseqs easy-search "$query" "$TARGET" "$output" "$tmp_dir" $PARAMS > "$log" 2>&1; then
        echo "  Успешно: $output"
    else
        echo "  Ошибка: см. $log"
    fi
done

echo "Все файлы обработаны. Результаты: $RESULTS_DIR"

#!/bin/bash

set -uo pipefail

# Папки с входными fasta-файлами (перебираем все)
QUERY_DIRS=(
    "human_category"
    "category_fastas_selection_criteria"
    "matrisom_and_tissue_regeneration_proteins"
)

# База, по которой ищем — primary-only, согласована с OrthoFinder
# (см. build_primary_merged.py). Раньше была all_birds_proteins_tagged_merged.faa
# (все изоформы) — это расходилось с primary-набором OrthoFinder.
TARGET="all_birds_proteins_primary_tagged_merged.faa"

# Общая папка для результатов
RESULTS_DIR="category_easy_search"
TMP_ROOT="${RESULTS_DIR}/tmp"
LOG_DIR="${RESULTS_DIR}/logs"

# Параметры easy-search
THREADS=8
PARAMS="-s 7 --min-seq-id 0.3 -e 1e-3 --alt-ali 30 --threads ${THREADS}"

mkdir -p "$RESULTS_DIR" "$TMP_ROOT" "$LOG_DIR"
shopt -s nullglob

# Цикл по всем папкам и всем fasta-файлам без исключений
for QUERY_DIR in "${QUERY_DIRS[@]}"; do
    if [ ! -d "$QUERY_DIR" ]; then
        echo "ПРОПУСК: папки нет — $QUERY_DIR"
        continue
    fi
    echo "=== Папка запросов: $QUERY_DIR ==="
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
done

echo "Все файлы обработаны. Результаты: $RESULTS_DIR"

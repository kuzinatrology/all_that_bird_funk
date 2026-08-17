#!/usr/bin/env python3
"""
Строит primary-only версию merged-протеома, чтобы mmseqs и OrthoFinder работали
на ОДНОМ наборе транскриптов.

Логика: собираем множество accession'ов из ortho_birds/.../primary_transcripts/
(по одной самой длинной изоформе на ген — их отобрал OrthoFinder-скрипт), затем
стримом фильтруем all_birds_proteins_tagged_merged.faa, оставляя только записи,
чей accession попал в это множество. NCBI-accession (NP_/XP_/...) уникален
глобально, поэтому фильтровать по нему безопасно без учёта видового тега.

Запуск из /home/poroshina:
    python dn_ds_pipeline_easy_search/build_primary_merged.py
"""
import glob

PRIMARY_GLOB = "ortho_birds/faa/primary_transcripts/*.faa"
MERGED = "dn_ds_pipeline_easy_search/all_birds_proteins_tagged_merged.faa"
OUT    = "dn_ds_pipeline_easy_search/all_birds_proteins_primary_tagged_merged.faa"


def acc_primary(header):        # ">NP_0012.1 desc [Species]" -> "NP_0012.1"
    return header[1:].split(None, 1)[0]


def acc_merged(header):         # ">Species|NP_0012.1 desc [...]" -> "NP_0012.1"
    after = header.split("|", 1)[1] if "|" in header else header[1:]
    return after.split(None, 1)[0]


allowed = set()
nfiles = 0
for f in sorted(glob.glob(PRIMARY_GLOB)):
    nfiles += 1
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"):
                allowed.add(acc_primary(line))
print(f"primary files: {nfiles}, уникальных accession: {len(allowed)}", flush=True)

total = kept = 0
with open(MERGED) as inp, open(OUT, "w") as out:
    write = False
    for line in inp:
        if line.startswith(">"):
            total += 1
            write = acc_merged(line) in allowed
            if write:
                kept += 1
                out.write(line)
        elif write:
            out.write(line)

print(f"merged записей: {total}", flush=True)
print(f"оставлено (primary): {kept}", flush=True)
print(f"primary-accession, не найденных в merged: {len(allowed) - kept}", flush=True)
print(f"OUT: {OUT}", flush=True)

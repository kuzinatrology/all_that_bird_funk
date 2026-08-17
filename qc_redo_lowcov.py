# QC-перегон: удаление последовательностей с покрытием <30% собственной CDS
# и пересчёт затронутых генов по всей цепочке раздела 4:
#   MACSE alignSequences -> exportAlignment -> HmmCleaner -> reportMaskAA2NT
#   -> фильтры видов/колонок (04e_filtered, порог колонок 0.5).
# Фоновый аналог ячейки 4.5 ноутбука + перегон. Запуск на lyon:
#   cd /home/poroshina/dn_ds_pipeline_easy_search
#   nohup /home/poroshina/.conda/envs/based/bin/python qc_redo_lowcov.py \
#       > qc_redo_lowcov.log 2>&1 < /dev/null &
# Всё инкрементально: повторный запуск доделывает недоделанное.
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from Bio import SeqIO

os.chdir("/home/poroshina/dn_ds_pipeline_easy_search")
OUTPUT_DIR = Path("pipeline_results")

MACSE_BIN  = "/home/poroshina/.conda/envs/based/bin/macse"
PERL       = "/usr/bin/perl"
HMMCLEANER = "/home/poroshina/perl5/bin/HmmCleaner.pl"
MACSE_ENV = dict(os.environ,
                 PATH="/home/poroshina/.conda/envs/based/bin:" + os.environ["PATH"])
HMM_ENV = dict(os.environ,
               PERL5LIB="/home/poroshina/perl5/lib/perl5",
               PATH=os.environ["PATH"] + ":/home/poroshina/.conda/envs/hmmcleaner/bin")

CDS_DIR    = OUTPUT_DIR / "03_nucleotides/gene_specific_cds"
PROT_DIR   = OUTPUT_DIR / "02_paired_proteins"
ALIGN_DIR  = OUTPUT_DIR / "04_alignmentss"
EXPORT_DIR = OUTPUT_DIR / "04b_exported"
STAT_DIR   = EXPORT_DIR / "stats"
CLEAN_DIR  = OUTPUT_DIR / "04c_hmmcleaner"
MASKED_DIR = OUTPUT_DIR / "04d_masked_NT"
DETAIL_DIR = MASKED_DIR / "mask_detail"
COLMAP_DIR = MASKED_DIR / "colmaps"
PAD_DIR    = MASKED_DIR / "aa_padded"
FILT_DIR   = OUTPUT_DIR / "04e_filtered"
FMAP_DIR   = FILT_DIR / "colmaps"
SUMMARY    = OUTPUT_DIR / "11_summary"
QC_CSV     = SUMMARY / "qc_lowcov_removed.csv"

SUF_NT  = "_aligned_NT.fasta"
SUF_AA  = "_aligned_AA.fasta"
SUF_HMM = "_aligned_AA_hmm.fasta"
GAPS = set("-.")

COV_THR      = 0.3   # вид доносит до 04d < 30% кодонов собственной CDS -> удалить из входов
MAX_SEQ_LOSS = 0.5   # фильтр видов (шаг 4.6)
MIN_COL_OCC  = 0.5   # фильтр колонок (шаг 4.6, решение 2026-08-15)
WORKERS      = 8

def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


# ---------- 1. детекция низкопокрытых ----------
def detect_lowcov():
    rows = []
    for ntp in sorted(MASKED_DIR.glob(f"*{SUF_NT}")):
        base = ntp.name[: -len(SUF_NT)]
        cds_p = CDS_DIR / f"{base}.fasta"
        if not cds_p.exists():
            continue
        orig = {r.id: len(r.seq) // 3 for r in SeqIO.parse(cds_p, "fasta")}
        for r in SeqIO.parse(ntp, "fasta"):
            if r.id not in orig or not orig[r.id]:
                continue
            s = str(r.seq)
            nongap = sum(1 for k in range(0, len(s), 3) if s[k:k+3] != "---")
            cov = nongap / orig[r.id]
            if cov < COV_THR:
                rows.append(dict(base=base, seq_id=r.id, coverage=round(cov, 3),
                                 orig_codons=orig[r.id], left_codons=nongap))
    return pd.DataFrame(rows, columns=["base", "seq_id", "coverage",
                                       "orig_codons", "left_codons"])


def remove_ids(fasta_path, drop_ids, by_acc=False):
    # в 02 имена видов с подвидом (Melospiza_melodia_melodia), в 03 без него —
    # общий ключ только accession после '|'
    if not fasta_path.exists():
        return 0
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    key = (lambda r: r.id.split("|")[-1]) if by_acc else (lambda r: r.id)
    keep = [r for r in recs if key(r) not in drop_ids]
    if len(keep) == len(recs):
        return 0
    with open(fasta_path, "w") as f:
        SeqIO.write(keep, f, "fasta")
    return len(recs) - len(keep)


def cleanup_gene(base):
    prot_base = base.replace("_bird_only_nucleotide", "_bird_only")
    victims = [
        ALIGN_DIR / f"{base}{SUF_NT}", ALIGN_DIR / f"{base}{SUF_AA}",
        EXPORT_DIR / f"{base}{SUF_NT}", EXPORT_DIR / f"{base}{SUF_AA}",
        STAT_DIR / f"{base}_stats.csv",
        CLEAN_DIR / f"{base}{SUF_AA}", CLEAN_DIR / f"{base}{SUF_HMM}",
        CLEAN_DIR / f"{base}{SUF_HMM.replace('.fasta', '.log')}",
        CLEAN_DIR / f"{base}{SUF_HMM.replace('.fasta', '.score')}",
        MASKED_DIR / f"{base}{SUF_NT}",
        COLMAP_DIR / f"{base}_colmap.csv",
        DETAIL_DIR / f"{base}_mask_detail.fasta",
        PAD_DIR / f"{base}{SUF_HMM}",
        FILT_DIR / f"{base}{SUF_NT}",
        FMAP_DIR / f"{base}_colmap.csv",
    ]
    n = 0
    for p in victims:
        if p.exists():
            p.unlink()
            n += 1
    return prot_base, n


log("=== 1. Детекция последовательностей с покрытием <", COV_THR, "===")
det = detect_lowcov()
log(f"найдено: {len(det)} последовательностей в {det.base.nunique()} генах")

if len(det):
    if QC_CSV.exists():
        old = pd.read_csv(QC_CSV)
        det = (pd.concat([old, det])
                 .drop_duplicates(subset=["base", "seq_id"])
                 .reset_index(drop=True))
    det.to_csv(QC_CSV, index=False)
    log(f"список зафиксирован: {QC_CSV} (всего записей: {len(det)})")

    log("=== 2. Удаление из входов (03 + 02) и зачистка производных ===")
    n03 = n02 = 0
    affected = sorted(det.base.unique())
    for base in affected:
        drop = set(det.loc[det.base == base, "seq_id"])
        n03 += remove_ids(CDS_DIR / f"{base}.fasta", drop)
        prot_base, _ = cleanup_gene(base)
        n02 += remove_ids(PROT_DIR / f"{prot_base}.fasta",
                          {i.split("|")[-1] for i in drop}, by_acc=True)
    log(f"удалено записей: из 03 — {n03}, из 02 — {n02}; "
        f"производные файлы {len(affected)} генов зачищены")
else:
    affected = []
    log("удалять нечего — идём доделывать недостающие производные")


# ---------- 3. MACSE alignSequences для генов без выравнивания ----------
def align_one(cds_p):
    base = cds_p.stem
    out_nt = ALIGN_DIR / f"{base}{SUF_NT}"
    out_aa = ALIGN_DIR / f"{base}{SUF_AA}"
    r = subprocess.run([MACSE_BIN, "-prog", "alignSequences", "-seq", str(cds_p),
                        "-out_NT", str(out_nt), "-out_AA", str(out_aa)],
                       capture_output=True, text=True, env=MACSE_ENV)
    ok = out_nt.exists() and out_nt.stat().st_size > 0
    return base, ok, (r.stderr or "")[-200:]

todo = [CDS_DIR / f"{b}.fasta" for b in affected
        if not (ALIGN_DIR / f"{b}{SUF_NT}").exists()]
log(f"=== 3. MACSE alignSequences: {len(todo)} генов ===")
t0 = time.perf_counter()
done_cnt = 0
def align_v(p):
    global done_cnt
    r = align_one(p)
    done_cnt += 1
    if done_cnt % 10 == 0:
        log(f"  MACSE: {done_cnt}/{len(todo)}")
    return r
with ThreadPoolExecutor(WORKERS) as ex:
    res = list(ex.map(align_v, todo))
bad = [x for x in res if not x[1]]
log(f"MACSE: {len(res) - len(bad)} ok, {len(bad)} упало, "
    f"{(time.perf_counter() - t0)/60:.1f} мин")
for b in bad[:10]:
    log("  FAIL:", b[0], b[2])


# ---------- 4. exportAlignment для генов без 04b ----------
def export_one(nt_path):
    base = nt_path.name[: -len(SUF_NT)]
    out_nt = EXPORT_DIR / f"{base}{SUF_NT}"
    if out_nt.exists() and out_nt.stat().st_size > 0:
        return base, "skip"
    r = subprocess.run([
        MACSE_BIN, "-prog", "exportAlignment", "-align", str(nt_path),
        "-codonForInternalStop", "---", "-codonForFinalStop", "---",
        "-codonForInternalFS", "---", "-codonForExternalFS", "---",
        "-charForRemainingFS", "-",
        "-out_NT", str(out_nt),
        "-out_AA", str(EXPORT_DIR / f"{base}{SUF_AA}"),
        "-out_stat_per_seq", str(STAT_DIR / f"{base}_stats.csv"),
    ], capture_output=True, text=True, env=MACSE_ENV)
    ok = out_nt.exists() and out_nt.stat().st_size > 0
    return base, "ok" if ok else f"FAIL rc={r.returncode}"

nt_files = sorted(ALIGN_DIR.glob(f"*{SUF_NT}"))
log(f"=== 4. exportAlignment: проверяю {len(nt_files)} генов ===")
with ThreadPoolExecutor(WORKERS) as ex:
    res = list(ex.map(export_one, nt_files))
cnt = pd.Series([s.split()[0] for _, s in res]).value_counts()
log("export:", dict(cnt))


# ---------- 5. HmmCleaner инкрементально ----------
import shutil
CLEAN_DIR.mkdir(parents=True, exist_ok=True)
all_aa = sorted(EXPORT_DIR.glob(f"*{SUF_AA}"))
todo = []
for p in all_aa:
    base = p.name[: -len(SUF_AA)]
    if not (CLEAN_DIR / f"{base}{SUF_HMM}").exists():
        shutil.copy(p, CLEAN_DIR / p.name)
        todo.append(CLEAN_DIR / p.name)
log(f"=== 5. HmmCleaner: {len(todo)} генов ===")

def clean_one(path):
    r = subprocess.run([PERL, HMMCLEANER, str(path), "--large", "--specificity"],
                       capture_output=True, text=True, env=HMM_ENV)
    base = path.name[: -len(SUF_AA)]
    out = CLEAN_DIR / f"{base}{SUF_HMM}"
    return path.name, r.returncode == 0 and out.exists() and out.stat().st_size > 0

with ThreadPoolExecutor(WORKERS) as ex:
    res = list(ex.map(clean_one, todo))
bad = [x for x in res if not x[1]]
log(f"HmmCleaner: {len(res) - len(bad)} ok, {len(bad)} упало")


# ---------- 6. reportMaskAA2NT инкрементально (с доращиванием хвоста) ----------
for d in (MASKED_DIR, DETAIL_DIR, COLMAP_DIR, PAD_DIR):
    d.mkdir(parents=True, exist_ok=True)

def mask_one(hmm_path):
    base   = hmm_path.name[: -len(SUF_HMM)]
    in_nt  = EXPORT_DIR / f"{base}{SUF_NT}"
    out_nt = MASKED_DIR / f"{base}{SUF_NT}"
    colmap = COLMAP_DIR / f"{base}_colmap.csv"
    if not in_nt.exists():
        return dict(base=base, status="no_NT")
    if out_nt.exists() and out_nt.stat().st_size > 0 and colmap.exists():
        return dict(base=base, status="skip")
    O = {r.id: str(r.seq) for r in SeqIO.parse(in_nt, "fasta")}
    H = {r.id: str(r.seq) for r in SeqIO.parse(hmm_path, "fasta")}
    ids = list(O)
    L = len(O[ids[0]]) // 3
    if set(H) != set(O):
        return dict(base=base, status="ID_MISMATCH")
    LH = len(H[ids[0]])
    if LH > L:
        return dict(base=base, status="AA_LONGER_THAN_NT")
    if LH < L:
        H = {i: s + "-" * (L - LH) for i, s in H.items()}
        aa_in = PAD_DIR / hmm_path.name
        with open(aa_in, "w") as f:
            for i in ids:
                f.write(f">{i}\n{H[i]}\n")
    else:
        aa_in = hmm_path
    r = subprocess.run([
        MACSE_BIN, "-prog", "reportMaskAA2NT",
        "-align", str(in_nt), "-align_AA", str(aa_in), "-mask_AA", "-",
        "-min_NT_to_keep_seq", "0", "-min_homology_to_keep_seq", "0",
        "-min_internal_homology_to_keep_seq", "0", "-min_seq_to_keep_site", "0",
        "-out_NT", str(out_nt),
        "-out_mask_detail", str(DETAIL_DIR / f"{base}_mask_detail.fasta"),
    ], capture_output=True, text=True, env=MACSE_ENV)
    if not (out_nt.exists() and out_nt.stat().st_size > 0):
        return dict(base=base, status=f"FAIL rc={r.returncode}")
    M = {r_.id: str(r_.seq) for r_ in SeqIO.parse(out_nt, "fasta")}
    empty = set(k for k in range(L) if all(H[i][k] in GAPS for i in ids))
    kept = [k for k in range(L) if k not in empty]
    dropped = set(O) - set(M)
    # сброс вида допустим, только если он замаскирован на 100%
    drop_ok = all(all(H[i][k] in GAPS for k in range(L)) for i in dropped)
    ok_width = all(len(M[i]) == 3 * len(kept) for i in M)
    pd.DataFrame({"new_codon_idx": range(len(kept)),
                  "old_codon_idx": kept}).to_csv(colmap, index=False)
    return dict(base=base,
                status="ok" if (drop_ok and ok_width) else "GEOM_MISMATCH")

hmm_files = sorted(CLEAN_DIR.glob(f"*{SUF_HMM}"))
log(f"=== 6. reportMaskAA2NT: проверяю {len(hmm_files)} генов ===")
with ThreadPoolExecutor(WORKERS) as ex:
    rows = list(ex.map(mask_one, hmm_files))
rm = pd.DataFrame(rows)
log("reportMaskAA2NT:", dict(rm.status.value_counts()))


# ---------- 7. фильтры видов и колонок (полный перегон, порог 0.5) ----------
FILT_DIR.mkdir(parents=True, exist_ok=True)
FMAP_DIR.mkdir(exist_ok=True)

def filter_one(nt_path):
    base = nt_path.name[: -len(SUF_NT)]
    out_nt = FILT_DIR / f"{base}{SUF_NT}"
    O = {r.id: str(r.seq) for r in SeqIO.parse(EXPORT_DIR / nt_path.name, "fasta")}
    D = {r.id: str(r.seq) for r in SeqIO.parse(nt_path, "fasta")}
    L0 = len(next(iter(O.values()))) // 3
    L1 = len(next(iter(D.values()))) // 3
    keep_ids, dropped = [], []
    for i, s in D.items():
        own0 = sum(1 for k in range(0, len(O[i]), 3) if O[i][k:k+3] != "---")
        own1 = sum(1 for k in range(0, len(s), 3) if s[k:k+3] != "---")
        if own0 and 1 - own1 / own0 > MAX_SEQ_LOSS:
            dropped.append(i)
        else:
            keep_ids.append(i)
    if not keep_ids:
        return dict(base=base, status="all_seqs_dropped")
    n_sp = len(keep_ids)
    keep_cols = [k for k in range(L1)
                 if sum(D[i][3*k:3*k+3] != "---" for i in keep_ids) >= MIN_COL_OCC * n_sp]
    with open(out_nt, "w") as f:
        for i in keep_ids:
            f.write(f">{i}\n" + "".join(D[i][3*k:3*k+3] for k in keep_cols) + "\n")
    cm = pd.read_csv(COLMAP_DIR / f"{base}_colmap.csv")
    to04b = dict(zip(cm.new_codon_idx, cm.old_codon_idx))
    pd.DataFrame({"final_codon_idx": range(len(keep_cols)),
                  "d04_codon_idx": keep_cols,
                  "b04_codon_idx": [to04b[k] for k in keep_cols],
                  }).to_csv(FMAP_DIR / f"{base}_colmap.csv", index=False)
    empty_left = [i for i in keep_ids
                  if all(D[i][3*k:3*k+3] == "---" for k in keep_cols)]
    return dict(base=base, status="ok" if not empty_left else "EMPTY_SEQ_LEFT",
                orig_codons=L0, final_codons=len(keep_cols),
                pct_of_orig=round(100 * len(keep_cols) / L0, 1),
                n_sp_final=n_sp, dropped_50pct=len(dropped))

nt_files = sorted(MASKED_DIR.glob(f"*{SUF_NT}"))
log(f"=== 7. Фильтры видов/колонок (порог {MIN_COL_OCC}): {len(nt_files)} генов ===")
with ThreadPoolExecutor(WORKERS) as ex:
    rows = list(ex.map(filter_one, nt_files))
ft = pd.DataFrame(rows)
log("фильтры:", dict(ft.status.value_counts()))
ok = ft[ft.status == "ok"]
log(f"колонок осталось от изначальной ширины 04b: медиана {ok.pct_of_orig.median():.1f}%, "
    f"среднее {ok.pct_of_orig.mean():.1f}%")
ft.to_csv(SUMMARY / "filter_seq_col_stats.csv", index=False)
log("DONE")

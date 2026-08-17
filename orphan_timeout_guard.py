# Сторож таймаута для осиротевших hyphy.
#
# После пересменки раннеров часть заданий BUSTED-E осталась без надзора: их
# родитель убит, значит шестичасовой лимит применить некому. Этот сторож
# доводит политику до конца — раз в 5 минут смотрит на сироты (hyphy c PPID=1)
# и убивает тех, кто перевалил за 6 часов ОТ СВОЕГО старта.
#
# Убитому гену: пустой json удаляется, в _busted_status.tsv дописывается строка
# TIMEOUT — чтобы ген не потерялся молча, а попал в список гигантов.
#
# Живые дети текущего раннера (PPID != 1) не трогаются: за ними следит он сам.
# Выходит, когда сирот не осталось.
#
# Запуск на lyon:
#   cd /home/poroshina/dn_ds_pipeline_easy_search
#   nohup /home/poroshina/.conda/envs/based/bin/python orphan_timeout_guard.py \
#       > orphan_timeout_guard.log 2>&1 < /dev/null &
import os
import signal
import subprocess
import time
from pathlib import Path

OUT_DIR = Path("/home/poroshina/dn_ds_pipeline_easy_search/pipeline_results/08_busted")
STATUS_TSV = OUT_DIR / "_busted_status.tsv"
SUF_ALN = "_binomial_NT.fasta"
TIMEOUT_SEC = 6 * 3600
POLL_SEC = 300


def log(*a):
    print(time.strftime("[%F %T] ") + " ".join(str(x) for x in a), flush=True)


def orphans():
    """hyphy-задания BUSTED без родителя (PPID=1) -> (pid, секунд_живёт, ген)."""
    out = subprocess.run(["ps", "-eo", "pid,ppid,etimes,args", "--no-headers"],
                         capture_output=True, text=True).stdout
    res = []
    for line in out.splitlines():
        parts = line.split(None, 3)
        if len(parts) < 4:
            continue
        pid, ppid, etimes, args = int(parts[0]), int(parts[1]), int(parts[2]), parts[3]
        if "/bin/hyphy" not in args or "busted" not in args or ppid != 1:
            continue
        gene = None
        toks = args.split()
        if "--alignment" in toks:
            name = Path(toks[toks.index("--alignment") + 1]).name
            if name.endswith(SUF_ALN):
                gene = name[: -len(SUF_ALN)]
        res.append((pid, etimes, gene))
    return res


log("сторож запущен: лимит 6 ч, опрос раз в", POLL_SEC, "с")
while True:
    orp = orphans()
    if not orp:
        log("сирот не осталось — сторож больше не нужен, выхожу")
        break

    overdue = [o for o in orp if o[1] > TIMEOUT_SEC]
    if overdue:
        for pid, etimes, gene in overdue:
            log(f"TIMEOUT: {gene} (pid {pid}, {etimes/3600:.2f} ч) — убиваю")
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            if gene:
                j = OUT_DIR / f"{gene}_busted.json"
                if j.exists() and j.stat().st_size == 0:
                    j.unlink()
                with open(STATUS_TSV, "a") as f:
                    f.write(f"{gene}\tTIMEOUT\t-\t{etimes}\t\n")
    else:
        oldest = max(o[1] for o in orp)
        log(f"сирот: {len(orp)}, самой старой {oldest/3600:.2f} ч, "
            f"до лимита {(TIMEOUT_SEC - oldest)/3600:.2f} ч")

    time.sleep(POLL_SEC)

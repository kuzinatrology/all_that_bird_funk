# BUSTED-E на всех генах из 05d_binomial.
#
#   hyphy CPU=1 ENV=TOLERATE_NUMERICAL_ERRORS=1; busted
#       --alignment 05d_binomial/<ген>_binomial_NT.fasta
#       --tree      07_labeled_trees/<ген>_labeled.nwk
#       --branches Foreground --error-sink Yes
#
# Порядок: сначала 24 гена из priority_genes.txt (в порядке файла), потом
# остальные по алфавиту.
#
# Ядра: жёсткий потолок MAX_WORKERS, ровно 1 процесс на ядро (CPU=1 +
# OMP/OPENBLAS/MKL/NUMEXPR=1). Раз в POLL_SEC секунд измеряется реальная
# загрузка машины по /proc/stat; всё, что занято сверх наших задач, считается
# чужим, и наш потолок опускается на эту величину. Уже запущенные задачи не
# убиваются (это часы счёта), но все идут под nice, поэтому планировщик и так
# отдаёт ядро чужому процессу. Когда чужие уходят — потолок сам поднимается.
#
# Таймаут TIMEOUT_SEC на ген; недосчитанный json удаляется, ген помечается
# TIMEOUT. Тайминги пишутся по каждому гену — по ним потом отсечём гигантов.
#
# Запуск на lyon:
#   cd /home/poroshina/dn_ds_pipeline_easy_search
#   nohup /home/poroshina/.conda/envs/based/bin/python busted_e_run.py \
#       > busted_e_run.log 2>&1 < /dev/null &
#
# Перезапуск безопасен: ген с непустым json в 08_busted пропускается, то есть
# продолжаем с места обрыва. На старте папка пустая (старая уехала в backup),
# поэтому в этом прогоне пересчитываются ВСЕ гены без исключения.
import os
import shutil
import signal
import subprocess
import time
from pathlib import Path

BASE = Path("/home/poroshina/dn_ds_pipeline_easy_search")
os.chdir(BASE)

HYPHY = "/home/poroshina/.conda/envs/based/bin/hyphy"
ALIGN_DIR = BASE / "pipeline_results/05d_binomial"
TREE_DIR  = BASE / "pipeline_results/07_labeled_trees"
OUT_DIR   = BASE / "pipeline_results/08_busted"
PRIORITY  = BASE / "priority_genes.txt"

SUF_ALN  = "_binomial_NT.fasta"
SUF_TREE = "_labeled.nwk"

MAX_WORKERS = 18        # стартовый потолок; меняется на лету через файл ниже
MIN_WORKERS = 2         # ниже не опускаемся, иначе прогон не кончится
# Живое управление: записать число в этот файл — и потолок сменится на
# следующем опросе, без перезапуска.  echo 24 > .../08_busted/_max_workers
CAP_FILE_NAME = "_max_workers"
TOTAL_CORES = os.cpu_count() or 32
RESERVE     = 2         # ядра, которые всегда оставляем машине
NICE        = 10        # уступаем CPU чужим процессам
TIMEOUT_SEC = 6 * 3600
POLL_SEC    = 60        # как часто пересматриваем потолок
TICK        = 5

STATUS_TSV = OUT_DIR / "_busted_status.tsv"
RUN_LOG    = OUT_DIR / "_busted_run.log"


def log(*a):
    line = time.strftime("[%F %T] ") + " ".join(str(x) for x in a)
    print(line, flush=True)
    with open(RUN_LOG, "a") as f:
        f.write(line + "\n")


OUR_UID = os.getuid()
CLK_TCK = os.sysconf("SC_CLK_TCK")


def cpu_snapshot():
    """Счётчик процессорного времени по КАЖДОМУ процессу + владелец.

    Считаем именно так, а не через /proc/stat: нам нужна нагрузка ТОЛЬКО чужих
    пользователей. Свои процессы (и BUSTED-E, и любая другая работа AP) в расчёт
    уступания не идут — своей загрузкой она распоряжается сама.
    ps pcpu тут не годится: это среднее за всю жизнь процесса, а не текущее.
    """
    out = {}
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            uid = os.stat(f"/proc/{pid}").st_uid
            if uid == OUR_UID:
                continue                      # свои не считаем вовсе
            with open(f"/proc/{pid}/stat") as f:
                fields = f.read().rsplit(") ", 1)[-1].split()
            out[pid] = int(fields[11]) + int(fields[12])   # utime + stime, тики
        except (FileNotFoundError, ProcessLookupError, PermissionError, IndexError):
            continue
    return time.time(), out


def others_cores(prev, cur):
    """Сколько ядер съели ЧУЖИЕ пользователи за интервал между снимками."""
    dt = cur[0] - prev[0]
    if dt <= 0:
        return 0.0
    ticks = sum(v - prev[1].get(pid, v) for pid, v in cur[1].items())
    return max(0.0, ticks / CLK_TCK / dt)


def our_busted_processes():
    """Сколько заданий BUSTED-E крутится ОТ НАС всего — включая осиротевшие от
    прошлых раннеров. Иначе после пересменки сироты и новые дети складываются,
    и лимит в 18 незаметно превращается в 36."""
    out = subprocess.run(["ps", "-u", str(OUR_UID), "-o", "args", "--no-headers"],
                         capture_output=True, text=True).stdout
    return sum(1 for l in out.splitlines() if "/bin/hyphy" in l and "busted" in l)


def read_cap_file(current):
    """Живой потолок из файла. Кривое содержимое игнорируем, не падаем."""
    p = OUT_DIR / CAP_FILE_NAME
    if not p.exists():
        return current
    try:
        v = int(p.read_text().strip())
    except (ValueError, OSError):
        return current
    return max(MIN_WORKERS, min(v, 64))


def count_our_hyphy():
    """Сколько hyphy крутится ОТ НАС всего — включая осиротевшие от прошлых
    раннеров. Потолок считается по этому числу, иначе после перезапуска
    сироты и новые дети складываются и лимит по ядрам улетает вдвое."""
    out = subprocess.run(["ps", "-eo", "args", "--no-headers"],
                         capture_output=True, text=True).stdout
    return sum(1 for l in out.splitlines() if "/bin/hyphy" in l and "busted" in l)


def inflight_genes():
    """Гены, которые уже считает кто-то другой (осиротевшие hyphy от прошлого
    раннера). Их не берём в работу, иначе будем считать один ген дважды."""
    out = subprocess.run(["ps", "-eo", "args", "--no-headers"],
                         capture_output=True, text=True).stdout
    genes = set()
    for line in out.splitlines():
        if "/bin/hyphy" not in line or "--alignment" not in line:
            continue
        parts = line.split()
        try:
            aln = parts[parts.index("--alignment") + 1]
        except (ValueError, IndexError):
            continue
        name = Path(aln).name
        if name.endswith(SUF_ALN):
            genes.add(name[: -len(SUF_ALN)])
    return genes


def build_tasks():
    prio = []
    if PRIORITY.exists():
        for line in PRIORITY.read_text().splitlines():
            g = line.strip()
            if g:
                prio.append(f"{g}_bird_only_nucleotide")
    all_genes = sorted(p.name[: -len(SUF_ALN)] for p in ALIGN_DIR.glob(f"*{SUF_ALN}"))
    seen, ordered = set(), []
    for g in prio:                       # приоритетные первыми, в порядке файла
        if g in set(all_genes) and g not in seen:
            seen.add(g)
            ordered.append(g)
    n_prio = len(ordered)
    for g in all_genes:
        if g not in seen:
            seen.add(g)
            ordered.append(g)
    return ordered, n_prio


def start_job(gene):
    aln  = ALIGN_DIR / f"{gene}{SUF_ALN}"
    tree = TREE_DIR / f"{gene}{SUF_TREE}"
    out_json = OUT_DIR / f"{gene}_busted.json"
    out_log  = OUT_DIR / f"{gene}_busted.log"
    if not tree.exists():
        return None, "NO_TREE", out_json
    env = dict(os.environ,
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    cmd = [HYPHY, "CPU=1", "ENV=TOLERATE_NUMERICAL_ERRORS=1;", "busted",
           "--alignment", str(aln),
           "--tree", str(tree),
           "--branches", "Foreground",
           "--error-sink", "Yes",
           "--output", str(out_json)]
    fh = open(out_log, "w")
    fh.write(" ".join(cmd) + "\n\n")
    fh.flush()
    p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env,
                         start_new_session=True,          # свой pgid -> убиваем группой
                         preexec_fn=lambda: os.nice(NICE))
    return (p, fh, time.time(), out_json), "STARTED", out_json


def kill_group(p):
    try:
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)
    except ProcessLookupError:
        return
    for _ in range(60):
        if p.poll() is not None:
            return
        time.sleep(1)
    try:
        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
    except ProcessLookupError:
        pass


OUT_DIR.mkdir(parents=True, exist_ok=True)
if not STATUS_TSV.exists():
    STATUS_TSV.write_text("gene\tstatus\trc\tseconds\tjson\n")

tasks, n_prio = build_tasks()
busy = inflight_genes()
todo = []
done_before = 0
for g in tasks:
    j = OUT_DIR / f"{g}_busted.json"
    if j.exists() and j.stat().st_size > 0:
        done_before += 1
    elif g in busy:
        continue                      # считается прошлым раннером, не дублируем
    else:
        todo.append(g)

log("=" * 64)
log("BUSTED-E: старт")
log(f"выравнивания: {ALIGN_DIR}")
log(f"деревья:      {TREE_DIR}")
log(f"результаты:   {OUT_DIR}")
log(f"генов всего: {len(tasks)} | приоритетных первыми: {n_prio} | "
    f"уже готово (пропуск): {done_before} | считает прошлый раннер: {len(busy)} | "
    f"в работу: {len(todo)}")
if busy:
    log("  доверено прошлому раннеру: " + ", ".join(sorted(busy)))
log(f"ядра: потолок {MAX_WORKERS}, минимум {MIN_WORKERS}, всего в машине "
    f"{TOTAL_CORES}, резерв {RESERVE}, nice {NICE}")
log(f"потолок меняется на лету: echo N > {OUT_DIR / CAP_FILE_NAME}")
log(f"таймаут на ген: {TIMEOUT_SEC/3600:.1f} ч")
log("=" * 64)

running = {}          # gene -> (proc, filehandle, t0, json)
queue = list(todo)
hard_cap = read_cap_file(MAX_WORKERS)

# Потолок считаем ДО первого запуска, иначе раннер успевает выстрелить полным
# залпом поверх чужой нагрузки и наших же сирот, а опрос случится только через
# минуту (так и вышло 2026-08-16: 18 новых поверх 30 сирот = 48 на 32 ядрах).
prev_stat = cpu_snapshot()
time.sleep(3)
_cur = cpu_snapshot()
_others = others_cores(prev_stat, _cur)
prev_stat = _cur
_ours = our_busted_processes()
cap = max(0, min(hard_cap, max(MIN_WORKERS,
                               int(TOTAL_CORES - RESERVE - _others))) - _ours)
last_poll = time.time()
log(f"стартовый потолок: {cap} (чужая нагрузка {_others:.1f} ядер, "
    f"наших заданий уже крутится {_ours})")
n_ok = n_fail = n_timeout = n_notree = 0
t_start = time.time()

while queue or running:
    now = time.time()

    # пересмотр потолка по реальной загрузке машины
    if now - last_poll >= POLL_SEC:
        cur = cpu_snapshot()
        others = others_cores(prev_stat, cur)
        prev_stat, last_poll = cur, now

        new_max = read_cap_file(hard_cap)
        if new_max != hard_cap:
            log(f"потолок из файла {CAP_FILE_NAME}: {hard_cap} -> {new_max}")
            hard_cap = new_max

        # сколько заданий BUSTED-E допустимо держать ВСЕГО (с учётом сирот)
        cap_total = min(hard_cap, max(MIN_WORKERS,
                                      int(TOTAL_CORES - RESERVE - others)))
        ours_elsewhere = max(0, our_busted_processes() - len(running))
        new_cap = max(0, cap_total - ours_elsewhere)

        if new_cap != cap:
            extra = f", сирот наших {ours_elsewhere}" if ours_elsewhere else ""
            log(f"потолок воркеров: {cap} -> {new_cap} "
                f"(чужая нагрузка {others:.1f} ядер{extra})")
            cap = new_cap

    # запуск новых
    while queue and len(running) < cap:
        gene = queue.pop(0)
        job, state, out_json = start_job(gene)
        if state == "NO_TREE":
            n_notree += 1
            with open(STATUS_TSV, "a") as f:
                f.write(f"{gene}\tNO_TREE\t-\t0\t\n")
            continue
        running[gene] = job

    # проверка завершившихся и таймаутов
    for gene in list(running):
        p, fh, t0, out_json = running[gene]
        rc = p.poll()
        elapsed = time.time() - t0
        if rc is None:
            if elapsed > TIMEOUT_SEC:
                log(f"TIMEOUT {gene} ({elapsed/3600:.2f} ч) — убиваю")
                kill_group(p)
                fh.close()
                if out_json.exists():
                    out_json.unlink()
                n_timeout += 1
                with open(STATUS_TSV, "a") as f:
                    f.write(f"{gene}\tTIMEOUT\t-\t{elapsed:.0f}\t\n")
                del running[gene]
            continue
        fh.close()
        ok = rc == 0 and out_json.exists() and out_json.stat().st_size > 0
        status = "OK" if ok else "FAIL"
        if ok:
            n_ok += 1
        else:
            n_fail += 1
            log(f"FAIL {gene} rc={rc} (см. {gene}_busted.log)")
        with open(STATUS_TSV, "a") as f:
            f.write(f"{gene}\t{status}\t{rc}\t{elapsed:.0f}\t"
                    f"{out_json.name if ok else ''}\n")
        del running[gene]
        n_done = n_ok + n_fail + n_timeout
        if n_done % 25 == 0:
            el = time.time() - t_start
            rate = n_done / el if el else 0
            left = (len(queue) + len(running)) / rate / 3600 if rate else 0
            log(f"прогресс: {n_done}/{len(todo)} | ok={n_ok} fail={n_fail} "
                f"timeout={n_timeout} | в работе {len(running)}/{cap} | "
                f"осталось ~{left:.1f} ч")

    time.sleep(TICK)

log("=" * 64)
log(f"BUSTED-E завершён за {(time.time()-t_start)/3600:.2f} ч")
log(f"ok={n_ok} fail={n_fail} timeout={n_timeout} no_tree={n_notree}")
log(f"статусы и тайминги: {STATUS_TSV}")
log("DONE")

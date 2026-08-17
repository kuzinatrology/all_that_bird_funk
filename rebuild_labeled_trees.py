# Пересборка размеченных деревьев для ВСЕХ генов из 05d_binomial.
# Серверный аналог ячейки 5 ноутбука: обрезка видового дерева pruned_tree.nwk
# под виды каждого гена + разметка форграунда, БЕЗ перезаписи выравниваний.
#
# Понадобилось, потому что 2026-08-16 push.sh затёр 943 свежих дерева майскими
# копиями с ноутбука (см. METHODS_CHANGELOG). Здесь же — проверка, что для
# каждого гена множество листьев дерева ровно равно множеству видов в
# выравнивании: без этого HyPhy падает с "number of tree tips is not equal to
# the number of sequences".
import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from ete3 import Tree

BASE = Path("/home/poroshina/dn_ds_pipeline_easy_search")
os.chdir(BASE)

ALIGN_DIR = BASE / "pipeline_results/05d_binomial"
TREE_DIR  = BASE / "pipeline_results/07_labeled_trees"
SPECIES_TREE = BASE / "pruned_tree.nwk"
SUF_ALN  = "_binomial_NT.fasta"
SUF_TREE = "_labeled.nwk"

# форграунд берём из уже посчитанной ячейкой таблицы, чтобы не тащить сюда df
FG_SRC = BASE / "pipeline_results/11_summary/tree_labeling_stats.csv"

full_tree = Tree(str(SPECIES_TREE), format=1)
TREE_LEAVES = {l.name for l in full_tree.get_leaves()}

# восстанавливаем множество форграунд-видов из уцелевших свежих деревьев:
# берём любое дерево, записанное ячейкой, и читаем метки {Foreground}
fg_species = set()
for t in sorted(TREE_DIR.glob(f"*{SUF_TREE}")):
    txt = t.read_text()
    if "{Foreground}" not in txt:
        continue
    tr = Tree(str(t), format=1)
    names = [l.name for l in tr.get_leaves()]
    if all("{" in n for n in names):          # дерево от новой ячейки
        fg_species = {n.split("{")[0] for n in names if n.endswith("{Foreground}")}
        # добираем со всех свежих деревьев, чтобы не зависеть от состава одного гена
        for t2 in sorted(TREE_DIR.glob(f"*{SUF_TREE}"))[:400]:
            tr2 = Tree(str(t2), format=1)
            ns = [l.name for l in tr2.get_leaves()]
            if all("{" in n for n in ns):
                fg_species |= {n.split("{")[0] for n in ns if n.endswith("{Foreground}")}
        break
if not fg_species:
    raise SystemExit("не удалось восстановить список форграунд-видов из деревьев")
print(f"форграунд-видов восстановлено: {len(fg_species)}")
print("  ", ", ".join(sorted(fg_species)))

rows = []
for p in sorted(ALIGN_DIR.glob(f"*{SUF_ALN}")):
    gene = p.name[: -len(SUF_ALN)]
    species = [r.id for r in SeqIO.parse(p, "fasta")]
    keep = sorted(set(species) & TREE_LEAVES)
    if len(keep) != len(species):
        rows.append(dict(gene=gene, status="SPECIES_NOT_IN_TREE",
                         n_aln=len(species), n_leaves=0, n_fg=0))
        continue
    t = full_tree.copy()
    t.prune(keep, preserve_branch_length=True)
    n_fg = 0
    for leaf in t.get_leaves():
        if leaf.name in fg_species:
            leaf.name = f"{leaf.name}{{Foreground}}"
            n_fg += 1
        else:
            leaf.name = f"{leaf.name}{{Background}}"
    out = TREE_DIR / f"{gene}{SUF_TREE}"
    t.write(outfile=str(out), format=1)

    check = {l.name.split("{")[0] for l in Tree(str(out), format=1).get_leaves()}
    rows.append(dict(gene=gene,
                     status="ok" if check == set(species) else "MISMATCH",
                     n_aln=len(species), n_leaves=len(check), n_fg=n_fg))

df = pd.DataFrame(rows)
print()
print("генов:", len(df))
print(df.status.value_counts().to_string())
bad = df[df.status != "ok"]
if len(bad):
    print("ПРОБЛЕМНЫЕ:")
    print(bad.head(20).to_string(index=False))
ok = df[df.status == "ok"]
print(f"листьев на дерево: медиана {ok.n_leaves.median():.0f}")
print(f"форграунд-ветвей: медиана {ok.n_fg.median():.0f}, "
      f"минимум {ok.n_fg.min()}, максимум {ok.n_fg.max()}")
df.to_csv(BASE / "pipeline_results/11_summary/tree_rebuild_check.csv", index=False)
print("DONE")

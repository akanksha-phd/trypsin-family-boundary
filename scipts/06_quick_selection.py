#Step 6: quick-and-dirty selection test.

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO

FAMILY = "trypsin-like"   # or "active-S1"
MAX_DIST = 0.40           # a member is "covered" if within 60% identity of a selected one
NS = [3, 5, 8, 10, 15, 20, 30, 50]
N_RANDOM = 100
MIN_OVERLAP = 100         # minimum aligned columns to trust an identity value

t = pd.read_csv("results/triad_check.tsv", sep="\t")
ok = t.triad_intact if FAMILY == "active-S1" else (t.triad_intact & (t.S1_res == "D"))
fam = t[ok].copy()
print(f"family definition '{FAMILY}': {len(fam)} of {len(t)} hits kept (+ seed)")

aln = {r.id: np.array(list(str(r.seq))) for r in AlignIO.read("results/domains.aln.fasta", "fasta")}
ids = ["PRSS1_domain"] + fam.target.tolist()
X = np.vstack([aln[i] for i in ids])
n = len(ids)

D = np.ones((n, n))
for i in range(n):
    both = (X[i] != "-") & (X != "-")
    same = (X[i] == X) & both
    ov = both.sum(axis=1)
    D[i] = np.where(ov >= MIN_OVERLAP, 1 - same.sum(axis=1) / np.maximum(ov, 1), 1.0)
np.fill_diagonal(D, 0.0)
D = np.minimum(D, D.T)

meta = fam.set_index("target")
cls = np.array(["trypsin"] + meta.loc[ids[1:], "class"].tolist())
genus = np.array(["Homo"] + [str(o).split()[0] for o in meta.loc[ids[1:], "organism"]])


def metrics(sel):
    nearest = D[:, sel].min(axis=1)
    return (nearest <= MAX_DIST).mean(), nearest.mean(), len(set(cls[sel]))


def farthest_point(k, start):
    sel = [start]
    nearest = D[:, start].copy()
    while len(sel) < k:
        j = int(np.argmax(nearest))
        sel.append(j)
        nearest = np.minimum(nearest, D[:, j])
    return sel


def stratified(k, rng):
    groups = {}
    for i, g in enumerate(genus):
        groups.setdefault(g, []).append(i)
    order = list(groups)
    rng.shuffle(order)
    pools = {g: list(rng.permutation(v)) for g, v in groups.items()}
    sel = []
    while len(sel) < k and any(pools.values()):
        for g in order:
            if pools[g] and len(sel) < k:
                sel.append(int(pools[g].pop()))
    return sel


medoid = int(np.argmin(D.mean(axis=1)))
rng = np.random.default_rng(0)
rows = []
for k in NS:
    if k >= n:
        continue
    r = np.array([metrics(list(rng.choice(n, k, replace=False))) for _ in range(N_RANDOM)])
    rows.append(("random", k, *r.mean(axis=0)))
    s = np.array([metrics(stratified(k, np.random.default_rng(s_))) for s_ in range(N_RANDOM)])
    rows.append(("taxonomy-stratified", k, *s.mean(axis=0)))
    rows.append(("max-diversity", k, *metrics(farthest_point(k, medoid))))

res = pd.DataFrame(rows, columns=["strategy", "N", "coverage", "mean_dist", "classes"])
res.to_csv("results/quick_selection.tsv", sep="\t", index=False)
print(res.round(3).to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for s, g in res.groupby("strategy"):
    ax[0].plot(g.N, g.coverage, marker="o", label=s)
    ax[1].plot(g.N, g.classes, marker="o", label=s)
ax[0].set_xlabel("Sequences selected (N)")
ax[0].set_ylabel(f"Fraction of family within {int((1-MAX_DIST)*100)}% identity")
ax[1].set_xlabel("Sequences selected (N)")
ax[1].set_ylabel("Functional classes represented")
ax[0].legend()
plt.suptitle(f"Quick-selection strategies ({FAMILY} family, n={n})")
plt.tight_layout()
plt.savefig("figures/fig3_strategy_comparison.png", dpi=200)
print("saved figures/fig3_strategy_comparison.png and results/quick_selection.tsv")

"""Step 8: elbow / saturation analysis (pangenome-style).

(a) Feature retention vs identity: for all filtered hits, the fraction that still carry an
    intact catalytic triad and Asp at the S1 position, by identity bin. An elbow here is a
    biological argument for a boundary (where the defining features of the domain disappear).
(b) Rarefaction: sample N sequences from the trypsin-like family and count distinct clusters
    at 50% and 70% identity. A curve that flattens means more sequences add little new
    diversity. The Heaps'-law exponent (clusters ~ N^gamma) summarises openness:
    gamma near 0 = saturating ("closed"), gamma near 1 = every new sequence is new diversity.
Elbow = point of maximum distance from the straight line joining the curve's endpoints
(a simple "kneedle"). Treat elbows as descriptive, not as proof.
Run from the repo root:  python scripts/08_elbow_analysis.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO

MIN_OVERLAP = 100
N_DRAW = 20
THRESHOLDS = [0.7, 0.5]

# concave curve: furthest above the diagonal
def elbow(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    xn = (x - x.min()) / (x.max() - x.min())
    yn = (y - y.min()) / (y.max() - y.min())
    return int(np.argmax(yn - xn))   


# (a) feature retention vs identity (all filtered hits)
t = pd.read_csv("results/triad_check.tsv", sep="\t")
t["features_ok"] = t.triad_intact & (t.S1_res == "D")
bins = [0, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 101]
t["bin"] = pd.cut(t.pident, bins)
fa = t.groupby("bin", observed=True).agg(n=("target", "size"), frac_features=("features_ok", "mean"),
                                       frac_triad=("triad_intact", "mean")).reset_index()
fa.to_csv("results/feature_retention.tsv", sep="\t", index=False)
print("Feature retention by identity bin:")
print(fa.round(3).to_string(index=False))

# (b) rarefaction on the trypsin-like family
ok = t.triad_intact & (t.S1_res == "D")
fam = t[ok]
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


def n_clusters(sub, d):
    near = sub <= d
    left = np.ones(len(sub), bool)
    k = 0
    while left.any():
        counts = (near & left[None, :]).sum(axis=1)
        counts[~left] = -1
        r = int(np.argmax(counts))
        left &= ~near[r]
        k += 1
    return k


grid = sorted(set([5, 10, 20, 40, 80, 120, 160, 240, 320, 440, n]))
grid = [g for g in grid if g <= n]
rng = np.random.default_rng(0)
rows = []
for g in grid:
    for thr in THRESHOLDS:
        vals = []
        for _ in range(N_DRAW if g < n else 1):
            idx = rng.choice(n, g, replace=False)
            vals.append(n_clusters(D[np.ix_(idx, idx)], 1 - thr))
        rows.append((g, int(thr * 100), np.mean(vals), np.std(vals)))
rare = pd.DataFrame(rows, columns=["N", "identity_pct", "mean_clusters", "sd"])
rare.to_csv("results/rarefaction.tsv", sep="\t", index=False)

print(f"\nRarefaction (family n={n}):")
for thr in THRESHOLDS:
    g = rare[rare.identity_pct == int(thr * 100)]
    big = g[g.N >= 20]
    gamma = np.polyfit(np.log(big.N), np.log(big.mean_clusters), 1)[0]
    e = elbow(g.N, g.mean_clusters)
    print(f"  clusters @{int(thr*100)}% identity: {g.mean_clusters.iloc[-1]:.0f} at N={n}; "
          f"Heaps exponent gamma={gamma:.2f}; elbow near N={int(g.N.iloc[e])}")
print(g.round(1).to_string(index=False) if False else "")
print(rare.round(1).to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
x = np.arange(len(fa))
ax[0].plot(x, fa.frac_features, marker="o", label="triad intact + Asp at S1")
ax[0].plot(x, fa.frac_triad, marker="s", label="triad intact")
ax[0].set_xticks(x)
ax[0].set_xticklabels([str(b) for b in fa.bin], rotation=60, fontsize=7)
ax[0].set_xlabel("Identity to PRSS1 (%)")
ax[0].set_ylabel("Fraction of hits")
ax[0].legend(fontsize=8)
ax[0].set_title("Feature retention vs identity")
for thr in THRESHOLDS:
    g = rare[rare.identity_pct == int(thr * 100)]
    ax[1].errorbar(g.N, g.mean_clusters, yerr=g.sd, marker="o", label=f"clusters @{int(thr*100)}% id")
ax[1].set_xlabel("Sequences sampled (N)")
ax[1].set_ylabel("Distinct clusters")
ax[1].legend(fontsize=8)
ax[1].set_title("Rarefaction")
plt.tight_layout()
plt.savefig("figures/fig2_elbow_analysis.png", dpi=200)
print("\nsaved figures/fig2_elbow_analysis.png, results/feature_retention.tsv, results/rarefaction.tsv")

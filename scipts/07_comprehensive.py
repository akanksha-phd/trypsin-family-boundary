## Step 7: comprehensive-mode design tests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO

FAMILY = "trypsin-like"
IDENTITIES = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
KS = [10, 20, 50, 100, 200]
N_DRAW = 50
MIN_OVERLAP = 100
MARKERS = {"His63": (63, "H"), "Asp107": (107, "D"), "Ser200": (200, "S"),
           "S1 (194)": (194, "D"), "Arg122": (122, "R")}

t = pd.read_csv("results/triad_check.tsv", sep="\t")
ok = t.triad_intact if FAMILY == "active-S1" else (t.triad_intact & (t.S1_res == "D"))
fam = t[ok].copy()
aln = {r.id: np.array(list(str(r.seq))) for r in AlignIO.read("results/domains.aln.fasta", "fasta")}
ids = ["PRSS1_domain"] + fam.target.tolist()
X = np.vstack([aln[i] for i in ids])
n = len(ids)
print(f"family '{FAMILY}': {n} sequences (including seed)")
is_trypsin = np.array([True] + (fam["class"] == "trypsin").tolist())

# distance matrix: 1 - identity over aligned columns
D = np.ones((n, n))
for i in range(n):
    both = (X[i] != "-") & (X != "-")
    same = (X[i] == X) & both
    ov = both.sum(axis=1)
    D[i] = np.where(ov >= MIN_OVERLAP, 1 - same.sum(axis=1) / np.maximum(ov, 1), 1.0)
np.fill_diagonal(D, 0.0)
D = np.minimum(D, D.T)


def greedy_clusters(d):
    near = D <= d
    left = np.ones(n, bool)
    reps, label = [], -np.ones(n, int)
    while left.any():
        counts = (near & left[None, :]).sum(axis=1)
        counts[~left] = -1
        r = int(np.argmax(counts))
        members = near[r] & left
        label[members] = len(reps)
        reps.append(r)
        left &= ~members
    return reps, label


# (a) cluster counts
rows_a, cluster_reps = [], {}
for ident in IDENTITIES:
    reps, _ = greedy_clusters(1 - ident)
    cluster_reps[ident] = reps
    rows_a.append((int(ident * 100), len(reps)))
ca = pd.DataFrame(rows_a, columns=["identity_pct", "n_clusters"])
ca.to_csv("results/cluster_counts.tsv", sep="\t", index=False)
print("\nClusters by identity threshold:")
print(ca.to_string(index=False))

# (b) conservation preservation
ref = X[0]
ref_cols = np.where(ref != "-")[0]
letters = list("ACDEFGHIKLMNPQRSTVWY")
code = {a: i for i, a in enumerate(letters)}
enc = np.vectorize(lambda c: code.get(c, -1))(X[:, ref_cols])


def conservation(rows):
    sub = enc[rows]
    out = np.full(sub.shape[1], np.nan)
    for j in range(sub.shape[1]):
        col = sub[:, j]
        col = col[col >= 0]
        if len(col) < 5:
            continue
        p = np.bincount(col, minlength=20) / len(col)
        p = p[p > 0]
        out[j] = 1 - (-(p * np.log2(p)).sum()) / np.log2(20)
    return out


def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m], b[m])[0, 1])


full = conservation(np.arange(n))
rng = np.random.default_rng(0)
medoid = int(np.argmin(D.mean(axis=1)))


def farthest(k):
    sel, nearest = [medoid], D[:, medoid].copy()
    while len(sel) < k:
        j = int(np.argmax(nearest))
        sel.append(j)
        nearest = np.minimum(nearest, D[:, j])
    return sel


rows_b = []
for k in KS:
    if k >= n:
        continue
    r = np.mean([corr(full, conservation(rng.choice(n, k, replace=False))) for _ in range(N_DRAW)])
    rows_b.append(("random", k, r))
    rows_b.append(("farthest-point", k, corr(full, conservation(farthest(k)))))
for ident, reps in cluster_reps.items():
    k = len(reps)
    if k < 5 or k >= n:
        continue
    rows_b.append((f"cluster reps @{int(ident*100)}%", k, corr(full, conservation(reps))))
    r = np.mean([corr(full, conservation(rng.choice(n, k, replace=False))) for _ in range(N_DRAW)])
    rows_b.append((f"random, same N as @{int(ident*100)}%", k, r))
cb = pd.DataFrame(rows_b, columns=["strategy", "n_seqs", "r_vs_full_conservation"])
cb.to_csv("results/subsample_stability.tsv", sep="\t", index=False)
print("\nConservation profile preserved (Pearson r vs full family):")
print(cb.round(3).to_string(index=False))

# (c) conservation track and autolysis-site composition
full_seq = "".join(l.strip() for l in open("data/PRSS1.fasta") if not l.startswith(">"))
offset = full_seq.index("IVGG")
tryp = conservation(np.where(is_trypsin)[0])
xpos = np.arange(len(ref_cols)) + offset + 1
fig, ax = plt.subplots(figsize=(11, 3.8))
ax.plot(xpos, full, label=f"{FAMILY} family (n={n})", lw=1.2)
ax.plot(xpos, tryp, label=f"trypsin-named only (n={int(is_trypsin.sum())})", lw=1.2)
print("\nMarker checks and residue composition at each marked column:")
for name, (pos, aa) in MARKERS.items():
    k = pos - 1 - offset
    seed_aa = full_seq[pos - 1]
    if seed_aa != aa:
        print(f"  {name}: seed has {seed_aa}, expected {aa} -> marker skipped")
        continue
    ax.axvline(pos, color="grey", ls=":", lw=0.8)
    ax.text(pos, 1.02, name, rotation=90, va="bottom", ha="center", fontsize=8)
    for lab, rows in [("family", np.arange(n)), ("trypsin", np.where(is_trypsin)[0])]:
        col = enc[rows, k]
        col = col[col >= 0]
        top = pd.Series([letters[c] for c in col]).value_counts(normalize=True).head(3)
        print(f"  {name} [{lab}]: " + ", ".join(f"{a} {v:.0%}" for a, v in top.items()))
ax.set_ylim(0, 1.15)
ax.set_xlabel("PRSS1 precursor residue number")
ax.set_ylabel("Conservation (1 - entropy/max)")
ax.legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.savefig("figures/fig5_conservation_track.png", dpi=200)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(ca.identity_pct, ca.n_clusters, marker="o")
ax[0].invert_xaxis()
ax[0].set_yscale("log")
ax[0].set_xlabel("Clustering identity threshold (%)")
ax[0].set_ylabel("Number of clusters")
ax[0].set_title(f"Family redundancy (n={n})")
for s, g in cb.groupby("strategy"):
    ax[1].plot(g.n_seqs, g.r_vs_full_conservation, marker="o", label=s)
ax[1].set_xscale("log")
ax[1].set_xlabel("Sequences in subsample")
ax[1].set_ylabel("r vs full-family conservation")
ax[1].legend(fontsize=6)
plt.tight_layout()
plt.savefig("figures/fig4_clusters_and_stability.png", dpi=200)
print("\nsaved figures/fig4_clusters_and_stability.png, fig5_conservation_track.png")

# Step 2: homolog counts at different identity and E-value thresholds, plus figure 1.

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

cols = "query target pident qcov tcov evalue bits qstart qend tstart tend".split()
d = pd.read_csv("results/hits.tsv", sep="\t", names=cols)
print("total hits:", len(d))

# Drop partial alignments so short fragment matches do not count as homologs
d = d[d.qcov >= 0.7]
print("after qcov >= 0.7:", len(d))

ids = [95, 90, 80, 70, 60, 50, 40, 35, 30, 25, 20]
counts = [int((d.pident >= i).sum()) for i in ids]
for i, c in zip(ids, counts):
    print(f"identity >= {i}%: {c}")
for e in [1e-100, 1e-50, 1e-30, 1e-20, 1e-10, 1e-5]:
    print(f"evalue <= {e}: {int((d.evalue <= e).sum())}")

pd.DataFrame({"min_identity": ids, "n_homologs": counts}).to_csv(
    "results/count_vs_identity.tsv", sep="\t", index=False)

plt.figure(figsize=(6, 4))
plt.plot(ids, counts, marker="o")
plt.gca().invert_xaxis()
plt.xlabel("Minimum % identity to PRSS1 domain")
plt.ylabel("Homologs in Swiss-Prot (PF00089)")
plt.title("Homolog count vs identity threshold")
plt.tight_layout()
plt.savefig("figures/fig1_count_vs_identity.png", dpi=200)
print("saved figures/fig1_count_vs_identity.png")

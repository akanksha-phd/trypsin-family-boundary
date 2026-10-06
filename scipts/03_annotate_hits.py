# Step 3: attach protein names, organism and a functional class to every hit.

import re
import pandas as pd

rows = []
for line in open("data/trypsin_swissprot.fasta"):
    if not line.startswith(">"):
        continue
    h = line[1:].strip()
    parts = h.split("|")
    acc, rest = parts[1], parts[2]
    entry = rest.split()[0]
    desc = rest.split(" ", 1)[1].split(" OS=")[0]
    m = re.search(r"OS=(.*?) OX=", h)
    rows.append((acc, entry, desc, m.group(1) if m else ""))
meta = pd.DataFrame(rows, columns=["target", "entry", "desc", "organism"])


def classify(desc):
    s = desc.lower()
    if "neurotrypsin" in s:
        return "neurotrypsin"
    if "elastase" in s:
        return "elastase"
    if "chymotrypsin" in s or "chymase" in s:
        return "chymotrypsin-like"
    if "tryptase" in s or "mastin" in s:
        return "tryptase"
    if "kallikrein" in s or "prostate-specific" in s or "tonin" in s:
        return "kallikrein"
    if "venom" in s or "thrombin-like" in s or "fibrinogenase" in s:
        return "venom protease"
    if any(w in s for w in ["coagulation", "complement", "plasminogen", "prothrombin", "factor"]):
        return "coag/complement"
    if re.search(r"\btrypsin\b|trypsinogen", s) or s == "serine protease 1":
        return "trypsin"
    return "other S1"


meta["class"] = meta.desc.map(classify)

cols = "query target pident qcov tcov evalue bits qstart qend tstart tend".split()
d = pd.read_csv("results/hits.tsv", sep="\t", names=cols)
d = d[(d.qcov >= 0.7) & (d.evalue <= 1e-5)].merge(meta, on="target")
d["bin"] = pd.cut(d.pident, [0, 30, 40, 50, 60, 70, 80, 90, 101])
print(len(d), "hits kept (qcov >= 0.7, evalue <= 1e-5)")
print(pd.crosstab(d["bin"], d["class"]))
d.to_csv("results/hits_annotated.tsv", sep="\t", index=False)

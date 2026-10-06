# Step 5: check the catalytic triad and the S1 specificity-pocket residue in every aligned hit.

import pandas as pd
import requests
from Bio import AlignIO

aln = {r.id: str(r.seq) for r in AlignIO.read("results/domains.aln.fasta", "fasta")}
ref = aln["PRSS1_domain"]
ref_nogap = ref.replace("-", "")

j = requests.get("https://rest.uniprot.org/uniprotkb/P07477.json", timeout=60).json()
active = sorted(f["location"]["start"]["value"] for f in j["features"] if f["type"] == "Active site")
full = "".join(l.strip() for l in open("data/PRSS1.fasta") if not l.startswith(">"))
offset = full.index("IVGG")
idx = [p - 1 - offset for p in active]  # 0-based positions in the mature domain
s1 = idx[-1] - 6

print("UniProt active sites (precursor numbering):", active)
print("residues at those sites:", [ref_nogap[i] for i in idx], "| S1 residue:", ref_nogap[s1])
print("self-check: expect H D S and S1 = D")

# Map domain positions to alignment columns
colmap, k = {}, -1
for c, ch in enumerate(ref):
    if ch != "-":
        k += 1
        colmap[k] = c
cols_ = [colmap[i] for i in idx + [s1]]

rows = [(acc, *[s[c] for c in cols_]) for acc, s in aln.items() if acc != "PRSS1_domain"]
t = pd.DataFrame(rows, columns=["target", "His", "Asp", "Ser", "S1_res"])
t["triad_intact"] = (t.His == "H") & (t.Asp == "D") & (t.Ser == "S")

meta = pd.read_csv("results/hits_annotated.tsv", sep="\t")[["target", "pident", "desc", "organism", "class"]]
t = t.merge(meta, on="target")

print("\nTriad intact, by class:")
print(pd.crosstab(t["class"], t.triad_intact))
print("\nS1 residue, by class:")
print(pd.crosstab(t["class"], t.S1_res))
t.to_csv("results/triad_check.tsv", sep="\t", index=False)

#!/usr/bin/env bash
# Step 4: trim each hit to its aligned protease domain and align all domains with MAFFT

set -euo pipefail
cd "$(dirname "$0")/.."

python - <<'EOF'
import pandas as pd
from Bio import SeqIO

h = pd.read_csv("results/hits_annotated.tsv", sep="\t")

seqs = {r.id.split("|")[1]: str(r.seq) for r in SeqIO.parse("data/trypsin_swissprot.fasta", "fasta")}
q = str(next(SeqIO.parse("data/PRSS1_domain.fasta", "fasta")).seq)
with open("data/domains.fasta", "w") as o:
    o.write(">PRSS1_domain\n" + q + "\n")
    for r in h.itertuples():
        o.write(f">{r.target}\n{seqs[r.target][r.tstart - 1:r.tend]}\n")
print(len(h), "domains written to data/domains.fasta")
EOF

mafft --auto --thread 4 data/domains.fasta > results/domains.aln.fasta
echo "alignment written to results/domains.aln.fasta"

# END
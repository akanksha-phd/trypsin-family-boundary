#!/usr/bin/env bash
# Step 1: download the seed and the Pfam PF00089 Swiss-Prot reference set,
# trim the seed to the mature protease domain, and search with MMseqs2


set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data results tmp figures

# Seed: human PRSS1 (UniProt P07477)
curl -sS -o data/PRSS1.fasta "https://rest.uniprot.org/uniprotkb/P07477.fasta"



# Reference universe: reviewed (Swiss-Prot) entries carrying the Pfam trypsin domain PF00089
curl -sS -o data/trypsin_swissprot.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?query=(xref:pfam-PF00089)+AND+(reviewed:true)&format=fasta"
echo "Reference sequences: $(grep -c '>' data/trypsin_swissprot.fasta)"


# Trim the seed to the mature enzyme. The activation peptide ends with DDDDK,
# and the mature chain starts at IVGG. Check this against the UniProt feature table.
python - <<'EOF'
seq = "".join(l.strip() for l in open("data/PRSS1.fasta") if not l.startswith(">"))
i = seq.index("IVGG")
dom = seq[i:]
open("data/PRSS1_domain.fasta", "w").write(">PRSS1_domain\n" + dom + "\n")
print("precursor length:", len(seq), "| mature domain length:", len(dom))
EOF

# search so the whole decay from close relatives to noise is visible.
# --max-seqs 1000: the MMseqs2 default keeps only 300 hits per query, which would
# silently cap every count below.
mmseqs easy-search data/PRSS1_domain.fasta data/trypsin_swissprot.fasta results/hits.tsv tmp \
  --format-output "query,target,pident,qcov,tcov,evalue,bits,qstart,qend,tstart,tend" \
  -e 10 -s 7.5 --max-seqs 1000 -v 1

echo "Total hits: $(wc -l < results/hits.tsv)"

## END

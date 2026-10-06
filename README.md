# trypsin-family-boundary
Empirical family-boundary analysis for trypsin

It tests how to define an empirical, defensible boundary for the trypsin family, as input to "quick and dirty" and "comprehensive" homolog sampling strategies.

Data
Seed: human PRSS1 (UniProt P07477), mature protease domain (starts at IVGG, 224 aa).
Reference set: 899 reviewed Swiss-Prot entries carrying the Pfam trypsin domain PF00089. This is a small curated set, not the full natural family.
Run order

Environment: conda env create -f environment.yml && conda activate trypsin

Step	Script	Output
1	scripts/01_get_data_and_search.sh	data/, results/hits.tsv
2	scripts/02_count_vs_identity.py	figures/fig1_count_vs_identity.png
3	scripts/03_annotate_hits.py	results/hits_annotated.tsv
4	scripts/04_trim_and_align.sh	results/domains.aln.fasta
5	scripts/05_triad_check.py	results/triad_check.tsv

Notes and caveats
MMseqs2 keeps 300 hits per query by default; step 1 sets --max-seqs 1000 so counts are not capped.
Hits are filtered to query coverage >= 0.7 and e-value <= 1e-5 before counting homologs.
Class labels are derived from UniProt protein names and are approximate.
We set domain boundaries for the seed at the activation site, not from a Pfam/InterPro annotation.

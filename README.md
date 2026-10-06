# Trypsin family boundary analysis

The project tests how to define an empirical, defensible boundary for the trypsin family, and compares "quick and dirty" and
"comprehensive" strategies for sampling homologs. 

## Data
- Seed: human PRSS1 (UniProt P07477), mature protease domain (starts at `IVGG`, 224 aa).
- Reference set: 899 reviewed Swiss-Prot entries carrying the Pfam trypsin domain PF00089.
  This is a small curated set, not the full natural family.

## Run order
Environment: `conda env create -f environment.yml && conda activate trypsin`
Run every step from the repository root (the scripts locate their own folder).

| Step | Command | What it does | Main outputs |
|---|---|---|---|
| 1 | `bash scripts/01_get_data_search.sh` | Downloads the seed and reference set, trims the seed, searches with MMseqs2 | `data/`, `results/hits.tsv` |
| 2 | `python scripts/02_count_identity.py` | Homolog counts by identity and E-value | `figures/fig1_count_vs_identity.png`, `results/count_vs_identity.tsv` |
| 3 | `python scripts/03_annotate_hits.py` | Adds protein names, organism and a functional class to each hit | `results/hits_annotated.tsv` |
| 4 | `bash scripts/04_trim_and_align.sh` | Trims hits to their protease domain and aligns them with MAFFT | `data/domains.fasta`, `results/domains.aln.fasta` |
| 5 | `python scripts/05_triad_check.py` | Checks the catalytic triad and the S1 pocket residue in every hit | `results/triad_check.tsv` |
| 6 | `python scripts/06_quick_selection.py` | Compares random, taxonomy-stratified and farthest-point selection | `figures/fig3_strategy_comparison.png`, `results/quick_selection.tsv` |
| 7 | `python scripts/07_comprehensive.py` | Cluster counts, subsample stability, conservation track | `figures/fig4_clusters_and_stability.png`, `figures/fig5_conservation_track.png`, `results/cluster_counts.tsv`, `results/subsample_stability.tsv` |
| 8 | `python scripts/08_elbow_analysis.py` | Feature retention vs identity and pangenome-style rarefaction | `figures/fig2_elbow_analysis.png`, `results/feature_retention.tsv`, `results/rarefaction.tsv` |

## Notes and caveats
- MMseqs2 keeps 300 hits per query by default; step 1 sets `--max-seqs 1000` so counts are not capped.
- Hits are filtered to query coverage >= 0.7 and E-value <= 1e-5 before analysis (796 of 899 remain).
- Family tiers: "active S1" = intact catalytic triad; "trypsin-like" = intact triad plus Asp at the S1 position.
- Functional classes come from UniProt protein names and are approximate.
- Residue numbers in step 7 use PRSS1 precursor numbering (247 aa).
- The seed's domain boundary was set at the activation site, not from a Pfam/InterPro annotation.
- Results demonstrate the method on a curated reference set; they are not a real absolute homolog count.

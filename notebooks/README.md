# Notebooks

The notebooks that ran the computational experiments, as executed in Google Colab, with outputs stripped to keep the repository small. Their outputs are the archives in `../results/`.

| Notebook | Experiments | Manuscript |
|---|---|---|
| `Experiments_01-06B_Robustness_R1-R3.ipynb` | 1–4 (equilibria, bifurcations, fold map, monostable memory), 5 (recovery–rechallenge), 6 and 6B (dynamic lag vs equilibrium hysteresis; quasistatic convergence), R1–R3 (robustness: α–ε map, multivariate sample, functional form) | Same-signal effect without bistability, Figures S1–S5 |
| `Experiment_07_Reset_Ratchet.ipynb` | 07 | Protocol 1; Figure 2, S7–S9 |
| `Experiment_08_Tipping_Window.ipynb` | 08 | Protocol 2 (first version; superseded by 14B) |
| `Experiment_09_Structural_Memory.ipynb` | 09 | Protocol 3; Figure 5 |
| `Experiment_10_Daily_Home_Care.ipynb` | 10 | Daily brushing vs smooth plaque curve; Figure S10 |
| `Experiment_11_Sign_of_Adaptation.ipynb` | 11 | Protocol 4; Figure 4A–C, S13 |
| `Experiment_12_Episodic_Activity.ipynb` | 12 | Protocol 5; Figure S17 |
| `Experiment_13_Resolving_Niche.ipynb` | 13 | Resolving niche; Figure S16 |
| `Experiment_14_Revision.ipynb` | 14A–14F | Protocols 6 and 7, defined therapy, alternative ratchet criteria, activity matching; Figures 3, 4D–F, S6, S11 |

The first notebook was run top to bottom as one file; its later sections (6B, R2, R3) use objects defined in earlier ones, so it is kept whole rather than split.

**Still to add before the first release:** `Paper_Figures.ipynb` and `Revision_Figures.ipynb`, which regenerate the figures from the exported data (named in the Data Availability Statement).

# Notebooks

The notebooks that ran the computational experiments, as executed in Google Colab, with outputs stripped to keep the repository small. Their outputs are the archives in `../results/`. The last column gives the figures of the JPR paper that draw on each experiment (SI Table S7); experiments without an entry belong to the development of the model or to an earlier version of the manuscript and are kept for the record.

| Notebook | Experiments | JPR paper |
|---|---|---|
| `Experiments_01-06B_Robustness_R1-R3.ipynb` | 1–4 (equilibria, bifurcations, fold map, monostable memory), 5 (recovery–rechallenge), 6 and 6B (dynamic lag vs equilibrium hysteresis; quasistatic convergence), R1–R3 (robustness: α–ε map, multivariate sample, functional form) | S2 (Experiment 5); S3 and Table S3 (R1–R3) |
| `Experiment_07_Reset_Ratchet.ipynb` | 07 | — |
| `Experiment_08_Tipping_Window.ipynb` | 08 | — (superseded by 14B) |
| `Experiment_09_Structural_Memory.ipynb` | 09 | — |
| `Experiment_10_Daily_Home_Care.ipynb` | 10 | — |
| `Experiment_11_Sign_of_Adaptation.ipynb` | 11 | — (its repeated-debridement panel is superseded by 19B) |
| `Experiment_12_Episodic_Activity.ipynb` | 12 | — |
| `Experiment_13_Resolving_Niche.ipynb` | 13 | — |
| `Experiment_14_Revision.ipynb` | 14A–14F: rechallenge, defined therapy and window, ratchet criteria, activity matching, 300-set sample, dose-dependent check | 1D, 2A–D, 2F, 5A; S1F, S4, S5, S6 |
| `Experiment_15_Model_Revision.ipynb` | 15A–15G: intrinsically bistable site and sample, calibrated and matched rechallenge, safe challenges and stopping rules, single vs repeated instrumentation, Weibull functions | 2E, 2F, 4C, 5B; S4, S6, S8B |
| `Experiment_16_Revision_Figures.ipynb` | 16A: sign of adaptation with the defined therapy; 16B: panel data for revised figures | 3, 4C, 5B; S4–S6, S7A |
| `Experiment_17_Resolution_Community_Memory.ipynb` | 17A–17C: resolution rate, pro-resolving adjunct, community memory crossed with host memory | 4A, 4B, 4D–F, 5D; S9; Table S4 |
| `Experiment_19_Debridement_From_Treated_State.ipynb` | 19A–19C: debridement at every recall from treated and never-diseased sites, by sign of adaptation, and convergence over 300 visits | 5C; S7B, S8A |

The first notebook was run top to bottom as one file; its later sections (6B, R2, R3) use objects defined in earlier ones, so it is kept whole rather than split.

The figures themselves are drawn from the exports by the scripts in `../scripts/` (they replace the `Paper_Figures` and `Revision_Figures` notebooks of an earlier draft).

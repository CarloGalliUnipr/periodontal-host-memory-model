# Exported results

Each archive is the `…_EXPORT` folder written by the notebook of the same name, unchanged. Inside each are the data tables (CSV/NPZ), the parameters (`parameters.json`), the internal checks, a `README.txt` describing every file, and `software_versions.txt`. `SHA256SUMS.txt` lists checksums (`sha256sum -c SHA256SUMS.txt`).

| Archive | Content |
|---|---|
| `Experiment_06B_Quasistatic_Hysteresis_Convergence_EXPORT.zip` | Ramp hysteresis converging to the equilibrium loop as the ramp slows |
| `Experiment_07_Reset_Ratchet_EXPORT.zip` | Periodic regrowth and removal: conversion thresholds, reduced map, latency scaling, validation against LSODA |
| `Experiment_08_Tipping_Window_EXPORT.zip` | Critical challenge by history and recovery time (state-reset therapy) |
| `Experiment_09_Structural_Memory_EXPORT.zip` | Attachment loss: stage boundaries, natural history, treatment ladder, maintenance |
| `Experiment_10_Daily_Home_Care_EXPORT.zip` | Daily brushing nested in recall vs its daily mean |
| `Experiment_11_Sign_of_Adaptation_EXPORT.zip` | Sensitizing, tolerizing and dose-dependent hosts |
| `Experiment_12_Episodic_Activity_EXPORT.zip` | Random plaque peaks, onset without therapy, treatment on detection |
| `Experiment_13_Resolving_Niche_EXPORT.zip` | Niche partly resolving with inflammation |
| `Experiment_14_Revision_EXPORT.zip` | 14A rechallenge; 14B window with defined therapy and relapse ceiling; 14C ratchet criteria and passage time; 14D activity-matched clustering; 14E 300-set robustness sample; 14F dose-dependent monotonicity; `all_checks.csv` |
| `Robustness_R1_alpha_epsilon_memory_map_EXPORT.zip` | Same-signal effect over the α–ε plane |
| `Robustness_R2_multivariate_structural_sample_EXPORT.zip` | 1,000 sampled parameter sets, monostable subset |
| `Robustness_R3_functional_form_EXPORT.zip` | Hill vs Weibull response functions |

**Still to add before the first release:** the exports of Experiment_05 (recovery–rechallenge) and Experiment_06 (dynamic lag vs equilibrium hysteresis), which are in Google Drive as loose folders (`Biosemiotic_Perio/Experiment_05_Recovery_Rechallenge_Memory`, `…/Experiment_06_Dynamic_Lag_vs_Equilibrium_Hysteresis`). Zip each folder, add it here and update `SHA256SUMS.txt`.

GitHub's web uploader accepts files up to 25 MB; the largest archive here is 10.8 MB.

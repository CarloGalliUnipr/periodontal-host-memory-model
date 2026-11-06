# Exported results

Each archive is the `…_EXPORT` folder written by the notebook of the same name, unchanged. Inside each are the data tables (CSV/NPZ), the parameters (`parameters.json`), the internal checks, a `README.txt` describing every file, and `software_versions.txt`. `SHA256SUMS.txt` lists checksums (`sha256sum -c SHA256SUMS.txt`). To unpack the archives used by the figures of the paper, run `python -I ../scripts/prepare_figure_data.py`.

| Archive | Content |
|---|---|
| `Experiment_05_Recovery_Rechallenge_Memory_EXPORT.zip` | Same signal, different history: recovery of the threshold and rechallenge (files at the top level of the archive) |
| `Experiment_06_Dynamic_Lag_vs_Equilibrium_Hysteresis_EXPORT.zip` | Dynamic lag against equilibrium hysteresis |
| `Experiment_06B_Quasistatic_Hysteresis_Convergence_EXPORT.zip` | Ramp hysteresis converging to the equilibrium loop as the ramp slows |
| `Experiment_07_Reset_Ratchet_EXPORT.zip` | Periodic regrowth and removal: conversion thresholds, reduced map, latency scaling, validation against LSODA |
| `Experiment_08_Tipping_Window_EXPORT.zip` | Critical challenge by history and recovery time (state-reset therapy) |
| `Experiment_09_Structural_Memory_EXPORT.zip` | Attachment loss: stage boundaries, natural history, treatment ladder, maintenance |
| `Experiment_10_Daily_Home_Care_EXPORT.zip` | Daily brushing nested in recall vs its daily mean |
| `Experiment_11_Sign_of_Adaptation_EXPORT.zip` | Sensitizing, tolerizing and dose-dependent hosts |
| `Experiment_12_Episodic_Activity_EXPORT.zip` | Random plaque peaks, onset without therapy, treatment on detection |
| `Experiment_13_Resolving_Niche_EXPORT.zip` | Niche partly resolving with inflammation |
| `Experiment_14_Revision_EXPORT.zip` | 14A rechallenge; 14B window with defined therapy and relapse ceiling; 14C ratchet criteria and passage time; 14D activity-matched clustering; 14E 300-set robustness sample; 14F dose-dependent monotonicity; `all_checks.csv` |
| `Experiment_15_Model_Revision_EXPORT.zip` | 15A intrinsically bistable site; 15B 300-set sample around it; 15C calibrated experimental gingivitis; 15D matching on resting inflammation; 15E safe challenges and stopping rules; 15F single vs repeated instrumentation; 15G Weibull functions |
| `Experiment_16_Revision_Figures_EXPORT.zip` | 16A sign of adaptation with the defined therapy; `figures/` with the panel data of the revised figures |
| `Experiment_17_Resolution_Community_Memory_EXPORT.zip` | 17A resolution rate; 17B pro-resolving adjunct; 17C community memory crossed with host memory; reduction check against the package |
| `Experiment_19_Debridement_From_Treated_State_EXPORT.zip` | 19A reference site and 19B sign of adaptation under debridement at every recall, from treated and never-diseased sites; 19C convergence over 300 visits; `numbers_for_text.csv`; `all_checks.csv` |
| `Robustness_R1_alpha_epsilon_memory_map_EXPORT.zip` | Same-signal effect over the α–ε plane |
| `Robustness_R2_multivariate_structural_sample_EXPORT.zip` | 1,000 sampled parameter sets, monostable subset |
| `Robustness_R3_functional_form_EXPORT.zip` | Hill vs Weibull response functions |

GitHub's web uploader accepts files up to 25 MB; the largest archive here is 14.2 MB.

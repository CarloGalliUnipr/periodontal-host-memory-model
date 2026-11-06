# Changelog

## 1.1.0 — 2026-11-06

Release accompanying the version of the manuscript submitted to the *Journal of Periodontal Research*. The model package is unchanged (`periomemory` 1.0.0), so every notebook and test of 1.0.0 gives the same results.

- Notebooks and exports of Experiments 15–17 and 19:
  - 15: intrinsically bistable site and the sample around it, calibrated and matched rechallenge, safe challenges and stopping rules, single versus repeated instrumentation, Weibull functions;
  - 16: sign of adaptation with the defined therapy, and panel data for the revised figures;
  - 17: resolution rate, pro-resolving adjunct, and community memory crossed with host memory;
  - 19: debridement at every recall from treated sites.
- Experiment 19 replaces the repeated-debridement analyses of Experiment 15F (part F2) and Experiment 11 (`debridement_by_sign.csv`) in Figures 5C, S7B and S8A. Those analyses started from a never-diseased site (H = M = A = 0), whereas the paper describes treated sites. Experiment 19 starts from the treated state; from the never-diseased start it reproduces the earlier tables (largest difference about 1e-14).
- `scripts/`: figure scripts for the paper (`make_figures_JPR.py`, `make_figures_JPR_SI.py`, `make_figure_S22.py`, `style.py`), `prepare_figure_data.py`, and the numbers ledger `numbers_ledger_JPR.csv`.
- `results/`: the exports of Experiments 05, 15, 16, 17 and 19 added to `SHA256SUMS.txt` (the archives of Experiments 05 and 06 were already present but not listed).
- README, notebook and results indexes, and `CITATION.cff` updated for the JPR paper.

## 1.0.0 — 2026-10-20

Release accompanying the submitted manuscript.

- Model with adaptive response state A, attachment loss D, generalized target A*(H) = αH − βH² and resolving niche (c).
- Protocols: defined therapy (debridement + 4 units of plaque control), critical challenge, periodic plaque removal with conversion / inflamed-recall / damaging-recall criteria, reduced reset-to-reset map, episodic activity with treatment on detection.
- Notebooks for Experiments 01–14 and Robustness R1–R3; exports of 06B, 07–14 and R1–R3.
- Regression tests against exported values; reference trajectories for the interactive lab.

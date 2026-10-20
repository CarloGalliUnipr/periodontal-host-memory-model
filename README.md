# Periodontal host-memory model

Code, notebooks and exported data for

> Galli C. **Host-response memory under periodic plaque removal: a dynamical model of post-therapy vulnerability and relapse in periodontitis.** *2026*

[![DOI](https://zenodo.org/badge/1427527316.svg)](https://doi.org/10.5281/zenodo.23513019)


This repository is the **archived scientific record** of the paper. The interactive teaching site, where the model can be explored in the browser, lives separately at **[periodontal-dynamics-lab](https://github.com/CarloGalliUnipr/periodontal-dynamics-lab)**. Results should be cited from this repository (by DOI and version), not from the interactive site.

> **Not a clinical tool.** The model is phenomenological and uncalibrated; its variables are dimensionless and the time mapping (1 unit ≈ 1.7 days) is illustrative. It does not support choosing recall intervals or predicting patient outcomes.

## Contents

| Folder | What it holds |
|---|---|
| `src/periomemory/` | Python package: the final model (H, M, A, D), equilibria and folds, and the paper's protocols (therapy, critical challenge, periodic plaque removal, reduced map, episodic activity) |
| `tests/` | Regression tests: the package reproduces values exported by the notebooks (about 30 s) |
| `notebooks/` | The Google Colab notebooks that ran every experiment, outputs stripped ([index](notebooks/README.md)) |
| `results/` | The exported data of each experiment, as produced by the notebooks ([index](results/README.md)) |
| `reference/` | Reference trajectories for checking other implementations, such as the JavaScript engine of the lab |
| `examples/` | A short script reproducing the logic of Figure 3A |

## The model

A single periodontal site with host inflammatory activity *H*, subgingival community activity *M*, adaptive response state *A* and attachment loss *D*, forced by supragingival plaque *u*:

```
dH/dt = -H + F(M + u; θ0·e^(-A))
dM/dt = r·(-M + G(H) + γ·D·(c + (1 - c)·H))
dA/dt = ε·(A*(H) - A),        A*(H) = α·H - β·H²
dD/dt = δ·[H - hD]+ · (1 - D)

F(S; θ) = Sⁿ / (Sⁿ + θⁿ),     G(H) = H^q / (H^q + K^q)
```

The **instant host** sets *A* = *A\**(*H*) at all times. It shares every equilibrium and fold with the memory host, so any difference between the two under the same forcing is an effect of memory. Reference set: n = q = 4, K = 0.9, θ0 = 0.8, r = 1, α = 1. Integration: classical RK4, step 0.1.

## Quick start

```bash
git clone https://github.com/CarloGalliUnipr/periodontal-host-memory-model.git
cd periodontal-host-memory-model
pip install -e ".[test]"
pytest                       # reproduces exported values
python examples/quickstart.py
```

```python
import numpy as np
import periomemory as pm

p = pm.REFERENCE
print(pm.folds(p))                     # (H_up, u_up, H_down, u_down)
print(pm.stage_boundaries(p))          # D_irr, D_spont under perfect plaque control

m = pm.Model(**p, eps=0.01)            # memory time 1/eps
diseased = pm.tile(pm.eq_state(p, 0.20, "high"), 1)
treated, ok, retained = pm.therapy(m, diseased)   # debridement + ~1 week of plaque control
print(ok, retained)                    # H back to health, ~97% of A retained
```

![Peak response to a 21-day challenge after therapy](docs/quickstart_rechallenge.png)

## Reproducing the paper

Each figure and table of the paper comes from the notebooks in `notebooks/`, run in Google Colab (Python 3.13; versions in each export's `software_versions.txt`). Each notebook writes its outputs to a folder named `…_EXPORT` and zips it; those zips are in `results/`, unchanged. The package in `src/` is a cleaned copy of the model and protocol code of `Experiment_14_Revision`, the last notebook, which itself reproduced the thresholds and critical challenges of Experiments 07–13 (see `all_checks.csv` in its export). `tests/` verifies that the package returns the same numbers.

To rerun an experiment, open its notebook in Colab (or Jupyter, skipping any Google Drive mount cell) and run all cells. Experiment_14_Revision has a `FAST = True` switch for a smoke test.

## Citation

Please cite the article and this software (see `CITATION.cff`; GitHub shows a "Cite this repository" button). Each release is archived on Zenodo with its own DOI.

## Licence

Code and notebooks: MIT (see `LICENSE`). Exported data in `results/` and `reference/`: CC BY 4.0 (see `LICENSE-DATA.md`).

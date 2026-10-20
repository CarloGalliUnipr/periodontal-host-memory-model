# Reference trajectories

`reference_trajectories.json` holds trajectories computed by `periomemory` (regenerate with `python scripts/export_reference_trajectories.py`), for checking any other implementation of the model, in particular the JavaScript engine (`docs/model.mjs`) of the interactive [Periodontal Dynamics Lab](https://github.com/CarloGalliUnipr/periodontal-dynamics-lab).

Cases: periodic plaque regrowth and removal (T = 60) for a memory host at ε = 0.03 (peak 0.62) and ε = 0.01 (peak 0.575, the slow plateau of Figure 2A), the instant host, and therapy (debridement + 4 units of plaque control) followed by baseline plaque. Each case lists parameters, forcing, initial state and H, M, A every unit.

An implementation with the same RK4 step (0.1) should agree to rounding error. The lab engine bounds H and M to [0, 1] for display; this does not affect these cases, in which both stay inside that interval.

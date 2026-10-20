"""Quick start: a treated site and a naive site receive the same 21-day experimental-gingivitis
challenge at increasing times after therapy (reference set, baseline plaque 0.20, plateau 0.40).

    python examples/quickstart.py   ->  quickstart_rechallenge.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import periomemory as pm

REF = pm.REFERENCE
ub, U = 0.20, 0.40
rec_days = np.array([0, 15, 30, 60, 90, 120, 180, 270, 365], float)
N = len(rec_days)
t_on = pm.align(pm.days_to_units(rec_days))
dur = pm.align(pm.days_to_units(21.0))
u_fun = pm.schedule(ub, t_on, dur, np.full(N, U), ka=pm.K_ACCESS)


def peak_response(model, state):
    peak = np.zeros(N)
    def mon(t, H, M, A, D):
        win = (t > t_on + 1e-9) & (t <= t_on + dur + 1e-9)
        np.maximum(peak, np.where(win, H, 0.0), out=peak)
    model.run(state, u_fun, 0.0, float(t_on.max() + dur), monitor=mon)
    return peak


fig, ax = plt.subplots(figsize=(6, 4))
naive = pm.tile(pm.eq_state(REF, ub, "low"), N)
diseased = pm.tile(pm.eq_state(REF, ub, "high"), N)
for label, eps, instant, style in [("instant host", 0.03, True, "o--"),
                                   ("memory, eps = 0.03", 0.03, False, "o-"),
                                   ("memory, eps = 0.01", 0.01, False, "s-")]:
    m = pm.Model(**REF, eps=eps, instant=instant)
    treated, ok, _ = pm.therapy(m, diseased)      # debridement + ~1 week of plaque control
    ax.plot(rec_days, peak_response(m, treated), style, label=f"treated, {label}")
m = pm.Model(**REF, eps=0.01)
ax.plot(rec_days, peak_response(m, naive), "k:", label="naive, eps = 0.01")
ax.set_xlabel("days since therapy (illustrative scale)")
ax.set_ylabel("peak H during a 21-day challenge")
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("quickstart_rechallenge.png", dpi=150)
print("wrote quickstart_rechallenge.png")

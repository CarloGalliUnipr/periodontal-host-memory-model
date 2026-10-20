"""Export reference trajectories for cross-checking other implementations (e.g. the JavaScript
engine of the interactive Periodontal Dynamics Lab).

    python scripts/export_reference_trajectories.py  ->  reference/reference_trajectories.json

Each case gives the parameters, the forcing, the initial state and H, M, A sampled every unit.
An independent implementation with the same RK4 step (0.1) should agree to ~1e-12; with a
different step or solver, to the tolerance of that solver.
"""

import json
import os

import numpy as np

import periomemory as pm

REF = pm.REFERENCE
OUT = os.path.join(os.path.dirname(__file__), "..", "reference", "reference_trajectories.json")


def periodic_case(name, eps, instant, peak, T, ncyc, state0):
    m = pm.Model(**REF, eps=eps, instant=instant)
    st = pm.tile(state0, 1)
    t_all, H, M, A = [], [], [], []
    for c in range(ncyc):
        tr = m.simulate(st, pm.regrowth(peak, T), 0.0, T, every=10)
        keep = slice(None) if c == 0 else slice(1, None)
        t_all += list(tr["t"][keep] + c * T); H += list(tr["H"][keep, 0]); M += list(tr["M"][keep, 0]); A += list(tr["A"][keep, 0])
        st, _ = m.run(st, pm.regrowth(peak, T), 0.0, T)
    return dict(name=name, params=dict(REF, eps=eps, instant=instant, b=0.0, gamma=0.0, delta=0.0),
                forcing=dict(kind="periodic_regrowth", peak=peak, T=T, k=pm.K_REGROW, cycles=ncyc,
                             note="u(tau) = peak/(1-exp(-k T)) * (1-exp(-k tau)), tau = t mod T; plaque removed at every multiple of T"),
                initial_state=dict(H=state0[0], M=state0[1], A=state0[2]),
                dt=pm.DT, t=t_all, H=H, M=M, A=A)


def therapy_case(eps):
    ub = 0.20
    m = pm.Model(**REF, eps=eps)
    high = pm.eq_state(REF, ub, "high", gamma=0.0)
    st = pm.tile(high, 1)
    st = (st[0], np.zeros(1), st[2], st[3])                      # debridement: M -> 0
    tr1 = m.simulate(st, pm.const(0.0), 0.0, pm.T_PC, every=10)   # perfect plaque control
    st, _ = m.run(st, pm.const(0.0), 0.0, pm.T_PC)
    tr2 = m.simulate(st, pm.const(ub), 0.0, 200.0, every=10)      # back to baseline
    t = list(tr1["t"]) + list(tr2["t"][1:] + pm.T_PC)
    return dict(name=f"therapy_then_baseline_eps{eps}", params=dict(REF, eps=eps, instant=False, b=0.0, gamma=0.0, delta=0.0),
                forcing=dict(kind="piecewise_constant", segments=[dict(t0=0.0, t1=pm.T_PC, u=0.0), dict(t0=pm.T_PC, t1=pm.T_PC + 200.0, u=ub)]),
                initial_state=dict(H=float(high[0]), M=0.0, A=float(high[2]),
                                   note="diseased equilibrium at u = 0.20 with M set to 0 (debridement)"),
                dt=pm.DT, t=t, H=list(tr1["H"][:, 0]) + list(tr2["H"][1:, 0]),
                M=list(tr1["M"][:, 0]) + list(tr2["M"][1:, 0]), A=list(tr1["A"][:, 0]) + list(tr2["A"][1:, 0]))


def main():
    f = pm.folds(REF)
    cases = [
        periodic_case("periodic_peak0.62_T60_eps0.03", 0.03, False, 0.62, 60.0, 5, (0.0, 0.0, 0.0, 0.0)),
        periodic_case("periodic_peak0.575_T60_eps0.01", 0.01, False, 0.575, 60.0, 5, (0.0, 0.0, 0.0, 0.0)),
        periodic_case("periodic_peak0.62_T60_instant", 0.03, True, 0.62, 60.0, 3, (0.0, 0.0, 0.0, 0.0)),
        therapy_case(0.03),
    ]
    doc = dict(
        package="periomemory", version=pm.__version__,
        equations=["dH/dt = -H + F(M+u; th0*exp(-A))", "dM/dt = r*(-M + G(H))", "dA/dt = eps*(a*H - A)",
                   "F(S;theta) = S^n/(S^n+theta^n)", "G(H) = H^q/(H^q+K^q)", "instant host: A = a*H at every RK4 stage"],
        integrator="classical RK4, fixed step dt = 0.1, no clipping of H or M",
        static_fold=dict(H_up=float(f[0]), u_up=float(f[1]), H_down=float(f[2]), u_down=float(f[3])),
        cases=cases)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(doc, fh, indent=1, default=float)
    print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    main()

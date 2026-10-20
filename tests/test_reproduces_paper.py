"""Regression tests: the package reproduces values exported by the experiment notebooks.

Reference values are copied from the exports in results/ (file and column given for each).
Tolerances are tight because the package uses the same integrator, step and bisection as
Experiment_14_Revision; differences larger than floating-point noise indicate a change in the model.
"""

import numpy as np
import pytest

import periomemory as pm
from periomemory import REFERENCE as REF


# ----------------------------------------------------------------- equilibria and stage boundaries
def test_static_fold_reference():
    # Experiment_07 summary.csv: "threshold at T=60, static fold (peak)" = 0.4682392215
    _, u_up, _, _ = pm.folds(REF)
    assert u_up == pytest.approx(0.4682392215, abs=1e-8)


def test_stage_boundaries():
    # Experiment_09 summary.csv: D_irr (analytic) 0.078775596, D_spont (analytic) 0.585299026875
    sb = pm.stage_boundaries(REF)
    assert sb["D_irr"] == pytest.approx(0.078775596, abs=1e-6)
    assert sb["D_spont"] == pytest.approx(0.585299026875, abs=1e-6)


def test_single_state_set_is_monostable():
    assert pm.folds(pm.SINGLE_STATE) is None


def test_memory_and_instant_hosts_share_equilibria():
    # Every equilibrium of the memory host satisfies A = A*(H): both hosts settle on the same state.
    for u, branch in [(0.20, "low"), (0.20, "high"), (0.40, "low")]:
        st = pm.tile(pm.eq_state(REF, u, branch, gamma=0.0), 1)
        for inst in (False, True):
            m = pm.Model(**REF, eps=0.03, instant=inst)
            out, _ = m.run(st, pm.const(u), 0.0, 400.0)
            assert out[0][0] == pytest.approx(st[0][0], abs=1e-6)
            assert out[2][0] == pytest.approx(st[2][0], abs=1e-6)


# ----------------------------------------------------------------- critical challenge (Exp 08 / 14B)
def test_critical_challenge_by_history():
    # 14B_tipping_window/B1_critical_vs_duration.csv, host eps=0.03, duration 2/eps
    ub, eps = 0.20, 0.03
    m = pm.Model(**REF, eps=eps)
    hists = ["naive", "treated_protocol", "treated_reset"]
    states = [pm.init_history(m, REF, ub, h, eps, 1, False) for h in hists]
    st = tuple(np.concatenate([s[k] for s in states]) for k in range(4))
    crit, _ = pm.critical_challenge(m, st, ub, np.zeros(3), np.full(3, 2 / eps), 10 / eps)
    expected = [0.564890, 0.294977, 0.272769]
    assert crit == pytest.approx(expected, abs=2e-6)


def test_critical_challenge_instant_host():
    # B1_critical_vs_duration.csv: instant, naive, duration 20 -> 0.47634888
    m = pm.Model(**REF, eps=0.1, instant=True)
    st = pm.init_history(m, REF, 0.20, "naive", 0.1, 1, True)
    crit, _ = pm.critical_challenge(m, st, 0.20, np.zeros(1), np.full(1, 20.0), 200.0)
    assert crit[0] == pytest.approx(0.47634888, abs=2e-6)


# ----------------------------------------------------------------- rechallenge after therapy (14A)
def test_rechallenge_after_therapy():
    # 14A_rechallenge/rechallenge_reference.csv: baseline 0.20, U = 0.40, host eps=0.01
    ub, U, eps = 0.20, 0.40, 0.01
    rec_days = np.array([0.0, 60.0, 90.0])
    m = pm.Model(**REF, eps=eps, gamma=pm.GAMMA)
    N = len(rec_days)
    t_on = pm.align(pm.days_to_units(rec_days)); dur = pm.align(pm.days_to_units(21.0))
    u_fun = pm.schedule(ub, t_on, dur, np.full(N, U), ka=pm.K_ACCESS)

    def peak_H(state):
        peak = np.zeros(N)
        def mon(t, H, M, A, D):
            win = (t > t_on + 1e-9) & (t <= t_on + dur + 1e-9)
            np.maximum(peak, np.where(win, H, 0.0), out=peak)
        m.run(state, u_fun, 0.0, float(t_on.max() + dur), monitor=mon)
        return peak

    naive = pm.tile(pm.eq_state(REF, ub, "low", D=0.0), N)
    treated, ok, retained = pm.therapy(m, pm.tile(pm.eq_state(REF, ub, "high", D=0.0), N))
    assert ok.all()
    assert retained[0] == pytest.approx(0.970567, abs=1e-6)
    assert peak_H(naive) == pytest.approx([0.060207] * 3, abs=1e-6)
    assert peak_H(treated) == pytest.approx([0.989615, 0.948759, 0.527179], abs=1e-6)


# ----------------------------------------------------------------- relapse ceiling (14B4)
def test_relapse_ceiling_eps001():
    # 14B_tipping_window/B4_relapse_ceiling_summary.csv: eps=0.01, treated_protocol -> 0.265
    eps = 0.01
    m = pm.Model(**REF, eps=eps)
    ubs = np.round(np.arange(0.07, 0.4651, 0.005), 3)
    H, M, A = [], [], []
    for u in ubs:
        s, ok, _ = pm.therapy(m, pm.tile(pm.eq_state(REF, u, "high", gamma=0.0), 1))
        H.append(s[0][0]); M.append(s[1][0]); A.append(s[2][0])
    st = (np.array(H), np.array(M), np.array(A), np.zeros(len(ubs)))
    out, _ = m.run(st, pm.const(ubs), 0.0, 10 / eps)
    rel = out[0] > 0.5
    ceiling = next(ubs[i] for i in range(len(ubs)) if rel[i:].all())
    assert ceiling == pytest.approx(0.265, abs=1e-9)


# ----------------------------------------------------------------- periodic plaque removal (Exp 07 / 14C)
@pytest.mark.parametrize("host, eps, instant, expected", [
    ("instant", 0.1, True, 0.50),
    ("eps=0.01", 0.01, False, 0.58),
])
def test_ratchet_threshold_T60(host, eps, instant, expected):
    # 14C_ratchet/C1_thresholds.csv and Experiment_07 summary.csv, T = 60
    T = 60.0
    peaks = np.round(np.arange(0.30, 0.9001, 0.01), 3)
    ncyc = 20 if instant else int(max(40, np.ceil(40 / (eps * T))))
    assert pm.ratchet_threshold(REF, eps, instant, T, peaks, ncyc) == pytest.approx(expected, abs=1e-9)


def test_ratchet_inflamed_recall_precedes_conversion():
    # 14C_ratchet/C1_grid.csv, eps=0.01, T=60, peak 0.58: first inflamed recall 8, first damaging 13
    conv, inf, dam = pm.ratchet(REF, 0.01, False, 60.0, [0.58], 67)
    assert inf[0] == 8 and dam[0] == 13
    assert conv[0] > inf[0]


# ----------------------------------------------------------------- dose-dependent host (14F)
def test_dose_dependent_host():
    # 14F_dose_dependent_check/dose_dependent_monotonicity.csv: naive 0.36, treated_reset 0.38
    a, b, ub, eps = 1.0, 1.4, 0.25, 0.03
    grid = np.round(np.arange(0.26, 4.0001, 0.01), 3); N = len(grid)
    m = pm.Model(**pm.SECOND, b=b, eps=eps)
    lo = pm.eq_state(pm.SECOND, ub, "low", gamma=0.0, a=a, b=b)
    hi = pm.eq_state(pm.SECOND, ub, "high", gamma=0.0, a=a, b=b)
    firsts = []
    for A0 in (lo[2], hi[2]):
        st = pm.tile((lo[0], lo[1], A0, 0.0), N)
        st, _ = m.run(st, pm.const(grid), 0.0, 2 / eps)
        st, _ = m.run(st, pm.const(ub), 0.0, 10 / eps)
        conv = st[0] > 0.5
        assert conv[np.argmax(conv):].all()          # monotone in challenge size
        firsts.append(grid[np.argmax(conv)])
    assert firsts == pytest.approx([0.36, 0.38], abs=1e-9)

"""Protocols used in the paper: therapy, critical challenge, periodic plaque removal (ratchet),
reduced reset-to-reset map and episodic activity.

These functions are taken from Experiment_14_Revision (the final revision analyses), which
reproduced the conversion thresholds and critical challenges of Experiments 07-13.
"""

import numpy as np

from .equilibria import eq_state, tile
from .model import Model, const, schedule, align
from .parameters import DT, GAMMA, H_D, K_REGROW, T_PC

TRAPZ = getattr(np, "trapezoid", None) or np.trapz


# ----------------------------------------------------------------- therapy
def therapy(model, start_state, t_pc=T_PC):
    """Debridement (M -> 0) followed by perfect supragingival plaque control (u = 0) for t_pc.

    Returns (state at the end, success flag H < 0.1, fraction of A retained).
    """
    H, M, A, D = [np.array(x, float) for x in start_state]
    A0 = A.copy()
    M = np.zeros_like(M)
    (H, M, A, D), _ = model.run((H, M, A, D), const(0.0), 0.0, t_pc)
    with np.errstate(invalid="ignore", divide="ignore"):
        retained = np.where(A0 != 0, A / A0, np.nan)
    return (H, M, A, D), H < 0.1, retained


def init_history(model, p, ub, history, eps, N, instant, t_pc=None):
    """Initial states at baseline plaque ub after one of the paper's histories.

    'naive': low equilibrium. 'gingivitis': u = 0.45 for 5/eps from naive (resolved gingivitis).
    'treated_protocol': diseased equilibrium, then therapy(). 'treated_reset': healthy H, M with
    the A of the diseased state (Experiment 08 definition).
    """
    low = eq_state(p, ub, "low", gamma=0.0)
    if history == "naive":
        return tile(low, N)
    if history == "gingivitis":
        st, _ = model.run(tile(low, N), const(0.45), 0.0, 5 / eps if not instant else 100.0)
        return st
    high = eq_state(p, ub, "high", gamma=0.0)
    if history == "treated_reset":
        return tile((low[0], low[1], high[2], 0.0), N)
    if history == "treated_protocol":
        st, ok, _ = therapy(model, tile(high, N), T_PC if t_pc is None else t_pc)
        if not np.all(ok):
            raise RuntimeError("therapy failed")
        return st
    raise ValueError(history)


# ----------------------------------------------------------------- bisection on challenge size
def bisect(converts, lo, hi, n_iter):
    """Vectorized bisection on amplitude. converts(U) -> bool array.

    Returns (critical value, nan if hi does not convert; flag for cases already converting at lo).
    """
    lo, hi = np.array(lo, float), np.array(hi, float)
    c_hi = converts(hi); c_lo = converts(lo)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        c = converts(mid)
        hi = np.where(c, mid, hi); lo = np.where(c, lo, mid)
    crit = np.where(c_hi, hi, np.nan)
    crit = np.where(c_lo, np.asarray(lo, float) * 0 + np.nan, crit)
    return crit, c_lo


def critical_challenge(model, states, ub, t_on, dur, settle, hi=4.0, n_iter=16):
    """Smallest step challenge (from baseline ub, starting at t_on, lasting dur) that leaves H > 0.5
    after `settle` units back at baseline."""
    t_on, dur = align(t_on), align(dur)
    t_end = float((t_on + dur).max() + settle)
    def converts(U):
        st, _ = model.run(states, schedule(ub, t_on, dur, U), 0.0, t_end)
        return st[0] > 0.5
    N = len(t_on)
    return bisect(converts, np.full(N, ub), np.full(N, hi), n_iter)


def half_life(x, y):
    """First x at which y (starting at 1) falls to 0.5, by linear interpolation."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    for i in range(1, len(y)):
        if y[i] <= 0.5 and y[i - 1] > 0.5:
            return x[i-1] + (y[i-1] - 0.5) * (x[i] - x[i-1]) / (y[i-1] - y[i])
    return np.nan


# ----------------------------------------------------------------- periodic plaque removal
def ratchet(p, eps, instant, T, peaks, ncyc, k=K_REGROW):
    """Naive sites under periodic regrowth and removal (interval T, given peaks), followed for ncyc cycles.

    Returns, per peak: conversion cycle (first cycle from which every later cycle is entrapped,
    H > 0.5 throughout), first inflamed recall (H > 0.5 at removal) and first damaging recall (H > 0.8).
    Cycles are counted from 1; nan if never.
    """
    peaks = np.asarray(peaks, float); N = len(peaks)
    umax = peaks / (1 - np.exp(-k * T))
    m = Model(**p, eps=eps, instant=instant)
    st = tuple(np.zeros(N) for _ in range(4))
    u_fun = lambda t: umax * (1 - np.exp(-k * t))
    ent = np.zeros((ncyc, N), bool); hrec = np.zeros((ncyc, N))
    for c in range(ncyc):
        hmin = np.full(N, 9.0)
        def mon(t, H, M, A, D): np.minimum(hmin, H, out=hmin)
        st, _ = m.run(st, u_fun, 0.0, T, monitor=mon)
        ent[c] = hmin > 0.5; hrec[c] = st[0]
    conv = np.full(N, np.nan)
    tail = np.flip(np.logical_and.accumulate(np.flip(ent, 0), 0), 0)
    for i in range(N):
        if tail[-1, i]:
            conv[i] = np.argmax(tail[:, i]) + 1
    def first(mask):
        out = np.full(N, np.nan)
        for i in range(N):
            if mask[:, i].any(): out[i] = np.argmax(mask[:, i]) + 1
        return out
    return conv, first(hrec > 0.5), first(hrec > 0.8)


def ratchet_threshold(p, eps, instant, T, peaks, ncyc, k=K_REGROW):
    """Lowest peak on the grid above which every peak converts (the tolerated-peak threshold)."""
    peaks = np.asarray(peaks, float)
    conv, _, _ = ratchet(p, eps, instant, T, peaks, ncyc, k=k)
    ok = ~np.isnan(conv)
    for i in range(len(ok)):
        if ok[i:].all():
            return peaks[i]
    return np.nan


def theta_map(p, T, peak, A_grid, ncyc=30, k=K_REGROW):
    """Reduced reset-to-reset map: mean H over one interval (last of ncyc) with memory frozen at A,
    from a naive start, and the entrapped flag. a * Theta(A) is the cycle-mean adaptation target."""
    N = len(A_grid); umax = peak / (1 - np.exp(-k * T))
    m = Model(**p, eps=0.0)
    st = (np.zeros(N), np.zeros(N), np.array(A_grid, float), np.zeros(N))
    u_fun = lambda t: umax * (1 - np.exp(-k * t))
    for c in range(ncyc):
        acc = np.zeros(N); hmin = np.full(N, 9.0)
        def mon(t, H, M, A, D):
            acc[:] += H * DT; np.minimum(hmin, H, out=hmin)
        st, _ = m.run(st, u_fun, 0.0, T, monitor=mon)
    return acc / T, hmin > 0.5


# ----------------------------------------------------------------- variable plaque, episodic activity
def episodic(p, eps, instant, peaks, therapy_on=True, T=60.0, gamma=GAMMA, delta=0.001, hD=H_D, k=K_REGROW):
    """Sites x intervals of random plaque peaks (array N x I) with treatment on detection
    (M -> 0 when H > 0.5 at recall). Returns (attachment loss per interval, treated flags, final D)."""
    peaks = np.asarray(peaks, float); N, I = peaks.shape
    m = Model(**p, eps=eps, instant=instant, gamma=gamma, delta=delta, hD=hD)
    st = tuple(np.zeros(N) for _ in range(4))
    loss = np.zeros((N, I)); treated = np.zeros((N, I), bool)
    for i in range(I):
        umax = peaks[:, i] / (1 - np.exp(-k * T))
        D0 = st[3].copy()
        st, _ = m.run(st, lambda t: umax * (1 - np.exp(-k * t)), 0.0, T)
        loss[:, i] = st[3] - D0
        if therapy_on:
            tr = st[0] > 0.5; treated[:, i] = tr
            st = (st[0], np.where(tr, 0.0, st[1]), st[2], st[3])
    return loss, treated, st[3]

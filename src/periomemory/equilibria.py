"""Equilibria, folds and stage boundaries.

At equilibrium H = F(M + u; th0 e^{-A*(H)}) and M = G(H), so the forcing that holds
a given H is

    u(H) = th0 * exp(-A*(H)) * (H / (1 - H))**(1/n) - G(H).

The memory host and the instant host satisfy the same relation and therefore
share every equilibrium and fold. With attachment loss D held fixed, equilibria
are those of the model at effective forcing u + gamma*D.
"""

import numpy as np
from scipy.optimize import brentq

from .model import hill
from .parameters import GAMMA

#: Dense grid of H on (0, 1) used to locate folds and bracket roots.
HG = np.linspace(1e-7, 1 - 1e-7, 200001)


def u_curve(p, H=HG, a=None, b=0.0):
    """Forcing u at which H is an equilibrium (parameter dict p)."""
    a = p["a"] if a is None else a
    return p["th0"] * np.exp(-(a * H - b * H * H)) * (H / (1 - H)) ** (1 / p["n"]) - hill(H, p["K"], p["q"])


def folds(p, a=None, b=0.0):
    """(H_up, u_up, H_down, u_down): loss of the low state (local max of u) and of the high state
    (local min of u). None if there is no fold (monostable)."""
    u = u_curve(p, a=a, b=b); du = np.diff(u)
    idx = np.where(np.sign(du[1:]) != np.sign(du[:-1]))[0] + 1
    if len(idx) < 2:
        return None
    i_up, i_dn = idx[0], idx[1]
    return HG[i_up], u[i_up], HG[i_dn], u[i_dn]


def eq_H(p, u, branch, a=None, b=0.0):
    """Equilibrium H at constant forcing u on the 'low' or 'high' branch (None if absent)."""
    a = p["a"] if a is None else a
    f = lambda H: u_curve(p, np.array([H]), a=a, b=b)[0] - u
    v = u_curve(p, a=a, b=b) - u
    i = np.where(np.sign(v[1:]) != np.sign(v[:-1]))[0]
    roots = [brentq(f, HG[j], HG[j + 1]) for j in i]
    if not roots:
        return None
    return roots[0] if branch == "low" else roots[-1]


def eq_state(p, u, branch, D=0.0, gamma=GAMMA, a=None, b=0.0):
    """Equilibrium (H, M, A, D) with D fixed (niche fully persistent): forcing u + gamma*D."""
    a = p["a"] if a is None else a
    H = eq_H(p, u + gamma * D, branch, a=a, b=b)
    if H is None:
        return None
    return H, hill(H, p["K"], p["q"]) + gamma * D, a * H - b * H * H, D


def stage_boundaries(p, gamma=GAMMA):
    """Attachment-loss boundaries under perfect plaque control (u = 0).

    D_irr = u_down / gamma: above it hygiene alone cannot restore a diseased site.
    D_spont = u_up / gamma: above it no low-inflammatory equilibrium remains.
    """
    f = folds(p)
    if f is None:
        return None
    _, u_up, _, u_dn = f
    return dict(D_irr=max(u_dn, 0.0) / gamma, D_spont=u_up / gamma, u_up=u_up, u_down=u_dn)


def tile(state, N):
    """Repeat a scalar state N times (for vectorized runs)."""
    return tuple(np.full(N, float(x)) for x in state)

"""Host-response memory model of a periodontal site.

State variables (all vectorized over sites):

    H  host inflammatory activity
    M  functional activity of the subgingival community
    A  adaptive response state (sets the response threshold theta = th0 * exp(-A))
    D  attachment loss (structural memory)

Equations::

    dH/dt = -H + F(M + u; th0 * exp(-A))
    dM/dt = r * (-M + G(H) + gamma * D * (c + (1 - c) * H))
    dA/dt = eps * (A*(H) - A),          A*(H) = a*H - b*H**2
    dD/dt = delta * max(H - hD, 0) * (1 - D)

    F(S; theta) = S**n / (S**n + theta**n),   G(H) = H**q / (H**q + K**q)

`u(t)` is supragingival plaque. The instant host sets A = A*(H) at every
integration stage; eps = 0 freezes A. Integration is fixed-step RK4 (dt = 0.1),
exactly as in the experiment notebooks (Experiment_14_Revision, cell "vectorized model").
"""

import numpy as np

from .parameters import DT, H_D


def hill(x, K, p):
    """Hill function x**p / (x**p + K**p), with x clipped at 0."""
    x = np.maximum(x, 0.0)
    xp = x ** p
    return xp / (xp + K ** p)


class Model:
    """H, M, A, D on vectors of sites. Parameters may be scalars or arrays of the vector length.

    Parameters
    ----------
    n, q, K, th0, r, a : response steepness, community steepness, community half-saturation,
        constitutive threshold, community rate, adaptation strength (alpha).
    b : quadratic term of the target A*(H) = a H - b H^2 (dose-dependent host; 0 by default).
    eps : memory rate (1/eps is the memory time).
    instant : if True, A follows A*(H) without delay (memory-free control with identical equilibria).
    gamma, c : support of the community by attachment loss, gamma*D*(c + (1-c)H); c = 1 is a fixed niche.
    delta, hD : rate and inflammation threshold of attachment loss.
    dt : RK4 step.
    """

    def __init__(self, n, q, K, th0, r, a, b=0.0, eps=0.03, instant=False,
                 gamma=0.0, c=1.0, delta=0.0, hD=H_D, dt=DT):
        self.n, self.q, self.K, self.th0, self.r, self.a, self.b = [
            np.asarray(v, float) for v in (n, q, K, th0, r, a, b)]
        self.eps, self.gamma, self.c, self.delta, self.hD = [
            np.asarray(v, float) for v in (eps, gamma, c, delta, hD)]
        self.instant = instant
        self.dt = float(dt)

    def target(self, H):
        """Target adaptive state A*(H) = a H - b H^2."""
        return self.a * H - self.b * H * H

    def threshold(self, A):
        """Response threshold theta = th0 * exp(-A)."""
        return self.th0 * np.exp(-A)

    def deriv(self, H, M, A, D, u):
        if self.instant:
            A = self.target(H)
        th = self.th0 * np.exp(-A)
        dH = -H + hill(M + u, th, self.n)
        dM = self.r * (-M + hill(H, self.K, self.q) + self.gamma * D * (self.c + (1.0 - self.c) * H))
        dA = np.zeros_like(H) if self.instant else self.eps * (self.target(H) - A)
        dD = self.delta * np.maximum(H - self.hD, 0.0) * (1.0 - D)
        return dH, dM, dA, dD

    def run(self, state, u_fun, t0, t1, monitor=None):
        """Integrate from t0 to t1 with forcing u_fun(t).

        state : tuple (H, M, A, D) of scalars or arrays.
        monitor : optional callback monitor(t, H, M, A, D) called after every step.
        Returns ((H, M, A, D), t_end).
        """
        DT_ = self.dt
        H, M, A, D = [np.array(x, float) for x in state]
        n = int(round((t1 - t0) / DT_)); t = t0
        for _ in range(n):
            u1, u2, u3 = u_fun(t), u_fun(t + DT_ / 2), u_fun(t + DT_)
            k1 = self.deriv(H, M, A, D, u1)
            k2 = self.deriv(H + DT_/2*k1[0], M + DT_/2*k1[1], A + DT_/2*k1[2], D + DT_/2*k1[3], u2)
            k3 = self.deriv(H + DT_/2*k2[0], M + DT_/2*k2[1], A + DT_/2*k2[2], D + DT_/2*k2[3], u2)
            k4 = self.deriv(H + DT_*k3[0], M + DT_*k3[1], A + DT_*k3[2], D + DT_*k3[3], u3)
            H = H + DT_/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
            M = M + DT_/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
            A = A + DT_/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
            D = D + DT_/6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3])
            t = t + DT_
            if self.instant:
                A = self.target(H)
            if monitor is not None:
                monitor(t, H, M, A, D)
        return (H, M, A, D), t

    def simulate(self, state, u_fun, t0, t1, every=1):
        """Integrate and return the trajectory as a dict of arrays (t, H, M, A, D, u, theta).

        `every` keeps one sample in `every` steps. Intended for plotting and cross-checks,
        not for large vectorized sweeps.
        """
        H0, M0, A0, D0 = [np.array(x, float) for x in state]
        rec = {k: [] for k in ("t", "H", "M", "A", "D", "u")}
        def keep(t, H, M, A, D):
            rec["t"].append(t); rec["H"].append(H.copy()); rec["M"].append(M.copy())
            rec["A"].append(A.copy()); rec["D"].append(D.copy()); rec["u"].append(np.asarray(u_fun(t), float).copy())
        A0_eff = self.target(H0) if self.instant else A0
        keep(t0, H0, M0, A0_eff, D0)
        count = [0]
        def mon(t, H, M, A, D):
            count[0] += 1
            if count[0] % every == 0:
                keep(t, H, M, A, D)
        self.run((H0, M0, A0, D0), u_fun, t0, t1, monitor=mon)
        out = {k: np.array(v) for k, v in rec.items()}
        out["theta"] = self.threshold(out["A"])
        return out


# ----------------------------------------------------------------- forcing helpers
def const(u):
    """Constant forcing."""
    return lambda t: u


def schedule(ub, t_on, dur, U, ka=None):
    """Baseline ub; from t_on for dur either a step to U (ka None) or plaque rising from ub towards U at rate ka."""
    ub, t_on, dur, U = [np.asarray(x, float) for x in (ub, t_on, dur, U)]
    def u(t):
        on = (t >= t_on - 1e-9) & (t < t_on + dur - 1e-9)
        val = U if ka is None else ub + (U - ub) * (1.0 - np.exp(-ka * np.maximum(t - t_on, 0.0)))
        return np.where(on, val, ub)
    return u


def regrowth(peak, T, k=None):
    """Plaque regrowing within one recall interval, u(tau) = u_max (1 - exp(-k tau)), tau in [0, T).

    `peak` is the plaque reached just before removal; u_max = peak / (1 - exp(-k T)).
    Run one interval at a time with Model.run(state, regrowth(peak, T), 0.0, T): the
    removal at each recall is the restart of tau at 0.
    """
    from .parameters import K_REGROW
    k = K_REGROW if k is None else k
    umax = np.asarray(peak, float) / (1 - np.exp(-k * T))
    return lambda t: umax * (1 - np.exp(-k * t))


def align(x, dt=DT):
    """Round times to the integration grid."""
    return np.round(np.asarray(x, float) / dt) * dt

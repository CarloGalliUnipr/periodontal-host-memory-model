"""Supporting-Information Figures S1-S9 of the JPR manuscript.

Usage: python -I make_figures_JPR_SI.py <data folder> <output folder>
<data folder> holds E14, E15, E16, E17, E19 and r1, r2, r3 (unpacked exports; see prepare_figure_data.py).
Figures S1 and S2 are recomputed from the model (periomemory 1.0.0); the others are drawn from the exports.
Figures S7B and S8A (debridement at every recall of treated sites) are read from Experiment 19 (parts 19B and 19A).
Figure S9 is drawn by make_figure_S22.py (renamed); run it separately.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import periomemory as pm
from style import *

DATA, OUT = sys.argv[1], sys.argv[2]
P = lambda *a: os.path.join(DATA, *a)
F16 = P("E16", "figures")
setup()
EPSC = {"eps=0.1": "#9ec5f4", "eps=0.03": AQUA, "eps=0.01": BLUE, "eps=0.003": "#104281", "instant": ORANGE}


def days_axis(a):
    a.set_xticks([0, 90, 180, 270, 365]); a.set_xlabel("days after therapy")


# ============================================================== Figure S1: equilibria and the two hosts after therapy
def figureS1():
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.7)); ax = axs.ravel(); src = {}
    # A host alone
    a = ax[0]; H = np.linspace(0.001, 0.999, 2000); rows = []
    for al, col in [(0, GREY), (0.5, "#6da7ec"), (1, BLUE)]:
        S = 0.8 * np.exp(-al * H) * (H / (1 - H)) ** 0.25
        a.plot(S, H, color=col, label=f"α = {al:g}")
        rows.append(pd.DataFrame({"alpha": al, "H": H, "S": S}))
    a.set_xlim(0, 2); a.set_xlabel("signal S"); a.set_ylabel("steady-state inflammation H*")
    a.legend(loc="lower right"); panel(a, "A", "Host alone (n = 4)"); src["A"] = pd.concat(rows)
    # B full loop, reference set
    a = ax[1]; rows = []
    for al, col in [(0, GREY), (0.5, "#6da7ec"), (1, BLUE)]:
        u = pm.u_curve(pm.REFERENCE, a=al)
        Hg = pm.equilibria.HG
        du = np.gradient(u)
        stable = du > 0
        k = slice(None, None, 50)
        uu, hh, ss = u[k], Hg[k], stable[k]
        for st, ls in [(True, "-"), (False, "--")]:
            y = np.where(ss == st, hh, np.nan)
            a.plot(uu, y, ls, color=col, lw=1.2)
        rows.append(pd.DataFrame({"alpha": al, "H": hh, "u": uu, "stable": ss}))
    a.axvspan(-0.6, 0, color="#f4f3ef", lw=0)
    a.set_xlim(-0.6, 1.0); a.set_xlabel("plaque u"); a.set_ylabel("steady-state inflammation H*")
    a.legend(handles=[Line2D([], [], color=c, label=f"α = {v:g}") for v, c in [(0, GREY), (0.5, "#6da7ec"), (1, BLUE)]] +
             [Line2D([], [], color=INK2, ls="--", label="unstable")], loc="lower right")
    panel(a, "B", "Reference set, full loop"); src["B"] = pd.concat(rows)
    # C bistable range against alpha
    a = ax[2]; rows = []
    for name, p, col in [("reference", pm.REFERENCE, BLUE), ("single-state", pm.SINGLE_STATE, GREY), ("second", pm.SECOND, INK)]:
        als = np.round(np.arange(0, 2.01, 0.05), 3); w = []
        for al in als:
            f = pm.folds(p, a=al)
            if f is None or f[1] <= 0:
                w.append(0.0)
            else:
                w.append(max(f[1] - max(f[3], 0.0), 0.0))
            rows.append(dict(set=name, alpha=al, width=w[-1], u_up=None if f is None else f[1], u_down=None if f is None else f[3]))
        a.plot(als, w, color=col, label=name)
    a.set_xlabel("adaptation strength α"); a.set_ylabel("plaque range with two stable states")
    a.legend(loc="upper left"); panel(a, "C", "Bistability and adaptation"); src["C"] = pd.DataFrame(rows)
    # D, E therapy time courses
    p = pm.REFERENCE; ub = 0.2; rows = []
    for host, eps, inst in [("eps=0.01", 0.01, False), ("eps=0.03", 0.03, False), ("instant", 0.03, True)]:
        m = pm.Model(**p, eps=eps, instant=inst)
        st = pm.eq_state(p, ub, "high")
        H, M, A, D = [np.array([x], float) for x in st]
        M = M * 0
        pc = pm.align(4.0)
        u_fun = lambda t: np.where(t < pc, 0.0, ub)
        tr = m.simulate((H, M, A, D), u_fun, 0.0, pm.days_to_units(365.0), every=5)
        d = pm.units_to_days(tr["t"])
        ax[3].plot(d, tr["H"][:, 0], color=EPSC[host]); ax[4].plot(d, tr["A"][:, 0], color=EPSC[host])
        rows.append(pd.DataFrame({"host": host, "days": d, "H": tr["H"][:, 0], "A": tr["A"][:, 0], "M": tr["M"][:, 0]}))
    for a, lab, lt in [(ax[3], "inflammation H", "D"), (ax[4], "response state A", "E")]:
        a.set_xscale("symlog", linthresh=10); a.set_xticks([0, 5, 10, 30, 100, 365]); a.set_xticklabels(["0", "5", "10", "30", "100", "365"])
        a.set_xlabel("days after debridement"); a.set_ylabel(lab); a.axvspan(0, 4 * UNIT_DAYS, color="#eef4fc", lw=0)
    panel(ax[3], "D", "Inflammation after therapy"); panel(ax[4], "E", "Response state after therapy")
    ax[3].legend(handles=[Line2D([], [], color=BLUE, label="memory 167 d"), Line2D([], [], color=AQUA, label="memory 56 d"),
                          Line2D([], [], color=ORANGE, label="instant host")], loc="upper right")
    src["D_E"] = pd.concat(rows)
    # F memory retained after therapy (14A)
    a = ax[5]
    r = pd.read_csv(P("E14", "14A_rechallenge", "rechallenge_reference.csv"))
    g = r[(r.baseline == 0.2) & (r.history == "treated_D0") & (r.host != "instant")].groupby("host").A_retained_after_therapy.first()
    a.bar([0, 1], [g["eps=0.01"], g["eps=0.03"]], color=[BLUE, AQUA], width=0.55)
    for i, v in enumerate([g["eps=0.01"], g["eps=0.03"]]):
        a.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=6.5, color=INK)
    a.set_xticks([0, 1]); a.set_xticklabels(["memory 167 d", "memory 56 d"]); a.set_ylim(0, 1.1)
    a.set_ylabel("A retained after therapy"); panel(a, "F", "Memory kept by therapy")
    src["F"] = g.reset_index()
    fig.tight_layout(w_pad=1.0, h_pad=1.3)
    save(fig, OUT, "FigureS1", src)


# ============================================================== Figure S2: same signal, different history (single-state set)
def figureS2():
    p = pm.SINGLE_STATE; eps = 0.1; m = pm.Model(**p, eps=eps)
    st0 = pm.eq_state(p, 0.05, "low", gamma=0.0)
    Fh = lambda S, th: S ** p["n"] / (S ** p["n"] + th ** p["n"])
    trA = m.simulate(st0, pm.const(0.05), 0, 200); trB = m.simulate(st0, pm.const(0.40), 0, 200)
    sA = tuple(np.array(trA[k][-1]) for k in "HMAD"); sB = tuple(np.array(trB[k][-1]) for k in "HMAD")
    thA, thB = p["th0"] * np.exp(-sA[2]), p["th0"] * np.exp(-sB[2]); S = np.sqrt(thA * thB)
    rA = m.simulate(sA, pm.const(0.05), 200, 360); rB = m.simulate(sB, pm.const(0.05), 200, 360)
    fig, axs = plt.subplots(1, 4, figsize=(7.0, 2.2)); src = {}
    a = axs[0]; Sg = np.linspace(0, 2, 400)
    a.plot(Sg, Fh(Sg, thA), "--", color=GREY, label=f"history A, θ = {thA:.2f}")
    a.plot(Sg, Fh(Sg, thB), color=BLUE, label=f"history B, θ = {thB:.2f}")
    a.axvline(S, color=MUTED, lw=0.6, ls=":")
    a.annotate("", (S, Fh(S, thB)), (S, Fh(S, thA)), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8))
    a.text(S + 0.05, 0.5, f"ΔR = {Fh(S, thB) - Fh(S, thA):.2f}\nat S* = {S:.2f}", fontsize=6, color=INK)
    a.set_xlabel("signal S"); a.set_ylabel("host response F(S; θ)"); a.legend(loc="center right", bbox_to_anchor=(1.0, 0.42), fontsize=5.4)
    panel(a, "A", "Same signal"); src["A"] = pd.DataFrame({"S": Sg, "F_A": Fh(Sg, thA), "F_B": Fh(Sg, thB)})
    a = axs[1]
    t = np.concatenate([trA["t"], rA["t"]])
    thAt = p["th0"] * np.exp(-np.concatenate([trA["A"], rA["A"]])); thBt = p["th0"] * np.exp(-np.concatenate([trB["A"], rB["A"]]))
    a.plot(t, thAt, "--", color=GREY, label="A: u = 0.05"); a.plot(t, thBt, color=BLUE, label="B: u = 0.40, then 0.05")
    a.axvline(200, color=MUTED, lw=0.6, ls=":")
    a.set_xlabel("time (model units)"); a.set_ylabel("threshold θ"); a.legend(loc="center", bbox_to_anchor=(0.55, 0.62), fontsize=5.4)
    panel(a, "B", "Threshold"); src["B"] = pd.DataFrame({"t": t, "theta_A": thAt, "theta_B": thBt})
    a = axs[2]
    T = rA["t"] - 200
    dR = Fh(S, p["th0"] * np.exp(-rB["A"])) - Fh(S, p["th0"] * np.exp(-rA["A"]))
    dth = p["th0"] * np.exp(-rA["A"]) - p["th0"] * np.exp(-rB["A"])
    a.plot(T, dR / dR[0], color=BLUE, label="excess response ΔR")
    a.plot(T, dth / dth[0], color=AQUA, label="threshold deficit")
    a.plot(T, np.exp(-eps * T), ":", color=INK2, label="passive decay e^(−εt)")
    hl = T[np.argmax(dR <= dR[0] / 2)]; ht = T[np.argmax(dth <= dth[0] / 2)]
    a.set_xlim(0, 80); a.set_xlabel("recovery time (units)"); a.set_ylabel("fraction of initial value")
    a.legend(loc="center right", fontsize=5.4); panel(a, "C", "Recovery")
    src["C"] = pd.DataFrame({"T": T, "dR": dR, "dtheta": dth, "half_life_dR": hl, "half_life_theta": ht})
    a = axs[3]
    dist = np.sqrt(sum((rA[k] - rB[k]) ** 2 for k in "HMA"))
    a.semilogy(T, dist, color=INK)
    a.set_xlabel("recovery time (units)"); a.set_ylabel("distance between states"); panel(a, "D", "Same final state")
    src["D"] = pd.DataFrame({"T": T, "distance": dist})
    fig.tight_layout(w_pad=0.8)
    save(fig, OUT, "FigureS2", src)
    return dict(thA=thA, thB=thB, S=S, dR=float(dR[0]), half_dR=hl, half_th=ht, dist_end=float(dist[-1]))


# ============================================================== Figure S3: robustness of the same-signal effect
def figureS3():
    fig, axs = plt.subplots(1, 4, figsize=(7.0, 2.3)); src = {}
    r1 = pd.read_csv(P("r1", "figR1_01_same_signal_memory_map_data.csv"))
    h1 = pd.read_csv(P("r1", "figR1_02_memory_half_life_map_data.csv"))
    for a, df, col, lab, letter, ttl in [(axs[0], r1, "delta_R_protocol", "ΔR", "A", "Size"),
                                         (axs[1], h1, "apparent_memory_half_life", "half-life (units)", "B", "Persistence")]:
        d = df[df.stable_monostable.astype(bool)]
        piv = d.pivot_table(index="alpha", columns="epsilon", values=col)
        al = piv.index.values; ep = piv.columns.values
        im = a.pcolormesh(np.log10(ep), al, piv.values, cmap="Blues", shading="nearest")
        cb = fig.colorbar(im, ax=a, fraction=0.05, pad=0.03); cb.ax.tick_params(labelsize=6); cb.set_label(lab, fontsize=6.5)
        a.set_xlabel("log₁₀ ε"); a.set_ylabel("adaptation strength α"); panel(a, letter, ttl)
        src[letter] = d
    r2 = pd.read_csv(P("r2", "figR2_03_delta_R_vs_alpha_data.csv"))
    a = axs[2]
    sc = a.scatter(r2.alpha, r2.delta_R_same_signal, c=np.log10(r2.epsilon), cmap="Blues", s=5, lw=0)
    a.set_xlabel("adaptation strength α"); a.set_ylabel("ΔR"); panel(a, "C", "434 sampled sets")
    cb = fig.colorbar(sc, ax=a, fraction=0.05, pad=0.03); cb.ax.tick_params(labelsize=6); cb.set_label("log₁₀ ε", fontsize=6.5)
    src["C"] = r2
    r3 = pd.read_csv(P("r3", "figR3_05_common_subset_delta_R_data.csv"))
    w = r3.pivot_table(index="sample_id", columns="variant", values="delta_R_same_signal")
    a = axs[3]
    a.scatter(w["HillF_HillG"], w["WeibullF_WeibullG"], s=5, color=BLUE, lw=0)
    lim = [0, max(w.max()) * 1.05]; a.plot(lim, lim, ":", color=MUTED, lw=0.7)
    a.set_xlabel("ΔR, Hill functions"); a.set_ylabel("ΔR, Weibull functions"); panel(a, "D", "Hill vs Weibull")
    src["D"] = w.reset_index()
    fig.tight_layout(w_pad=0.6)
    save(fig, OUT, "FigureS3", src)


# ============================================================== Figure S4: robustness of the post-therapy window
def figureS4():
    ab = pd.read_csv(os.path.join(F16, "FigureS6", "FigureS6_A_B_data.csv"))
    e = pd.read_csv(os.path.join(F16, "FigureS6", "FigureS6_E_data.csv"))
    f = pd.read_csv(os.path.join(F16, "FigureS6", "FigureS6_F_data.csv"))
    ref = ab[ab["sample"] == "reference neighbourhood"]
    thd = ab[ab["sample"] == "third-set neighbourhood"]
    fig, axs = plt.subplots(1, 4, figsize=(7.0, 2.3))
    bins = np.linspace(0.3, 1.0, 36)
    for a, col_ref, col_th, letter, ttl, xl in [(axs[0], "ratio_treated_naive", "ratio_treated_naive", "A", "Treated / naive", "critical challenge, treated / naive"),
                                               (axs[1], "ratio_treated_instant", "ratio_treated_instant", "B", "Treated / instant host", "critical challenge, treated / instant")]:
        a.hist(ref[col_ref].dropna(), bins=bins, color=BLUE, alpha=0.85, label=f"around reference (n = {ref[col_ref].notna().sum()})")
        a.hist(thd[col_th].dropna(), bins=bins, color=VIOLET, alpha=0.75, label=f"around third set (n = {thd[col_th].notna().sum()})")
        a.axvline(1, color=MUTED, lw=0.6, ls=":")
        a.set_xlabel(xl); a.set_ylabel("parameter sets"); panel(a, letter, ttl)
    axs[0].legend(loc="upper left", fontsize=5.4)
    a = axs[2]
    sc = a.scatter(e.a, e.ratio_treated_naive, c=e.K, cmap="Blues", s=6, lw=0)
    cb = fig.colorbar(sc, ax=a, fraction=0.05, pad=0.03); cb.ax.tick_params(labelsize=6); cb.set_label("K", fontsize=6.5)
    a.set_xlabel("adaptation strength α"); a.set_ylabel("critical challenge, treated / naive"); panel(a, "C", "Third-set sample")
    a = axs[3]; x = np.arange(len(f)); wb = 0.36
    a.bar(x - wb / 2, f.hill, wb, color=BLUE, label="Hill"); a.bar(x + wb / 2, f.weibull, wb, color="#9ec5f4", label="Weibull")
    a.set_xticks(x); a.set_xticklabels(["ref.\n56 d", "ref.\n167 d", "third\n56 d", "third\n167 d"], fontsize=6)
    a.axhline(1, color=MUTED, lw=0.6, ls=":"); a.set_ylim(0, 1.1); a.set_ylabel("critical challenge, treated / naive")
    a.legend(loc="upper left", fontsize=5.6); panel(a, "D", "Hill vs Weibull")
    fig.tight_layout(w_pad=0.6)
    save(fig, OUT, "FigureS4", {"A_B": ab, "C": e, "D": f})


# ============================================================== Figure S5: after therapy, definitions and the size of the window
def figureS5():
    F3 = os.path.join(F16, "Figure3"); S11 = os.path.join(F16, "FigureS11")
    d = pd.read_csv(os.path.join(F3, "Figure3_D_data.csv")); v = pd.read_csv(os.path.join(F3, "Figure3_E_data.csv"))
    sa = pd.read_csv(os.path.join(S11, "FigureS11_A_data.csv")); sb = pd.read_csv(os.path.join(S11, "FigureS11_B_data.csv"))
    sc = pd.read_csv(os.path.join(S11, "FigureS11_C_data.csv")); sd = pd.read_csv(os.path.join(S11, "FigureS11_D_data.csv"))
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.7)); ax = axs.ravel()
    a = ax[0]
    for (h, hist), col, ls, lab in [(("eps=0.03", "naive"), GREY, "--", "naive"), (("eps=0.03", "gingivitis"), AQUA, ":", "resolved gingivitis"),
                                    (("eps=0.03", "treated_protocol"), AQUA, "-", "treated"), (("instant", "naive"), ORANGE, "-", "instant host")]:
        g = d[(d.host == h) & (d.history == hist)].sort_values("duration_days")
        y = g.critical if hist != "naive" or h == "instant" else g.critical
        if h == "instant":
            y = g.critical_instant if g.critical.isna().all() else g.critical
        ok = y.notna()
        a.plot(g.duration_days[ok], y[ok], ls, color=col, marker="o", ms=2.6, label=lab)
    a.axhline(FOLD_REF, color=MUTED, lw=0.6, ls=":"); a.axhline(0.2, color=INK2, lw=0.6)
    a.set_xscale("log"); a.set_xlabel("challenge duration (days)"); a.set_ylabel("critical challenge")
    a.legend(loc="center right", bbox_to_anchor=(1.0, 0.45), fontsize=5.4); panel(a, "A", "Duration (memory 56 d)")
    a = ax[1]
    for h in ["eps=0.1", "eps=0.03", "eps=0.01"]:
        g = v[(v.host == h)].sort_values("eps_T_rec")
        a.plot(g.eps_T_rec, g.V, "-o", color=EPSC[h], ms=2.6, label={"eps=0.1": "memory 17 d", "eps=0.03": "memory 56 d", "eps=0.01": "memory 167 d"}[h])
    x = np.linspace(0, 8, 200); a.plot(x, np.exp(-x), ":", color=INK2, label="passive decay")
    a.set_xlim(0, 6); a.set_xlabel("memory times since therapy"); a.set_ylabel("normalized vulnerability V")
    a.legend(loc="upper right", fontsize=5.6); panel(a, "B", "Recovery of vulnerability")
    a = ax[2]
    for h in ["eps=0.03", "eps=0.01"]:
        for hist, ls in [("treated_protocol", "-"), ("treated_reset", "--")]:
            g = sa[(sa.host == h) & (sa.history == hist)].sort_values("eps_T_rec")
            a.plot(g.eps_T_rec, g.critical, ls, color=EPSC[h])
        a.axhline(sa[sa.host == h].critical_naive.iloc[0], color=EPSC[h], lw=0.6, ls=":")
    a.set_xlabel("memory times since therapy"); a.set_ylabel("critical challenge, treated")
    a.legend(handles=[Line2D([], [], color=INK2, label="defined therapy"), Line2D([], [], color=INK2, ls="--", label="instantaneous reset"),
                      Line2D([], [], color=INK2, ls=":", lw=0.6, label="naive")], loc="lower right", fontsize=5.6)
    panel(a, "C", "Definition of therapy")
    a = ax[3]
    sb = sb.copy(); sb["eps"] = pd.to_numeric(sb.host.str[4:], errors="coerce")
    for kind, ls, lab in [("treated_protocol", "-", "defined therapy"), ("treated_reset", "--", "instantaneous reset")]:
        g = sb[(sb.kind == kind) & sb.eps.notna()].sort_values("eps")
        ok = g.ceiling.notna()
        a.plot(g.eps[ok], g.ceiling[ok], ls, color=BLUE, marker="o", ms=3, label=lab)
    a.axhline(FOLD_REF, color=ORANGE, lw=1, label="static threshold")
    a.set_xscale("log"); a.set_xlabel("memory rate ε"); a.set_ylabel("relapse ceiling (baseline plaque)")
    a.legend(loc="center left", fontsize=5.4); panel(a, "D", "Relapse ceiling")
    a = ax[4]
    fr = [0.5, 0.6, 0.7, 0.8, 0.9]
    for h in ["eps=0.1", "eps=0.03", "eps=0.01"]:
        g = sc[(sc.host == h) & (sc.history == "treated_protocol")].iloc[0]
        a.plot(fr, [g[f"window_at_{f}_naive_memory_units"] for f in fr], "-o", color=EPSC[h], ms=3)
    a.set_xlabel("relevant challenge, fraction of naive"); a.set_ylabel("window length (memory times)")
    a.legend(handles=[Line2D([], [], color=EPSC[h], label=l) for h, l in [("eps=0.1", "memory 17 d"), ("eps=0.03", "memory 56 d"), ("eps=0.01", "memory 167 d")]], loc="upper left", fontsize=5.4)
    panel(a, "E", "Length of the window")
    a = ax[5]
    for h in ["eps=0.1", "eps=0.03", "eps=0.01"]:
        g = sd[sd.host == h].sort_values("plaque_control_days")
        ok = g.therapy_success.astype(bool)
        a.plot(g.plaque_control_days[ok], g.critical_treated[ok], "-o", color=EPSC[h], ms=3)
        a.plot(g.plaque_control_days[~ok], g.critical_treated[~ok], "x", color=EPSC[h], ms=4.5)
    a.set_xscale("log"); a.set_xticks([3.3, 6.7, 16.7, 33.3, 66.7]); a.set_xticklabels(["3", "7", "17", "33", "67"]); a.minorticks_off()
    a.set_xlabel("plaque control (days)"); a.set_ylabel("critical challenge, treated")
    a.text(3.5, 0.56, "× inflammation not yet resolved", fontsize=5.4, color=INK2)
    panel(a, "F", "Length of plaque control")
    fig.tight_layout(w_pad=1.0, h_pad=1.3)
    save(fig, OUT, "FigureS5", {"A": d, "B": v, "C": sa, "D": sb, "E": sc, "F": sd})


# ============================================================== Figure S6: rechallenge design
def figureS6():
    F3 = os.path.join(F16, "Figure3"); S11 = os.path.join(F16, "FigureS11"); S12 = os.path.join(F16, "FigureS12")
    fig, axs = plt.subplots(3, 3, figsize=(7.0, 6.9)); ax = axs.ravel(); src = {}

    def calib(a, df, letter, title):
        for h in ["eps=0.01", "eps=0.03", "instant"]:
            g = df[(df.host == h) & (df.history == "treated")].sort_values("days")
            if g.empty:
                g = df[(df.host == h) & (df.history != "naive")].sort_values("days")
            a.plot(g.days, g.peak_H, "-o", color=EPSC[h], ms=2.8)
            pp = g[g.persistent.astype(bool)]
            a.plot(pp.days, pp.peak_H, "o", mfc="white", mec=EPSC[h], ms=4, mew=0.8)
        a.axhline(0.3, color=GREY, ls="--", lw=1)
        days_axis(a); a.set_ylim(0, 1.05); a.set_ylabel("peak inflammation H"); panel(a, letter, title)
    A = pd.read_csv(os.path.join(F3, "Figure3_A_data.csv")); calib(ax[0], A, "A", "Calibrated challenge, r = 1"); src["A"] = A
    E = pd.read_csv(os.path.join(S11, "FigureS11_E_data.csv")); calib(ax[1], E, "B", "Calibrated, r = 0.03"); src["B"] = E
    Fd = pd.read_csv(os.path.join(S11, "FigureS11_F_data.csv")); calib(ax[2], Fd, "C", "Third set, r = 0.03"); src["C"] = Fd
    ax[0].legend(handles=[Line2D([], [], color=BLUE, label="memory 167 d"), Line2D([], [], color=AQUA, label="memory 56 d"),
                          Line2D([], [], color=ORANGE, label="instant host"), Line2D([], [], color=GREY, ls="--", label="naive")],
                 loc="center right", fontsize=5.6)
    a = ax[3]
    B = pd.read_csv(os.path.join(S12, "FigureS12_D_data.csv")); src["D"] = B
    for (h, hist), ls in [(("eps=0.01", "treated_D0"), "-"), (("eps=0.03", "treated_D0"), "-"), (("eps=0.01", "treated_Dres"), "--"), (("eps=0.03", "treated_Dres"), "--")]:
        g = B[(B.host == h) & (B.history == hist)].sort_values("recovery_days")
        y = g.lowest_U_persistent.fillna(0.72)
        a.plot(g.recovery_days, y, ls, color=EPSC[h], marker="o", ms=2.8)
        c = g[g.lowest_U_persistent.isna()]
        a.plot(c.recovery_days, np.full(len(c), 0.72), "^", mfc="white", mec=EPSC[h], ms=4)
    a.axhline(0.55, color=ORANGE, lw=1)
    days_axis(a); a.set_ylim(0.2, 0.76); a.set_ylabel("smallest converting plateau")
    a.legend(handles=[Line2D([], [], color=INK2, label="D = 0"), Line2D([], [], color=INK2, ls="--", label="D = 0.10"),
                      Line2D([], [], color=ORANGE, label="instant host")], loc="upper right", fontsize=5.4)
    a.text(5, 0.21, "lines end where no plateau up to 0.70 converts", fontsize=5.2, color=INK2)
    panel(a, "D", "Conversion boundary")
    GH = pd.read_csv(os.path.join(S11, "FigureS11_G_H_data.csv")); src["E_F"] = GH
    for a, col, ylab, letter, ttl in [(ax[4], "U_safe", "largest plateau without conversion", "E", "Safe challenge"),
                                      (ax[5], "excess_at_U_safe", "peak, treated − naive", "F", "Signal at the safe plateau")]:
        for h in ["eps=0.03", "eps=0.01"]:
            for dur, ls in [(7, ":"), (14, "--"), (21, "-")]:
                g = GH[(GH.host == h) & (GH.duration_days == dur)].sort_values("days")
                a.plot(g.days, g[col], ls, color=EPSC[h])
                if col == "U_safe":
                    c = g[g.lowest_converting_U.isna()]
                    a.plot(c.days, c[col], "o", mfc="white", mec=EPSC[h], ms=3.2, mew=0.7)
        days_axis(a); a.set_ylabel(ylab); panel(a, letter, ttl)
    ax[4].legend(handles=[Line2D([], [], color=INK2, ls=s, label=f"{d}-day") for d, s in [(7, ":"), (14, "--"), (21, "-")]],
                 loc="lower right", fontsize=5.6)
    a = ax[6]
    Bx = pd.read_csv(os.path.join(S12, "FigureS12_B_data.csv")); src["G"] = Bx
    for h in ["eps=0.01", "eps=0.03", "instant"]:
        for hist, ls in [("treated_D0", "-"), ("treated_Dres", "--")]:
            g = Bx[(Bx.host == h) & (Bx.history == hist)].sort_values("recovery_days")
            a.plot(g.recovery_days, g.peak_H, ls, color=EPSC[h])
            pp = g[g.persistent.astype(bool)]
            a.plot(pp.recovery_days, pp.peak_H, "o", mfc="white", mec=EPSC[h], ms=3.6, mew=0.8)
    days_axis(a); a.set_ylim(0, 1.05); a.set_ylabel("peak inflammation H"); panel(a, "G", "Residual attachment loss")
    a.legend(handles=[Line2D([], [], color=INK2, label="D = 0"), Line2D([], [], color=INK2, ls="--", label="D = 0.10")], loc="center right", fontsize=5.4)
    s12 = pd.read_csv(os.path.join(S12, "FigureS12_E_data.csv")); s12f = pd.read_csv(os.path.join(S12, "FigureS12_F_data.csv"))
    a = ax[7]
    for eps, df in [(0.01, s12f), (0.03, s12)]:
        g = df[(df.eps == eps)].drop_duplicates("days").sort_values("days")
        a.plot(g.days, g.H_rest_treated, "-o", color=BLUE if eps == 0.01 else AQUA, ms=2.8)
    a.axhline(float(s12.H_rest_naive.iloc[0]), color=GREY, ls="--")
    days_axis(a); a.set_ylabel("resting inflammation H"); panel(a, "H", "Resting inflammation")
    a.legend(handles=[Line2D([], [], color=BLUE, label="treated, memory 167 d"), Line2D([], [], color=AQUA, label="treated, memory 56 d"), Line2D([], [], color=GREY, ls="--", label="naive")], loc="upper right", fontsize=5.4)
    src["H"] = pd.concat([s12, s12f])
    a = ax[8]
    for U, ls in [(0.4, "-")]:
        g = s12f[(s12f.U == U)].sort_values("days")
        a.plot(g.days, g.peak_treated, ls, color=BLUE, marker="o", ms=2.8, label="treated")
        a.plot(g.days, g.peak_naive_matched_same_increment, ls, color=INK, marker="s", ms=2.8, label="naive, matched H, same rise")
        a.plot(g.days, g.peak_naive_matched_same_U, ":", color=INK, label="naive, matched H, same plateau")
        a.plot(g.days, g.peak_naive_standard, "--", color=GREY, label="naive")
    days_axis(a); a.set_ylim(0, 1.05); a.set_ylabel("peak inflammation H"); a.legend(loc="center right", bbox_to_anchor=(1.0, 0.66), fontsize=5.0)
    panel(a, "I", "Matching (memory 167 d)")
    fig.tight_layout(w_pad=0.9, h_pad=1.2)
    save(fig, OUT, "FigureS6", src)


# ============================================================== Figure S7: sign of adaptation
def figureS7():
    S13 = os.path.join(F16, "FigureS13")
    a13 = pd.read_csv(os.path.join(S13, "FigureS13_A_data.csv"))
    # treated start (Experiment 19, part 19B); Experiment 11 started from a never-diseased site
    c13 = pd.read_csv(P("E19", "19B_sign_every_recall", "debridement_by_sign.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.5))
    a = axs[0]; g = a13[a13.host.isna()].sort_values("alpha")
    a.fill_between(g.alpha, np.maximum(g.u_down, 0), g.u_up, color="#eef4fc", lw=0)
    a.plot(g.alpha, g.u_up, color=BLUE, label="loss of the low state"); a.plot(g.alpha, g.u_down, "--", color=BLUE, label="loss of the high state")
    a.axhline(0.25, color=INK2, lw=0.6, ls=":")
    r = a13[a13.host.notna()]
    a.plot(r.alpha, r.relapse_ceiling, "D", color=INK, ms=4, label="relapse ceiling, treated")
    a.axvspan(0.275, g.alpha.max(), color="#f4f3ef", lw=0, hatch="////", ec=MUTED)
    a.axhline(0, color=MUTED, lw=0.5)
    a.set_xlabel("adaptation strength α"); a.set_ylabel("plaque u"); a.legend(loc="lower left", fontsize=5.6)
    panel(a, "A", "Two stable states (second set)")
    a = axs[1]
    cols = {-0.5: VIOLET, -0.25: "#a49af0", 0.25: BLUE, 0.5: "#104281"}
    for al, col in cols.items():
        for h, ls in [("eps=0.03", "-"), ("instant", "--")]:
            gg = c13[(c13.a == al) & (c13.host == h)].sort_values("visit")
            a.plot(gg.visit, gg.relapse_free_fraction, ls, color=col, lw=1.1)
    a.set_ylim(0, 1.3); a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); a.set_xlabel("recall visit with debridement"); a.set_ylabel("fraction of interval free of inflammation")
    a.legend(handles=[Line2D([], [], color=c, label=f"α = {k:g}") for k, c in cols.items()] +
             [Line2D([], [], color=INK2, label="memory 56 d"), Line2D([], [], color=INK2, ls="--", label="instant host")],
             loc="upper center", fontsize=5.4, ncol=3)
    panel(a, "B", "Repeated debridement by sign")
    fig.tight_layout(w_pad=1.5)
    save(fig, OUT, "FigureS7", {"A": a13, "B": c13})


# ============================================================== Figure S8: repeated instrumentation
def figureS8():
    # treated start (Experiment 19, part 19A); Experiment 15F started from a never-diseased site
    f2 = pd.read_csv(P("E19", "19A_reference_every_recall", "F2_debridement_every_recall.csv"))
    f15 = pd.read_csv(os.path.join(F16, "FigureS15", "FigureS15_F_data.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.5))
    a = axs[0]
    for h, col in [("eps=0.003", "#104281"), ("eps=0.01", BLUE), ("instant", ORANGE)]:
        for r, ls in [(1.0, "-"), (0.03, "--")]:
            g = f2[(f2.host == h) & (f2.r == r) & (f2.peak == 0.6)].sort_values("visit")
            a.plot(g.visit, g.free_fraction, ls, color=col)
    a.set_ylim(0, 1.05); a.set_xlabel("recall visit with debridement"); a.set_ylabel("fraction of interval free of inflammation")
    a.legend(handles=[Line2D([], [], color="#104281", label="memory 556 d"), Line2D([], [], color=BLUE, label="memory 167 d"),
                      Line2D([], [], color=ORANGE, label="instant host"), Line2D([], [], color=INK2, label="r = 1"),
                      Line2D([], [], color=INK2, ls="--", label="r = 0.03")], loc="upper left", fontsize=5.6, ncol=2)
    panel(a, "A", "Debridement at every recall, peak 0.60")
    a = axs[1]
    for h, col in [("eps=0.01", BLUE), ("instant", ORANGE)] + ([("eps=0.03", AQUA)] if "eps=0.03" in set(f15.host) else []):
        for arm, ls in [("single", "--"), ("repeated", "-")]:
            g = f15[(f15.host == h) & (f15.arm == arm)].sort_values("D0")
            if len(g):
                a.plot(g.D0, g.inflammation_free_fraction, ls, color=col, marker="o", ms=2.6)
    a.set_ylim(0, 1.05); a.set_xlabel("attachment loss at therapy D₀"); a.set_ylabel("fraction of follow-up free of inflammation")
    a.legend(handles=[Line2D([], [], color=BLUE, label="memory 167 d"), Line2D([], [], color=AQUA, label="memory 56 d"), Line2D([], [], color=ORANGE, label="instant host"),
                      Line2D([], [], color=INK2, ls="--", label="debridement once"), Line2D([], [], color=INK2, label="at every recall")],
             loc="lower left", fontsize=5.4)
    panel(a, "B", "Single vs repeated (r = 0.03, peak 0.30)")
    fig.tight_layout(w_pad=1.5)
    save(fig, OUT, "FigureS8", {"A": f2, "B": f15})


if __name__ == "__main__":
    figureS1(); s2 = figureS2(); print(s2); figureS3(); figureS4(); figureS5(); figureS6(); figureS7(); figureS8()
    print("wrote", OUT)

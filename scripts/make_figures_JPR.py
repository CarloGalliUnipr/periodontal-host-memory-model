"""Main-text Figures 1-5 of the JPR manuscript, drawn from the experiment exports.

Usage: python -I make_figures_JPR.py <data folder> <output folder>
<data folder> holds E14, E15, E16, E17, E19 (the unpacked exports of Experiments 14-17 and 19).
Figure 1C is simulated with periomemory 1.0.0; every other panel is read from the exports.
Figure 5C (repeated debridement of a treated site) is read from Experiment 19 (part 19A).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import periomemory as pm
from style import *

DATA, OUT = sys.argv[1], sys.argv[2]
E14, E15, E16, E17, E19 = [os.path.join(DATA, k) for k in ("E14", "E15", "E16", "E17", "E19")]
setup()

REF_RC = pd.read_csv(os.path.join(E14, "14A_rechallenge", "rechallenge_reference.csv"))
REF_RC = REF_RC[(REF_RC.baseline == 0.2) & (~REF_RC.no_challenge.astype(bool))
                & REF_RC.history.isin(["naive", "treated_D0"])].copy()
REF_RC["days"] = REF_RC.recovery_days
REF_RC["treated"] = REF_RC.history != "naive"
SS_RC = pd.read_csv(os.path.join(E14, "14A_rechallenge", "rechallenge_single_state.csv"))
SS_RC["days"] = SS_RC.recovery_days; SS_RC["treated"] = SS_RC.history != "naive"
SS_RC["persistent"] = SS_RC.H_after_settling > 0.5
TH_RC = pd.read_csv(os.path.join(E15, "15A_intrinsic_bistability", "A3_rechallenge.csv"))
TH_RC["treated"] = TH_RC.history != "naive"

DAYS_SHOWN = [0, 30, 90, 180, 365]


def u_at(df, level=0.2):
    """U0.2: the challenge plateau at which peak inflammation reaches `level`, by linear interpolation
    between tested plateaus. Returns (value, flag): flag '=' if interpolated, '≤' if the lowest plateau
    tested already reached the level (value is an upper bound), '>' if no plateau tested reached it (nan)."""
    g = df.sort_values("U")
    U, pk = g.U.to_numpy(float), g.peak_H.to_numpy(float)
    above = np.nonzero(pk >= level)[0]
    if len(above) == 0:
        return np.nan, ">"
    i = above[0]
    if i == 0:
        return U[0], "≤"
    return U[i - 1] + (level - pk[i - 1]) * (U[i] - U[i - 1]) / (pk[i] - pk[i - 1]), "="


# ============================================================== Figure 1
def figure1():
    fig = plt.figure(figsize=(7.0, 5.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], height_ratios=[1, 1], hspace=0.45, wspace=0.28)
    # ---------------- A: schematic
    ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(-0.03, 1.02); ax.set_ylim(0, 1); ax.axis("off")
    panel(ax, "A", "Model")
    nodes = {"u": (0.10, 0.80, "plaque\nu"), "M": (0.10, 0.36, "subgingival\ncommunity M"),
             "H": (0.50, 0.58, "host\ninflammation H"), "A": (0.88, 0.58, "response\nstate A")}
    for k, (x, y, s) in nodes.items():
        fc = "#eef4fc" if k in "HA" else "#f4f3ef"
        ec = BLUE if k == "A" else INK2
        ax.add_patch(FancyBboxPatch((x - 0.115, y - 0.09), 0.23, 0.18, boxstyle="round,pad=0.01,rounding_size=0.03",
                                    fc=fc, ec=ec, lw=0.9 if k == "A" else 0.6))
        ax.text(x, y, s, ha="center", va="center", fontsize=6.1, color=INK)

    def arrow(p, q, rad=0.0, color=INK2, ls="-", style="-|>", lw=0.9):
        ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f"arc3,rad={rad}", arrowstyle=style,
                                     mutation_scale=7, color=color, lw=lw, ls=ls))
    arrow((0.225, 0.78), (0.375, 0.65)); arrow((0.225, 0.38), (0.375, 0.52))
    arrow((0.46, 0.485), (0.225, 0.30), rad=-0.25)
    ax.text(0.36, 0.25, "G(H)", fontsize=6.2, color=INK2)
    ax.text(0.50, 0.72, "F(M + u; θ)", fontsize=6, color=INK2, ha="center", va="center")
    arrow((0.62, 0.60), (0.76, 0.60)); ax.text(0.69, 0.635, "rate ε", fontsize=5.8, color=INK2, ha="center")
    arrow((0.86, 0.485), (0.56, 0.485), rad=-0.45, color=BLUE, style="-|>")
    ax.text(0.71, 0.27, "θ = θ₀·e^(−A)", fontsize=6.2, color=BLUE, ha="center")
    ax.text(0.0, 0.10, "Debridement sets M to 0 and leaves A unchanged.\n"
            "α > 0: activity lowers the threshold (training); α < 0: raises it (tolerance).\n"
            "Instant host: A = αH at every moment (same equilibria, no memory).",
            fontsize=5.9, color=INK2, va="center")

    # ---------------- B: protocol
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "B", "Therapy and graded rechallenge")
    t = np.linspace(-30, 200, 2000)
    t_ch, dur = 110.0, 21.0
    U = 0.45
    u = np.where(t < 0, 0.2, 0.2)
    u = np.where((t >= 0) & (t < 7), 0.0, u)
    on = (t >= t_ch) & (t < t_ch + dur)
    u = np.where(on, 0.2 + (U - 0.2) * (1 - np.exp(-0.5 / UNIT_DAYS * (t - t_ch))), u)
    ax.plot(t, u, color=INK, lw=1.2)
    ax.axvspan(-30, 0, color="#f4f3ef", lw=0); ax.axvspan(t_ch, t_ch + dur, color="#eef4fc", lw=0)
    ax.axvline(0, color=INK2, lw=0.8, ls="--")
    ax.text(-15, 0.40, "estab-\nlished\ndisease", ha="center", fontsize=6, color=INK2)
    ax.text(2, 0.52, "debridement,\nthen ~7 d of\nplaque control", ha="left", va="top", fontsize=6, color=INK2)
    ax.annotate("", (t_ch - 2, 0.30), (12, 0.30), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.7))
    ax.text((12 + t_ch) / 2, 0.32, "recovery interval\n0–365 days", ha="center", fontsize=6, color=INK2)
    ax.text(t_ch + dur / 2, 0.52, "21-day\nchallenge\nto plateau U", ha="center", va="top", fontsize=6, color=INK2)
    ax.text(t_ch + dur + 4, 0.215, "baseline\nplaque", fontsize=6, color=INK2, va="bottom")
    ax.set_xlim(-30, 200); ax.set_ylim(-0.02, 0.56)
    ax.set_xticks([]); ax.set_yticks([0, 0.2, U]); ax.set_yticklabels(["0", "0.2", "U"])
    ax.set_ylabel("supragingival plaque u"); ax.set_xlabel("time (not to scale)")

    # ---------------- C: simulated time courses
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "C", "Same challenge, 90 days after therapy")
    p = pm.REFERENCE; ub, Uc = 0.2, 0.35
    rec = pm.align(pm.days_to_units(90.0)); dur_u = pm.align(pm.days_to_units(21.0)); tail = pm.days_to_units(45.0)
    src = []
    for label, eps, instant, hist, col, ls in [("naive site, host with memory", 0.01, False, "naive", GREY, "--"),
                                             ("treated site, host with memory", 0.01, False, "treated", BLUE, "-"),
                                             ("treated site, instant host", 0.01, True, "treated", ORANGE, "-")]:
        m = pm.Model(**p, eps=eps, instant=instant)
        if hist == "naive":
            st = pm.eq_state(p, ub, "low")
        else:
            st, ok, _ = pm.therapy(m, pm.eq_state(p, ub, "high"))
        u_fun = pm.schedule(ub, rec, dur_u, Uc, ka=pm.K_ACCESS)
        tr = m.simulate(st, u_fun, 0.0, float(rec + dur_u + tail), every=2)
        d = pm.units_to_days(tr["t"] - rec)
        k = d >= -10
        ax.plot(d[k], tr["H"][k], ls, color=col, label=label)
        src.append(pd.DataFrame({"series": label, "days_from_challenge_start": d[k], "H": tr["H"][k],
                                 "A": tr["A"][k], "u": tr["u"][k]}))
    ax.axvspan(0, 21, color="#eef4fc", lw=0)
    ax.text(10.5, 0.31, "challenge", ha="center", fontsize=6, color=INK2)
    ax.set_xlim(-10, 45); ax.set_ylim(0, 0.46)
    ax.set_xlabel("days from start of challenge"); ax.set_ylabel("inflammation H")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.99), fontsize=5.8)
    src_c = pd.concat(src)

    # ---------------- D: readout
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "D", "Readout: curve shift")
    g = REF_RC[REF_RC.host == "eps=0.01"]
    rows = []
    for treated, d, col, ls, lab in [(False, 90, GREY, "--", "naive site"), (True, 90, BLUE, "-", "treated site, 90 days")]:
        gg = g[(g.treated == treated) & (g.days == d)].sort_values("U")
        ax.plot(gg.U, gg.peak_H, ls, color=col, label=lab, zorder=2)
        pers = gg.persistent.astype(bool)
        ax.plot(gg.U[~pers], gg.peak_H[~pers], "o", color=col, ms=3.5, zorder=3)
        ax.plot(gg.U[pers], gg.peak_H[pers], "o", mfc="white", mec=col, ms=4, mew=0.9, zorder=3)
        v, _ = u_at(gg); rows.append((lab, v)); rows_df = gg
    ax.axhline(0.2, color=MUTED, lw=0.6, ls=":")
    (_, un), (_, ut) = rows
    ax.annotate("", (ut, 0.2), (un, 0.2), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9, mutation_scale=8))
    ax.text((un + ut) / 2, 0.235, "U₀.₂ shifted\n%.2f → %.2f" % (un, ut), ha="center", fontsize=6, color=INK)
    ax.plot([], [], "o", mfc="white", mec=INK2, label="site did not return to health")
    ax.set_xlabel("challenge plateau U"); ax.set_ylabel("peak inflammation H")
    ax.set_ylim(0, 1.05); ax.set_xlim(0.23, 0.72)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62))
    save(fig, OUT, "Figure1", {"C": src_c, "D": g[g.days == 90]})


# ============================================================== Figure 2
def curves(ax, df, host, title, letter, floor_note=True, treated_col=None):
    g = df[df.host == host]
    nv = g[~g.treated & (g.days == 0)].sort_values("U")
    ax.plot(nv.U, nv.peak_H, "--", color=GREY, lw=1.2, zorder=4)
    for i, d in enumerate(DAYS_SHOWN):
        gg = g[g.treated & (g.days == d)].sort_values("U")
        col = treated_col or SEQ[i]
        ax.plot(gg.U, gg.peak_H, "-", color=col, lw=1.1, zorder=3)
        pers = gg.persistent.astype(bool)
        ax.plot(gg.U[~pers], gg.peak_H[~pers], "o", color=col, ms=2.6, zorder=3)
        ax.plot(gg.U[pers], gg.peak_H[pers], "o", mfc="white", mec=col, ms=3.2, mew=0.7, zorder=3)
        if treated_col:
            break
    pers = nv.persistent.astype(bool)
    ax.plot(nv.U[pers], nv.peak_H[pers], "o", mfc="white", mec=GREY, ms=3.2, mew=0.7, zorder=4)
    ax.axhline(0.2, color=MUTED, lw=0.5, ls=":")
    ax.set_ylim(0, 1.05); panel(ax, letter, title)
    ax.set_xlabel("challenge plateau U"); ax.set_ylabel("peak inflammation H")


def figure2():
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 5.5))
    ax = axs.ravel()
    curves(ax[0], REF_RC, "eps=0.01", "Reference, memory 167 d", "A")
    curves(ax[1], REF_RC, "eps=0.03", "Reference, memory 56 d", "B")
    curves(ax[2], REF_RC, "instant", "Reference, instant host", "C", treated_col=ORANGE)
    ax[2].text(0.30, 0.55, "treated sites at every\ninterval = naive sites", fontsize=6, color=INK2)
    curves(ax[3], SS_RC[SS_RC.host != "eps=0.1"], "eps=0.01", "Single-state, memory 167 d", "D")
    curves(ax[4], TH_RC, "eps=0.03", "Third set, memory 56 d", "E")
    leg = [Line2D([], [], color=GREY, ls="--", label="naive")] + \
          [Line2D([], [], color=SEQ[i], label=(f"treated, {d} d" if i == 0 else f"{d} d")) for i, d in enumerate(DAYS_SHOWN)] + \
          [Line2D([], [], ls="none", marker="o", mfc="white", mec=INK2, label="open: did not return to health")]
    fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=7, fontsize=6.2)

    # ---------------- F: U0.2 relative to naive
    a = ax[5]; rows = []
    sets = [("reference", REF_RC, "-", "o", ["eps=0.01", "eps=0.03", "instant"]),
            ("third", TH_RC, "--", "s", ["eps=0.01", "eps=0.03", "instant"]),
            ("single-state", SS_RC, ":", "^", ["eps=0.01", "eps=0.03", "instant"])]
    for name, df, ls, mk, hosts in sets:
        for h in hosts:
            g = df[df.host == h]
            un, _ = u_at(g[~g.treated & (g.days == 0)])
            xs, ys, fl = [], [], []
            for d in sorted(g.days.unique()):
                if name == "single-state" and d == 0:   # day 0 is the end of conditioning
                    continue
                v, f = u_at(g[g.treated & (g.days == d)])
                xs.append(d); ys.append(v / un); fl.append(f)
                rows.append(dict(set=name, host=h, days=d, U02_treated=v, flag=f, U02_naive=un, ratio=v / un))
            xs, ys, fl = np.array(xs), np.array(ys), np.array(fl)
            if name == "single-state" and h == "instant":   # day 0 is the end of conditioning, not therapy
                pass
            a.plot(xs, ys, ls, color=HOST[h], lw=1.1)
            ok = fl == "="
            a.plot(xs[ok], ys[ok], mk, color=HOST[h], ms=3)
            a.plot(xs[~ok], ys[~ok], "v", mfc="white", mec=HOST[h], ms=3.6, mew=0.8)
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.set_xticks([0, 90, 180, 270, 365]); a.set_ylim(0.4, 1.06)
    a.set_xlabel("days after therapy"); a.set_ylabel("U₀.₂ treated / naive")
    panel(a, "F", "Shift and its recovery")
    hl = [Line2D([], [], color=BLUE, label="F: memory 167 d"), Line2D([], [], color=AQUA, label="memory 56 d"),
          Line2D([], [], color=ORANGE, label="instant host"),
          Line2D([], [], color=INK2, ls="-", marker="o", ms=3, label="reference"),
          Line2D([], [], color=INK2, ls="--", marker="s", ms=3, label="third"),
          Line2D([], [], color=INK2, ls=":", marker="^", ms=3, label="single-state"),
          Line2D([], [], ls="none", marker="v", mfc="white", mec=INK2, label="F: at or below lowest U tested")]
    fig.legend(handles=hl, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=7, fontsize=6.2, columnspacing=1.0)
    fig.tight_layout(w_pad=1.0, h_pad=1.3, rect=(0, 0.09, 1, 1))
    u02 = pd.DataFrame(rows)
    save(fig, OUT, "Figure2", {"A_C": REF_RC, "D": SS_RC, "E": TH_RC, "F": u02})
    return u02


# ============================================================== Figure 3
def figure3():
    F = os.path.join(E16, "figures", "Figure4")
    a4 = pd.read_csv(os.path.join(F, "Figure4_A_data.csv"))
    b4 = pd.read_csv(os.path.join(F, "Figure4_B_data.csv"))
    c4 = pd.read_csv(os.path.join(F, "Figure4_C_data.csv"))
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.45))
    a = axs[0]
    a.axvspan(-0.55, 0, color="#f1effc", lw=0); a.axvspan(0, 0.55, color="#eef4fc", lw=0)
    s = a4[a4.set == "second"].sort_values("alpha")
    ok = s.ratio_treated.notna() & (s.alpha > -0.5)
    a.plot(s.alpha[ok], s.ratio_treated[ok], "-o", color=INK, ms=3.2, label="second set (memory 56 d)")
    nm = s[s.alpha == -0.5]   # conversion not monotone in challenge size: ratio of lowest converting challenges
    a.plot(nm.alpha, [0.470 / 0.425] * len(nm), "o", mfc="white", mec=INK, ms=4, mew=0.9)
    a.plot(s.alpha[~ok], [0.62] * (~ok).sum(), "x", color=INK, ms=5, mew=1.1)
    a.text(0.47, 0.645, "therapy\nfails", fontsize=5.8, color=INK2, ha="center")
    for eps, col in [(0.03, AQUA), (0.01, BLUE)]:
        t = a4[(a4.set == "third") & (a4.eps == eps) & a4.ratio.notna()].sort_values("alpha")
        a.plot(t.alpha, t.ratio, "--s", color=col, ms=3, label=f"third set (memory {round(1/eps*UNIT_DAYS)} d)")
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.text(-0.27, 1.125, "tolerance", ha="center", fontsize=6.5, color=VIOLET)
    a.text(0.27, 1.125, "training", ha="center", fontsize=6.5, color=BLUE)
    a.set_xlim(-0.55, 0.55); a.set_ylim(0.6, 1.15)
    a.set_xlabel("adaptation strength α"); a.set_ylabel("critical challenge, treated / naive")
    a.legend(loc="lower left", fontsize=5.6); panel(a, "A", "Sign of the shift")

    a = axs[1]
    cols = {-0.25: VIOLET, 0.0: GREY, 0.25: BLUE}
    for al, col in cols.items():
        g = b4[(b4.alpha == al)].sort_values("eps_T_rec")
        a.plot(g.eps_T_rec, g.ratio_treated, "-o", color=col, ms=3, label=f"α = {al:g}")
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.set_xlabel("memory times since therapy"); a.set_ylabel("critical challenge, treated / naive")
    a.legend(loc="lower right", fontsize=5.8); panel(a, "B", "Recovery (second set)")

    a = axs[2]
    for hist, col, lab in [("gingivitis", AQUA, "after resolved gingivitis"), ("treated_protocol", INK, "after therapy")]:
        g = c4[c4.history == hist].sort_values("eps_T_rec")
        a.plot(g.eps_T_rec, g.ratio, "-o", color=col, ms=3, label=lab)
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.set_ylim(0.95, 1.07)
    a.set_xlabel("memory times since conditioning"); a.set_ylabel("critical challenge, conditioned / naive")
    a.legend(loc="upper right", fontsize=5.8); panel(a, "C", "Dose-dependent host")
    fig.tight_layout(w_pad=1.2)
    save(fig, OUT, "Figure3", {"A": a4, "B": b4, "C": c4})


# ============================================================== Figure 4
def figure4():
    A17 = os.path.join(E17, "17A_resolution_rate"); C17 = os.path.join(E17, "17C_community_memory")
    hl = pd.read_csv(os.path.join(A17, "A2_half_life.csv")); th = pd.read_csv(os.path.join(A17, "A2_therapy.csv"))
    ih = pd.read_csv(os.path.join(A17, "A2_instant_host.csv")); rc = pd.read_csv(os.path.join(A17, "A2_relapse_ceiling.csv"))
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.9)); ax = axs.ravel()
    xt = [0.25, 0.5, 1, 2, 4]

    a = ax[0]
    for h in ["eps=0.01", "eps=0.03"]:
        g = hl[(hl.host == h) & (hl.r == 1.0)].sort_values("rho")
        a.plot(g.rho, g.ratio_immediately, "-o", color=HOST[h], ms=3)
        fail = th[(th.host == h) & (th.r == 1.0) & (~th.restored.astype(bool))]
        a.plot(fail.rho, np.full(len(fail), 0.43 if h == "eps=0.01" else 0.45), "x", color=HOST[h], ms=5, mew=1.1)
    g = ih[ih.r == 1.0].sort_values("rho")
    a.plot(g.rho, g.crit_instant_treated / g.crit_instant_naive, "-o", color=ORANGE, ms=3)
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.text(0.27, 0.40, "×  therapy fails", fontsize=5.8, color=INK2)
    a.set_xscale("log"); a.set_xticks(xt); a.set_xticklabels([f"{v:g}" for v in xt]); a.minorticks_off()
    a.set_ylim(0.38, 1.05); a.set_xlabel("resolution rate ρ"); a.set_ylabel("critical challenge, treated / naive")
    a.legend(handles=[Line2D([], [], color=BLUE, label="memory 167 d"), Line2D([], [], color=AQUA, label="memory 56 d"),
                      Line2D([], [], color=ORANGE, label="instant host")], loc="center right", fontsize=5.8)
    panel(a, "A", "Slow resolution")

    a = ax[1]
    for h in ["eps=0.01", "eps=0.03", "instant"]:
        g = rc[(rc.host == h) & (rc.r == 1.0)].sort_values("rho")
        failed = set(th[(th.host == h) & (th.r == 1.0) & (~th.restored.astype(bool))].rho)
        g = g[~g.rho.isin(failed)]
        y = g.ceiling_over_fold.values; ok = ~np.isnan(y)
        a.plot(g.rho[ok], y[ok], "-o", color=HOST[h], ms=3)
        a.plot(g.rho[~ok], np.ones((~ok).sum()), "o", mfc="white", mec=HOST[h], ms=3.6, mew=0.8)
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.text(0.26, 1.03, "open: no relapse below the static threshold", fontsize=5.4, color=INK2)
    a.text(0.26, 0.43, "memory hosts at ρ = 0.25: therapy fails", fontsize=5.4, color=INK2)
    a.set_xscale("log"); a.set_xticks(xt); a.set_xticklabels([f"{v:g}" for v in xt]); a.minorticks_off()
    a.set_ylim(0.4, 1.08); a.set_xlabel("resolution rate ρ"); a.set_ylabel("relapse ceiling / static threshold")
    panel(a, "B", "Relapse without challenge")

    a = ax[2]
    s12 = pd.read_csv(os.path.join(E16, "figures", "FigureS12", "FigureS12_F_data.csv"))
    g = s12[(s12.eps == 0.01) & (s12.U == 0.4)].sort_values("days")
    a.plot(g.days, g.peak_treated, "-o", color=BLUE, ms=3, label="treated")
    p = g[g.persistent_treated.astype(bool)]
    a.plot(p.days, p.peak_treated, "o", mfc="white", mec=BLUE, ms=4, mew=0.9)
    a.plot(g.days, g.peak_naive_matched_same_increment, "-s", color=INK, ms=3, label="naive, matched resting H")
    a.plot(g.days, g.peak_naive_standard, "--", color=GREY, label="naive")
    a.set_xticks([0, 90, 180, 270, 365]); a.set_ylim(0, 1.05)
    a.set_xlabel("days after therapy"); a.set_ylabel("peak inflammation H")
    a.legend(loc="center right", fontsize=5.6); panel(a, "C", "Matched resting bleeding")

    w = pd.read_csv(os.path.join(C17, "C2_window_and_matching.csv"))
    a = ax[3]
    g = w[(w.challenge == "2/0.01") & (w.beta == 0.5) & (w.r == 1.0) & (w.host == "instant")]
    comm = {"mu=0.01": "#1c5cab", "mu=0.03": "#3987e5", "mu=0.1": "#9ec5f4"}
    for c, col in comm.items():
        gg = g[g.community == c].sort_values("days")
        a.plot(gg.days, gg.peak_treated, "-o", color=col, ms=3,
               label=f"treated, community {round(1/float(c[3:])*UNIT_DAYS)} d")
    gg = g[g.community == "mu=0.01"].sort_values("days")
    p = gg[gg.persistent_treated.astype(bool)]
    a.plot(p.days, p.peak_treated, "o", mfc="white", mec=comm["mu=0.01"], ms=4.2, mew=0.9)
    a.plot(gg.days, gg.peak_matched, "s", mfc="none", mec=INK, ms=4.2, mew=0.7, label="naive, same community state")
    a.plot(gg.days, gg.peak_naive, "--", color=GREY, label="naive")
    a.set_xticks([0, 90, 180, 270, 365]); a.set_ylim(0, 1.05)
    a.set_xlabel("days after therapy"); a.set_ylabel("peak inflammation H")
    a.legend(loc="upper right", fontsize=5.4); panel(a, "D", "Community memory only")

    a = ax[4]
    sel = w[(w.challenge == "2/0.01") & (w.beta == 0.5) & (w.r == 0.03)]
    for (h, c), col, mk, lab in [(("eps=0.01", "instant"), BLUE, "o", "host only"),
                                 (("instant", "mu=0.01"), "#9ec5f4", "s", "community only"),
                                 (("eps=0.01", "mu=0.01"), INK, "o", "both")]:
        gg = sel[(sel.host == h) & (sel.community == c)].sort_values("days")
        a.plot(gg.days, gg.ratio_naive, "-" + mk, color=col, ms=3, label=lab)
        if c != "instant":
            a.plot(gg.days, gg.ratio_matched, ":", color=col, lw=1.1)
    a.axhline(1, color=MUTED, lw=0.5, ls=":")
    a.set_xticks([0, 90, 180, 270, 365]); a.set_ylim(0.35, 1.05)
    a.set_xlabel("days after therapy"); a.set_ylabel("critical challenge, treated / naive")
    a.legend(handles=a.get_legend_handles_labels()[0] + [Line2D([], [], color=INK2, ls=":", label="vs naive with same community")],
             loc="lower right", fontsize=5.4)
    panel(a, "E", "Critical challenge")

    a = ax[5]
    dc = pd.read_csv(os.path.join(C17, "C2_decomposition.csv"))
    sel2 = dc[(dc.challenge == "2/0.01") & (dc.host == "eps=0.01") & (dc.community == "mu=0.01") & (dc.r == 0.03)
              & dc.days.isin([0.0, 90.0])].sort_values(["beta", "days"])
    x = np.arange(len(sel2)); wb = 0.22
    a.bar(x - wb, sel2.deficit_host, wb, color=BLUE, label="host memory")
    a.bar(x, sel2.deficit_community, wb, color="#9ec5f4", label="community memory")
    a.bar(x + wb, sel2.deficit_both, wb, color=INK, label="both")
    a.plot(x + wb, sel2.deficit_host + sel2.deficit_community, "_", color=ORANGE, ms=10, mew=1.6, label="sum of the two")
    a.set_xticks(x); a.set_xticklabels([f"β {b:g}\nday {int(d)}" for b, d in zip(sel2.beta, sel2.days)])
    a.set_ylabel("deficit, 1 − treated/naive"); a.set_ylim(0, 0.85)
    a.legend(loc="upper center", fontsize=5.4, handlelength=1.2, ncol=2, columnspacing=0.8); panel(a, "F", "Host × community")
    fig.tight_layout(w_pad=1.0, h_pad=1.3)
    save(fig, OUT, "Figure4", {"A": hl.merge(th, on=["r", "host", "rho"], how="outer"), "A_instant": ih, "B": rc,
                               "C": g, "D": w[(w.challenge == "2/0.01") & (w.beta == 0.5) & (w.r == 1.0) & (w.host == "instant")],
                               "E": sel, "F": sel2})


# ============================================================== Figure 5
def figure5():
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.8)); ax = axs.ravel()
    a = ax[0]
    b4 = pd.read_csv(os.path.join(E14, "14B_tipping_window", "B4_relapse_ceiling_summary.csv"))
    g = b4[(b4.kind == "treated_protocol") & (b4.host != "instant")].copy()
    g["eps"] = g.host.str[4:].astype(float); g["memory_days"] = 1 / g.eps * UNIT_DAYS
    g = g.sort_values("memory_days")
    ok = g.ceiling.notna()
    a.plot(g.memory_days[ok], g.ceiling[ok] / FOLD_REF, "-o", color=BLUE, ms=3.5, label="host with memory")
    a.plot(g.memory_days[~ok], np.ones((~ok).sum()), "o", mfc="white", mec=BLUE, ms=4, mew=0.9, label="no relapse below static threshold")
    a.axhline(1, color=ORANGE, lw=1.2, label="instant host (static threshold)")
    a.set_xscale("log"); a.set_xticks([3, 10, 30, 100, 300, 1000]); a.set_xticklabels(["3", "10", "30", "100", "300", "1000"])
    a.minorticks_off(); a.set_ylim(0.4, 1.08)
    a.set_xlabel("memory time (days)"); a.set_ylabel("relapse ceiling / static threshold")
    a.legend(loc="lower left", fontsize=5.8); panel(a, "A", "Relapse without challenge")

    a = ax[1]
    s11 = pd.read_csv(os.path.join(E16, "figures", "FigureS11", "FigureS11_G_H_data.csv"))
    for dur, col in [(7, "#9ec5f4"), (14, "#3987e5"), (21, "#104281")]:
        gg = s11[(s11.host == "eps=0.01") & (s11.duration_days == dur)].sort_values("days")
        cens = gg.lowest_converting_U.isna()
        a.plot(gg.days, gg.U_safe, "-", color=col, label=f"{dur}-day challenge")
        a.plot(gg.days[~cens], gg.U_safe[~cens], "o", color=col, ms=3)
        a.plot(gg.days[cens], gg.U_safe[cens], "^", mfc="white", mec=col, ms=3.8, mew=0.8)
    a.set_xticks([0, 90, 180, 270, 365]); a.set_ylim(0.2, 0.75)
    a.set_xlabel("days after therapy"); a.set_ylabel("largest plateau without conversion")
    a.legend(handles=a.get_legend_handles_labels()[0] + [Line2D([], [], ls="none", marker="^", mfc="white", mec=INK2,
             label="no conversion up to 0.70 (naive sites: none)")], loc="lower right", fontsize=5.8)
    panel(a, "B", "Safe challenge size (memory 167 d)")

    a = ax[2]
    # treated start (Experiment 19, part 19A); Experiment 15F started from a never-diseased site
    f2 = pd.read_csv(os.path.join(E19, "19A_reference_every_recall", "F2_debridement_every_recall.csv"))
    for (h, col) in [("eps=0.01", BLUE), ("eps=0.003", "#104281"), ("instant", ORANGE)]:
        for r, ls in [(1.0, "-"), (0.03, "--")]:
            gg = f2[(f2.host == h) & (f2.r == r) & (f2.peak == 0.75)].sort_values("visit")
            a.plot(gg.visit, gg.free_fraction, ls, color=col)
    a.set_ylim(0, 0.6); a.set_xlabel("recall visit with debridement"); a.set_ylabel("fraction of interval free of inflammation")
    a.legend(handles=[Line2D([], [], color="#104281", label="memory 556 d"), Line2D([], [], color=BLUE, label="memory 167 d"),
                      Line2D([], [], color=ORANGE, label="instant host"), Line2D([], [], color=INK2, label="r = 1"),
                      Line2D([], [], color=INK2, ls="--", label="r = 0.03")], loc="upper right", fontsize=5.6, ncol=2)
    panel(a, "C", "Repeated debridement of a treated site")

    a = ax[3]
    b2 = pd.read_csv(os.path.join(E17, "17B_adjunct", "B2_challenge_with_adjunct.csv"))
    g0 = b2[b2.eps == 0.01]
    for hist, arm, col, ls, lab in [("treated", "none", INK, "-", "treated, no course"),
                                    ("treated", "concurrent lam=1.0", BLUE, "-", "treated, course during challenge"),
                                    ("treated", "early lam=1.0", AQUA, "--", "treated, course on days 0–28"),
                                    ("naive", "none", GREY, "--", "naive, no course"),
                                    ("naive", "concurrent lam=1.0", GREY, ":", "naive, course during challenge")]:
        gg = g0[(g0.history == hist) & (g0.arm == arm)].sort_values("days")
        a.plot(gg.days, gg.peak_H, ls, color=col, marker="o" if hist == "treated" else None, ms=3, label=lab)
        p = gg[gg.persistent.astype(bool)]
        a.plot(p.days, p.peak_H, "o", mfc="white", mec=col, ms=4.2, mew=0.9)
    a.set_xticks([0, 90, 180, 270, 365]); a.set_ylim(0, 1.05)
    a.set_xlabel("days after therapy"); a.set_ylabel("peak inflammation H")
    a.legend(loc="upper left", bbox_to_anchor=(0.06, 0.9), fontsize=5.4); panel(a, "D", "Pro-resolving course (memory 167 d)")
    fig.tight_layout(w_pad=1.2, h_pad=1.3)
    save(fig, OUT, "Figure5", {"A": g, "B": s11, "C": f2, "D": g0})


if __name__ == "__main__":
    figure1(); u02 = figure2(); figure3(); figure4(); figure5()
    print(u02.pivot_table(index=["set", "host"], columns="days", values="U02_treated").round(3))
    print("wrote", OUT)

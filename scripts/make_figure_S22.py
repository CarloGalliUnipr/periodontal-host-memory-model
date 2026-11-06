"""Figure S9 of the JPR Supporting Information (formerly Figure S22), from the
Experiment_17_Resolution_Community_Memory export.

Usage: python make_figure_S22.py <export folder> <output folder>
Writes FigureS22.pdf, FigureS22.png (600 dpi) and the source data of every panel.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D

EXP, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
A = os.path.join(EXP, "17A_resolution_rate"); B = os.path.join(EXP, "17B_adjunct"); C = os.path.join(EXP, "17C_community_memory")

BLUE, ORANGE, AQUA, GREY, INK = "#2a78d6", "#eb6834", "#1baf7a", "#8c8b86", "#0b0b0b"
SEQ = {"mu=0.01": "#1f5fae", "mu=0.03": "#4f93e3", "mu=0.1": "#9cc3f0"}
HOSTC = {"eps=0.01": BLUE, "eps=0.03": AQUA, "instant": ORANGE}
FOLD = 0.4682392215

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.titlesize": 8, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "legend.frameon": False, "legend.fontsize": 6.5,
                     "lines.linewidth": 1.4, "lines.markersize": 4})
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.9))
ax = axs.ravel()
src = {}

# ---------------------------------------------------------------- A: window against resolution rate
h = pd.read_csv(os.path.join(A, "A2_half_life.csv"))
th = pd.read_csv(os.path.join(A, "A2_therapy.csv"))
for host in ["eps=0.01", "eps=0.03"]:
    for r, ls in [(1.0, "-"), (0.03, "--")]:
        g = h[(h.host == host) & (h.r == r)].sort_values("rho")
        ax[0].plot(g.rho, g.ratio_immediately, ls, color=HOSTC[host], marker="o", ms=3)
        fail = th[(th.host == host) & (th.r == r) & (~th.restored)]
        ax[0].plot(fail.rho, np.full(len(fail), 0.42), "x", color=HOSTC[host], ms=5, mew=1.2)
src["A"] = h
ax[0].set_xscale("log"); ax[0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); ax[0].set_xticks([0.25, 0.5, 1, 2, 4]); ax[0].set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax[0].set_ylim(0.40, 0.60); ax[0].set_xlabel("resolution rate ρ"); ax[0].set_ylabel("treated / naive critical challenge")
ax[0].set_title("A  Window after therapy")
ax[0].text(0.27, 0.428, "therapy fails (r = 1)", fontsize=6, color=INK)

# ---------------------------------------------------------------- B: relapse ceiling against resolution rate
rc = pd.read_csv(os.path.join(A, "A2_relapse_ceiling.csv"))
for host in ["eps=0.01", "eps=0.03", "instant"]:
    for r, ls in [(1.0, "-"), (0.03, "--")]:
        g = rc[(rc.host == host) & (rc.r == r)].sort_values("rho")
        y = g.ceiling_over_fold.values
        ok = ~np.isnan(y)
        ax[1].plot(g.rho[ok], y[ok], ls, color=HOSTC[host], marker="o", ms=3)
        ax[1].plot(g.rho[~ok], np.ones((~ok).sum()), "o", mfc="white", mec=HOSTC[host], ms=3.5)
src["B"] = rc
ax[1].axhline(1, color=GREY, lw=0.6, ls=":")
ax[1].set_xscale("log"); ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); ax[1].set_xticks([0.25, 0.5, 1, 2, 4]); ax[1].set_xticklabels(["0.25", "0.5", "1", "2", "4"])
ax[1].set_ylim(0.25, 1.05); ax[1].set_xlabel("resolution rate ρ"); ax[1].set_ylabel("relapse ceiling / static threshold")
ax[1].set_title("B  Relapse without challenge")

# ---------------------------------------------------------------- C: adjunct, relapse ceiling against course start
cc = pd.read_csv(os.path.join(B, "B1_relapse_ceiling.csv"))
for eps, col in [(0.03, AQUA), (0.01, BLUE)]:
    base = float(cc[(cc.r == 1.0) & (cc.eps == eps) & (cc.arm == "none")].relapse_ceiling.iloc[0])
    ax[2].axhline(base, color=col, lw=0.8, ls=":")
    for lam, ls in [(1.0, "-"), (0.5, "--")]:
        g = cc[(cc.r == 1.0) & (cc.eps == eps) & (cc.lam == lam)].sort_values("start_days")
        ax[2].plot(g.start_days, g.relapse_ceiling, ls, color=col, marker="o", ms=3)
src["C"] = cc
ax[2].text(92, 0.358, "ε = 0.03", color=AQUA, fontsize=6, va="bottom", ha="right")
ax[2].text(92, 0.297, "ε = 0.01", color=BLUE, fontsize=6, va="bottom", ha="right")
ax[2].set_xticks([0, 30, 60, 90]); ax[2].set_xlabel("course start (days after therapy)")
ax[2].set_ylabel("relapse ceiling (baseline plaque)"); ax[2].set_title("C  Adjunct and relapse")

# ---------------------------------------------------------------- D: adjunct, challenge in the window
b2 = pd.read_csv(os.path.join(B, "B2_challenge_with_adjunct.csv"))
g0 = b2[b2.eps == 0.01]
for hist, arm, col, ls, mk in [("treated", "none", INK, "-", "o"), ("treated", "concurrent lam=1.0", BLUE, "-", "o"),
                               ("treated", "early lam=1.0", AQUA, "--", "s"), ("naive", "none", GREY, "--", None),
                               ("naive", "concurrent lam=1.0", GREY, ":", None)]:
    g = g0[(g0.history == hist) & (g0.arm == arm)].sort_values("days")
    ax[3].plot(g.days, g.peak_H, ls, color=col, marker=mk, ms=3)
    p = g[g.persistent]
    ax[3].plot(p.days, p.peak_H, "o", mfc="white", mec=col, ms=5.5, mew=0.9)
src["D"] = g0
ax[3].set_xticks([0, 90, 180, 270, 365]); ax[3].set_ylim(0, 1.7); ax[3].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax[3].set_xlabel("days after therapy"); ax[3].set_ylabel("peak inflammation H"); ax[3].set_title("D  Adjunct and challenge")

# ---------------------------------------------------------------- E: community memory alone, calibrated challenge
w = pd.read_csv(os.path.join(C, "C2_window_and_matching.csv"))
g = w[(w.challenge == "2/0.01") & (w.beta == 0.5) & (w.r == 1.0) & (w.host == "instant")]
for com, col in [("mu=0.01", SEQ["mu=0.01"]), ("mu=0.03", SEQ["mu=0.03"]), ("mu=0.1", SEQ["mu=0.1"])]:
    gg = g[g.community == com].sort_values("days")
    ax[4].plot(gg.days, gg.peak_treated, "-", color=col, marker="o", ms=3)
    if com == "mu=0.01":
        ax[4].plot(gg.days, gg.peak_matched, "s", mfc="none", mec=INK, ms=4.5, mew=0.7)
        p = gg[gg.persistent_treated]
        ax[4].plot(p.days, p.peak_treated, "o", mfc="white", mec=col, ms=5.5, mew=0.9)
gg = g[g.community == "mu=0.01"].sort_values("days")
ax[4].plot(gg.days, gg.peak_naive, "--", color=GREY)
src["E"] = g
ax[4].set_xticks([0, 90, 180, 270, 365]); ax[4].set_ylim(0, 1.05)
ax[4].set_xlabel("days after therapy"); ax[4].set_ylabel("peak inflammation H")
ax[4].set_title("E  Community memory alone")

# ---------------------------------------------------------------- F: decomposition
dc = pd.read_csv(os.path.join(C, "C2_decomposition.csv"))
sel = dc[(dc.challenge == "2/0.01") & (dc.host == "eps=0.01") & (dc.community == "mu=0.01") & (dc.r == 0.03)
         & dc.days.isin([0.0, 90.0])].sort_values(["beta", "days"])
x = np.arange(len(sel)); wbar = 0.22
ax[5].bar(x - wbar, sel.deficit_host, wbar, color=BLUE, label="host memory")
ax[5].bar(x, sel.deficit_community, wbar, color=SEQ["mu=0.1"], label="community memory")
ax[5].bar(x + wbar, sel.deficit_both, wbar, color=INK, label="both")
ax[5].plot(x + wbar, sel.deficit_host + sel.deficit_community, "_", color=ORANGE, ms=11, mew=1.6, label="sum of the two")
ax[5].set_xticks(x); ax[5].set_xticklabels([f"β={b}\nday {int(d)}" for b, d in zip(sel.beta, sel.days)])
ax[5].set_ylabel("deficit, 1 − treated/naive"); ax[5].set_ylim(0, 0.85)
ax[5].set_title("F  Host × community")
ax[5].legend(loc="upper center", ncol=2, fontsize=5.8, handlelength=1.2, columnspacing=0.8)
src["F"] = sel

# ---------------------------------------------------------------- shared legends
hosts = [Line2D([], [], color=BLUE, label="ε = 0.01"), Line2D([], [], color=AQUA, label="ε = 0.03"),
         Line2D([], [], color=ORANGE, label="instant host"), Line2D([], [], color=INK, ls="-", label="r = 1"),
         Line2D([], [], color=INK, ls="--", label="r = 0.03")]
ax[1].legend(handles=hosts, loc="lower right", fontsize=6, handlelength=1.8)
ax[2].legend(handles=[Line2D([], [], color=INK, ls="-", label="λ = 1"), Line2D([], [], color=INK, ls="--", label="λ = 0.5"),
                      Line2D([], [], color=INK, ls=":", lw=0.8, label="no course")], loc="upper left", bbox_to_anchor=(0.0, 0.66), fontsize=6, handlelength=1.6)
ax[3].legend(handles=[Line2D([], [], color=INK, label="treated"),
                      Line2D([], [], color=BLUE, label="treated + concurrent course"),
                      Line2D([], [], color=AQUA, ls="--", label="treated + course on days 0–28"),
                      Line2D([], [], color=GREY, ls="--", label="naive"),
                      Line2D([], [], color=GREY, ls=":", label="naive + concurrent course")],
             loc="upper left", bbox_to_anchor=(0.0, 1.02), fontsize=5.6, handlelength=1.5)
ax[4].legend(handles=[Line2D([], [], color=SEQ[k], label=f"treated, μ = {k[3:]}") for k in ["mu=0.01", "mu=0.03", "mu=0.1"]] +
             [Line2D([], [], ls="none", marker="s", mfc="none", mec=INK, label="naive, same B (μ = 0.01)"),
              Line2D([], [], color=GREY, ls="--", label="naive")], loc="upper right", fontsize=5.8, handlelength=1.8)

fig.tight_layout(w_pad=1.2, h_pad=1.4)
fig.savefig(os.path.join(OUT, "FigureS22.pdf"))
fig.savefig(os.path.join(OUT, "FigureS22.png"), dpi=600)
for k, v in src.items():
    v.to_csv(os.path.join(OUT, f"FigureS22_{k}_source.csv"), index=False)
print("wrote", OUT)

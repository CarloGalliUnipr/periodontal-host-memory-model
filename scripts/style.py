"""Shared plotting style for the JPR figures."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from periomemory import UNIT_DAYS

BLUE, ORANGE, AQUA, GREY, INK = "#2a78d6", "#eb6834", "#1baf7a", "#8c8b86", "#0b0b0b"
INK2, MUTED, VIOLET = "#4a4945", "#b5b3ad", "#6b5bd2"
SEQ = ["#104281", "#1c5cab", "#3987e5", "#6da7ec", "#9ec5f4"]
HOST = {"eps=0.01": BLUE, "eps=0.03": AQUA, "instant": ORANGE}
FOLD_REF = 0.4682392215

__all__ = ["BLUE", "ORANGE", "AQUA", "GREY", "INK", "INK2", "MUTED", "VIOLET", "SEQ", "HOST",
           "FOLD_REF", "UNIT_DAYS", "setup", "panel", "save"]


def setup():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.6,
                         "ytick.major.width": 0.6, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
                         "axes.labelsize": 7.2, "axes.titlesize": 7.8, "axes.titleweight": "bold",
                         "axes.titlelocation": "left", "legend.frameon": False, "legend.fontsize": 6.2,
                         "lines.linewidth": 1.3, "lines.markersize": 4, "savefig.bbox": "tight"})


def panel(ax, letter, title):
    ax.set_title(f"{letter}  {title}")


def save(fig, out, name, src=None):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, f"{name}.pdf"))
    fig.savefig(os.path.join(out, f"{name}.png"), dpi=600)
    fig.savefig(os.path.join(out, f"{name}.tif"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    for k, v in (src or {}).items():
        v.to_csv(os.path.join(out, f"{name}_{k}_source.csv"), index=False)
    plt.close(fig)

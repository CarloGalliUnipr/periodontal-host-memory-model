"""periomemory: host-response memory model of a periodontal site.

Reference implementation accompanying Galli C., "Host-response memory under periodic plaque
removal: a dynamical model of post-therapy vulnerability and relapse in periodontitis".
"""

__version__ = "1.0.0"

from .parameters import (REFERENCE, SECOND, SINGLE_STATE, DT, K_REGROW, K_ACCESS, GAMMA, H_D,
                         T_PC, UNIT_DAYS, days_to_units, units_to_days)
from .model import Model, hill, const, schedule, regrowth, align
from .equilibria import u_curve, folds, eq_H, eq_state, stage_boundaries, tile
from .protocols import (therapy, init_history, bisect, critical_challenge, half_life, ratchet,
                        ratchet_threshold, theta_map, episodic)

__all__ = [
    "REFERENCE", "SECOND", "SINGLE_STATE", "DT", "K_REGROW", "K_ACCESS", "GAMMA", "H_D", "T_PC",
    "UNIT_DAYS", "days_to_units", "units_to_days",
    "Model", "hill", "const", "schedule", "regrowth", "align",
    "u_curve", "folds", "eq_H", "eq_state", "stage_boundaries", "tile",
    "therapy", "init_history", "bisect", "critical_challenge", "half_life", "ratchet",
    "ratchet_threshold", "theta_map", "episodic",
]

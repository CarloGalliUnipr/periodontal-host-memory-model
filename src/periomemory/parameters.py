"""Reference parameter sets and protocol constants.

All quantities are dimensionless. The time mapping (one model unit ~ 1.7 days,
from anchoring a recall interval T = 60 to about 100 days) is illustrative and is
used only to report results in days.
"""

# ----------------------------------------------------------------- parameter sets
#: Bistable reference set used throughout the paper (Table S1).
REFERENCE = dict(n=4.0, q=4.0, K=0.9, th0=0.8, r=1.0, a=1.0)

#: Second set: bistable even without adaptation (used for the sign of adaptation).
SECOND = dict(n=4.0, q=4.0, K=0.5, th0=0.6, r=1.0, a=1.0)

#: Single-state (monostable) set: no result can depend on a second stable state.
SINGLE_STATE = dict(n=2.0, q=1.0, K=0.9, th0=1.0, r=1.0, a=1.0)

# ----------------------------------------------------------------- numerics
#: Fixed RK4 step used in every experiment.
DT = 0.1

# ----------------------------------------------------------------- forcing and protocols
#: Plaque regrowth rate, u(tau) = u_max (1 - exp(-k tau)).
K_REGROW = 0.05
#: Fast regrowth of accessible plaque (experimental-gingivitis challenge, home care).
K_ACCESS = 0.5
#: Support of the subgingival community by attachment loss (gamma D).
GAMMA = 0.8
#: Inflammation above which attachment loss accrues.
H_D = 0.8
#: Perfect supragingival plaque control after debridement (~7 days).
T_PC = 4.0

# ----------------------------------------------------------------- time mapping
#: Days per model unit (T = 60 units ~ 100 days).
UNIT_DAYS = 100.0 / 60.0


def days_to_units(d):
    """Convert days to model time units (illustrative mapping)."""
    return d / UNIT_DAYS


def units_to_days(t):
    """Convert model time units to days (illustrative mapping)."""
    return t * UNIT_DAYS

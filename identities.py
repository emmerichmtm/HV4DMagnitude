"""The three exact identities the report's four-prefix normal form rests on.

Each function here has a brute-force twin in ``tests/`` so that the identity is
checked rather than assumed.  Where the report's statement and the geometry
disagree about an endpoint, the geometry wins and the difference is documented.

Deviation from the report, recorded deliberately
------------------------------------------------
The report writes the easy residual as ``[x_j <= f_ij(x_i)]`` and states the
first branch of the staircase formula under ``B <= f(A)``.  The region actually
left uncovered by a union of grounded quadrants is ``{x_j < f_ij(x_i)}``, and
with that strict form the first branch holds under ``B < f(A)``.  Under
Lebesgue measure the two conventions agree.  Under the magnitude measure they
differ exactly when ``f`` attains the value 0 (the atom at ``x_j = 0`` then has
mass), so we implement the strict form and expose ``closed=`` to reproduce the
report's variant for comparison.
"""

from __future__ import annotations

import math
from typing import Callable, List, Sequence, Tuple

from core import INF, Measure1D, Step

__all__ = [
    "prefix_bound_ge", "prefix_bound_gt", "prefix_mass", "cum_closed", "cum_open",
    "source_prefix", "threshold_antiderivative", "meet_contract",
]


# --------------------------------------------------------------------------- #
# Prefix bounds with explicit open/closed semantics
# --------------------------------------------------------------------------- #


def prefix_bound_ge(f: Step, u: float) -> float:
    """Exclusive ``theta`` with ``{x : f(x) >= u} == [0, theta)`` (f nonincreasing)."""
    if f.vs[0] < u:
        return 0.0
    last = 0
    for j, v in enumerate(f.vs):
        if v >= u:
            last = j
        else:
            break
    return INF if last == len(f.vs) - 1 else f.xs[last]


def prefix_bound_gt(f: Step, u: float) -> float:
    """Exclusive ``theta`` with ``{x : f(x) > u} == [0, theta)`` (f nonincreasing)."""
    if not (f.vs[0] > u):
        return 0.0
    last = 0
    for j, v in enumerate(f.vs):
        if v > u:
            last = j
        else:
            break
    return INF if last == len(f.vs) - 1 else f.xs[last]


def prefix_mass(mu: Measure1D, theta: float, closed: bool) -> float:
    """Mass of ``[0, theta]`` or ``[0, theta)``.

    These differ only at ``theta == 0``, which is exactly where the magnitude
    atom lives -- hence the explicit flag instead of a single helper.
    """
    return mu.prefix_closed(theta) if closed else mu.prefix_open(theta)


def cum_closed(weight: Callable[[float], float], grid: Sequence[float],
               mu: Measure1D, t: float) -> float:
    """``int_{[0,t]} weight d mu`` for a grid-aligned piecewise-constant weight."""
    return _cum(weight, grid, mu, t, closed=True)


def cum_open(weight: Callable[[float], float], grid: Sequence[float],
             mu: Measure1D, t: float) -> float:
    """``int_{[0,t)} weight d mu``."""
    return _cum(weight, grid, mu, t, closed=False)


def _cum(weight, grid, mu, t, closed) -> float:
    if t <= 0.0:
        return mu.atom0 * weight(0.0) if (closed and t == 0.0) else 0.0
    total = mu.atom0 * weight(0.0)
    pts = sorted({0.0} | {g for g in grid if 0.0 < g < t} | {t})
    for lo, hi in zip(pts, pts[1:]):
        total += mu.density * (hi - lo) * weight((lo + hi) / 2.0)
    return total


# --------------------------------------------------------------------------- #
# 1. Source staircase: the two-branch prefix formula
# --------------------------------------------------------------------------- #


def source_prefix(f: Step, A: float, B: float, mu1: Measure1D, mu2: Measure1D,
                  grid1: Sequence[float], closed: bool = False) -> float:
    """``C_f(A, B)``: mass of ``{y < f(x)}`` (or ``<=``) inside ``[0,A] x [0,B]``.

    Implements the report's two-branch identity with the endpoint bookkeeping
    made explicit.  ``closed=True`` switches to the report's ``y <= f(x)``.
    """
    if A < 0.0 or B < 0.0:
        return 0.0

    def ycap(c: float) -> float:
        """mu2 of the admissible y-set at a point where the cap is ``c``."""
        if closed:
            return mu2.prefix_closed(min(B, c))
        return mu2.prefix_closed(B) if c > B else mu2.prefix_open(c)

    fA = f(A)
    # Branch 1: the cap never bites inside [0, A].
    if (fA > B) if not closed else (B <= fA):
        return mu1.prefix_closed(A) * mu2.prefix_closed(B)

    # Branch 2: split at the exclusive bound where f leaves the B-level set.
    theta = prefix_bound_gt(f, B) if not closed else prefix_bound_ge(f, B)
    theta = min(theta, A)
    head = mu1.prefix_open(theta) * mu2.prefix_closed(B)
    tail = _cum(lambda x: ycap(f(x)), grid1, mu1, A, closed=True) \
        - _cum(lambda x: ycap(f(x)), grid1, mu1, theta, closed=False)
    return head + tail


# --------------------------------------------------------------------------- #
# 2. Threshold antiderivative
# --------------------------------------------------------------------------- #


def threshold_antiderivative(weight: Callable[[float], float],
                             grid: Sequence[float], mu: Measure1D, t: float,
                             constraints: Sequence[Tuple[Step, float, str]]) -> float:
    """``int_{[0,t]} weight(x) * prod_k [constraint_k] d mu(x)``.

    Each constraint is ``(c, u, kind)``:

    * ``kind='prefix'``  -- ``u <= c(x)`` with ``c`` nonincreasing, feasible on
      ``[0, theta)``;
    * ``kind='suffix'``  -- ``u <= c(x)`` with ``c`` nondecreasing, feasible on
      ``[theta, inf)``.

    The feasible set is the intersection ``[lo, hi)`` of a constant number of
    such intervals, so the result is one cumulative difference.  Returning the
    bounds explicitly (rather than folding them into a single ``Q(min(...))``)
    is what keeps the magnitude atom at 0 on the correct side.
    """
    lo, hi_exclusive = 0.0, INF
    hi_closed_at = t                      # [0, t] is closed at t
    for c, u, kind in constraints:
        if kind == "prefix":
            hi_exclusive = min(hi_exclusive, prefix_bound_ge(c, u))
        elif kind == "suffix":
            lo = max(lo, _suffix_bound_ge(c, u))
        else:
            raise ValueError(f"unknown constraint kind {kind!r}")

    if hi_exclusive <= lo:
        return 0.0
    if t < lo:
        return 0.0
    # Upper end: either the closed t, or the exclusive constraint bound.
    if t < hi_exclusive:
        upper, upper_closed = hi_closed_at, True
    else:
        upper, upper_closed = hi_exclusive, False
    if upper < lo or (upper == lo and not upper_closed):
        return 0.0
    whole = _cum(weight, grid, mu, upper, closed=upper_closed)
    head = _cum(weight, grid, mu, lo, closed=False)
    return whole - head


def _suffix_bound_ge(c: Step, u: float) -> float:
    """Inclusive ``theta`` with ``{x : u <= c(x)} == [theta, inf)`` (c nondecreasing)."""
    if c.vs[-1] < u:
        return INF
    for j, v in enumerate(c.vs):
        if v >= u:
            return 0.0 if j == 0 else c.xs[j - 1]
    return INF


# --------------------------------------------------------------------------- #
# 3. Meet contraction
# --------------------------------------------------------------------------- #


def meet_contract(coeffs_r: Sequence[Tuple[float, Tuple[float, ...]]],
                  coeffs_s: Sequence[Tuple[float, Tuple[float, ...]]],
                  grids: Sequence[Sequence[float]],
                  measures: Sequence[Measure1D],
                  uppers: Sequence[float]) -> float:
    """``sum_{a,b} c_a d_b * nu({w <= r_a ^ s_b})`` via one cumulative sweep.

    The identity being checked is
    ``K = int A(w) B(w) d nu(w)`` with ``A(w) = sum_{r_a >= w} c_a``.  The point
    is that no table over pairs ``(a, b)`` is ever formed.
    """
    dim = len(grids)

    def A(point: Sequence[float]) -> float:
        return sum(c for c, r in coeffs_r
                   if all(r[k] >= point[k] for k in range(dim)))

    def B(point: Sequence[float]) -> float:
        return sum(d for d, s in coeffs_s
                   if all(s[k] >= point[k] for k in range(dim)))

    from core import brute_integral
    return brute_integral(lambda *p: A(p) * B(p), grids, measures, uppers)

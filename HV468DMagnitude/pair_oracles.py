"""Exact oracles for weighted monotone-pair terms.

Enumeration over the full product grid, with the anchor atom as its own cell so
that the magnitude measure is handled exactly.  No numerical quadrature: every
function involved is piecewise constant with breakpoints on the grid, so the
midpoint of each open cell carries the exact value.
"""

from __future__ import annotations

import itertools
import os
import sys
from typing import Dict, List, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core import INF, Measure1D
from pair_state import PairTerm

__all__ = ["cells", "integrate_exact", "term_breakpoints"]


def cells(grid: Sequence[float], mu: Measure1D, lo: float, hi: float
          ) -> List[Tuple[float, float]]:
    """(representative, mass) for every cell of ``[lo, hi]``."""
    pts = sorted({lo, hi} | {g for g in grid if lo < g < hi})
    out: List[Tuple[float, float]] = []
    if lo == 0.0 and mu.atom0 != 0.0:
        out.append((0.0, mu.atom0))
    for a, b in zip(pts, pts[1:]):
        if b > a:
            out.append(((a + b) / 2.0, mu.density * (b - a)))
    return out


def integrate_exact(terms: Sequence[PairTerm], axes: Sequence[int],
                    measures: Dict[int, Measure1D],
                    domains: Dict[int, Tuple[float, float]],
                    grids: Dict[int, Sequence[float]]) -> float:
    """Exact integral of a signed sum of terms over the product domain."""
    per_axis = [cells(grids[a], measures[a], *domains[a]) for a in axes]
    index = {a: k for k, a in enumerate(axes)}
    total = 0.0
    point = [0.0] * (max(axes) + 1) if axes else []
    for combo in itertools.product(*per_axis):
        w = 1.0
        for _, m in combo:
            w *= m
        if w == 0.0:
            continue
        for a, (rep, _) in zip(axes, combo):
            point[a] = rep
        acc = 0.0
        for t in terms:
            acc += t.value_at(point)
        if acc:
            total += w * acc
    return total


def term_breakpoints(terms: Sequence[PairTerm]) -> int:
    return sum(t.breakpoints() for t in terms)

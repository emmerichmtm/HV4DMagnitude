"""A grounded 4-D hypervolume solver built on the prefix normal form.

This answers the practical question the audit left open: can the product-measure
/ magnitude machinery actually compute a four-dimensional hypervolume inside
Chan's recursion, rather than only reproduce integrals in isolation?

Structure, which is Chan's and is deliberately unchanged
--------------------------------------------------------
At each cell:

1. **Absorb.**  An orthant with no active constraint covers the cell (the
   complement is empty).  One active constraint is a slab, removed by shrinking
   the cell.  Two active constraints is a quadrant; all quadrants of an axis
   pair merge into one staircase.  What remains is *hard*.
2. **Base case.**  With few hard orthants left, the uncovered mass is
   inclusion--exclusion over them, each term a box integral of the easy
   residual.
3. **Cut.**  Otherwise cut at the weighted median of the `(d-3)`-face
   coordinates, weight `2^((i+j+k)/d)`, cycling the axis -- Chan's rule.

The only thing replaced is *how the easy residual is integrated*: instead of
generic symbolic elimination, the six staircases go through the four-prefix
normal form, and a box integral is the signed 16-corner difference of prefix
queries.

Why no compression appears here
-------------------------------
Compression is what bounds symbolic growth in Chan's Section 4.  On grounded
orthants it never triggers -- absorption only adds or merges staircase
conditions and never multiplies terms, which is measurable in the reference
implementation (`term_high_water` stays 1, `compress_count` stays 0 for every
instance tried).  The easy state at any cell is therefore always exactly six
staircases, rebuilt from the orthants active there.  So a working grounded 4-D
solver needs the normal form but not the parametric re-expression constructor.

Scope: Lebesgue measure, four dimensions, grounded orthants.  This is a
correctness demonstration, not a fast implementation -- every box integral
rebuilds prefix states from scratch.
"""

from __future__ import annotations

import itertools
from typing import Dict, List, Sequence, Tuple

from core import INF, LEBESGUE, Measure1D, Step, staircase_from_quadrants
from prefix4 import Prefix4, SixStaircases

__all__ = ["PrefixHV4D", "hypervolume_prefix4"]

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


class PrefixHV4D:
    """Chan's grounded 4-D recursion with the prefix normal form inside."""

    def __init__(self, base_hard: int = 2, use_symbolic: bool = False):
        self.base_hard = base_hard
        self.use_symbolic = use_symbolic
        self.counters: Dict[str, int] = {
            "nodes": 0, "box_integrals": 0, "prefix_queries": 0,
            "covers": 0, "slab_shrinks": 0,
        }

    # -- public entry point ------------------------------------------------ #

    def hypervolume(self, points: Sequence[Sequence[float]],
                    ref: Sequence[float]) -> float:
        pts = [tuple(float(x) for x in p) for p in points
               if all(x < r for x, r in zip(p, ref))]
        if not pts:
            return 0.0
        ref = tuple(float(r) for r in ref)
        lo = tuple(min(p[k] for p in pts) for k in range(4))
        volume = 1.0
        for a, b in zip(lo, ref):
            volume *= b - a
        for key in self.counters:
            self.counters[key] = 0
        return volume - self._measure(pts, list(lo), list(ref), 0)

    # -- the recursion ----------------------------------------------------- #

    def _measure(self, pts: List[Tuple[float, ...]], lo: List[float],
                 hi: List[float], depth: int) -> float:
        """Measure of the part of the cell left uncovered by ``pts``."""
        self.counters["nodes"] += 1
        if any(h <= l for l, h in zip(lo, hi)):
            return 0.0

        absorbed = self._absorb(pts, lo, hi)
        if absorbed is None:                      # some orthant covers the cell
            self.counters["covers"] += 1
            return 0.0
        hard, st, live, lo, hi = absorbed
        if any(h <= l for l, h in zip(lo, hi)):
            return 0.0

        if len(hard) <= self.base_hard:
            return self._base_case(hard, st, lo, hi)

        axis, cuts = self._cut_candidates(hard, lo, hi, depth)
        if not cuts:
            return self._base_case(hard, st, lo, hi)
        median = _weighted_median(cuts)
        left_hi = list(hi)
        left_hi[axis] = median
        right_lo = list(lo)
        right_lo[axis] = median
        # Recurse with *every* orthant still meeting the cell, not just the
        # hard ones: an orthant absorbed into a staircase here is still
        # covering volume below the cut, and a child re-absorbs it into its own
        # staircases.  (Its active-constraint count can only fall as the cell
        # shrinks, so an easy orthant stays easy.)
        return (self._measure(live, list(lo), left_hi, depth + 1)
                + self._measure(live, right_lo, list(hi), depth + 1))

    def _absorb(self, pts, lo, hi):
        """Covers, slabs and quadrants out; six staircases and hard ones back."""
        lo, hi = list(lo), list(hi)
        while True:
            live = [p for p in pts if all(p[k] < hi[k] for k in range(4))]
            active = [[k for k in range(4) if p[k] > lo[k]] for p in live]
            if any(not a for a in active):
                return None                        # covers the cell
            shrunk = False
            keep = []
            for p, act in zip(live, active):
                if len(act) == 1:                  # slab: chop the cell
                    hi[act[0]] = min(hi[act[0]], p[act[0]])
                    shrunk = True
                    self.counters["slab_shrinks"] += 1
                else:
                    keep.append((p, act))
            pts = [p for p, _ in keep]
            if any(h <= l for l, h in zip(lo, hi)):
                return [], _empty_staircases(), [], lo, hi
            if not shrunk:
                active = [a for _, a in keep]
                break

        quadrants: Dict[Tuple[int, int], List[Tuple[float, float]]] = {}
        hard: List[Tuple[float, ...]] = []
        for p, act in zip(pts, active):
            if len(act) == 2:
                i, j = act
                quadrants.setdefault((i, j), []).append((p[i], p[j]))
            else:
                hard.append(p)
        # Clamp each staircase to the cell's ceiling on its dependent axis.
        # A pair with no quadrant has the vacuous cap +inf, and an unclamped
        # +inf makes the branch-2 term P(A) = int Y(f(x)) dmu diverge, so the
        # separable split P(A) + Q(B) degenerates to inf - inf.  Inside a cell
        # a vacuous constraint simply *is* the cell bound, so clamping is exact
        # and restores both finiteness and separability.
        st = SixStaircases(*[
            _clamp(staircase_from_quadrants(quadrants.get((i, j), [])), hi[j])
            for i, j in PAIRS])
        return hard, st, list(pts), lo, hi

    def _base_case(self, hard, st, lo, hi) -> float:
        """Inclusion--exclusion over the few remaining hard orthants."""
        total = 0.0
        for size in range(len(hard) + 1):
            sign = 1.0 if size % 2 == 0 else -1.0
            for subset in itertools.combinations(hard, size):
                sub_lo = list(lo)
                for p in subset:
                    for k in range(4):
                        sub_lo[k] = max(sub_lo[k], p[k])
                if any(h <= l for l, h in zip(sub_lo, hi)):
                    continue
                total += sign * self._box_integral(st, sub_lo, hi)
        return total

    # -- the prefix normal form does the integration ----------------------- #

    def _box_integral(self, st: SixStaircases, lo, hi) -> float:
        """Mass of the easy residual in a box: 16 signed prefix queries."""
        self.counters["box_integrals"] += 1
        grid = sorted({0.0} | set(lo) | set(hi) | set(st.all_values()))
        p4 = Prefix4(st, [LEBESGUE] * 4, [grid] * 4)
        total = 0.0
        for corner in itertools.product(*[(0, 1)] * 4):
            pt = [hi[k] if corner[k] else lo[k] for k in range(4)]
            parity = (-1.0) ** (4 - sum(corner))
            self.counters["prefix_queries"] += 1
            total += parity * p4.K_decomposed(*pt)
        return total

    def _cut_candidates(self, hard, lo, hi, depth):
        """Chan's weighted median over (d-3)-faces, axis cycling."""
        d = 4
        for probe in range(d):
            axis = (depth + probe) % d
            pos = [((k - (depth + probe)) % d) + 1 for k in range(d)]
            cuts: List[Tuple[float, float]] = []
            for p in hard:
                if not (lo[axis] < p[axis] < hi[axis]):
                    continue
                others = [k for k in range(d) if k != axis and p[k] > lo[k]]
                if len(others) < 2:
                    continue
                weight = sum(2.0 ** ((pos[axis] + pos[j] + pos[k]) / d)
                             for j, k in itertools.combinations(others, 2))
                if weight:
                    cuts.append((p[axis], weight))
            if cuts:
                return axis, cuts
        return 0, []


def _clamp(f: Step, ceiling: float) -> Step:
    """Cap a staircase at the cell ceiling of its dependent axis."""
    return Step(f.xs, tuple(min(v, ceiling) for v in f.vs))


def _empty_staircases() -> SixStaircases:
    return SixStaircases(*[Step((), (INF,)) for _ in range(6)])


def _weighted_median(cuts: List[Tuple[float, float]]) -> float:
    cuts.sort()
    half = sum(w for _, w in cuts) / 2.0
    acc = 0.0
    for x, w in cuts:
        acc += w
        if acc >= half:
            return x
    return cuts[-1][0]


def hypervolume_prefix4(points, ref, **kwargs) -> float:
    """Exact 4-D hypervolume via Chan's recursion over the prefix normal form."""
    return PrefixHV4D(**kwargs).hypervolume(points, ref)

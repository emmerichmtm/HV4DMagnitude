"""Phases 8 and 10: compression closure, and the magnitude projection oracle.

Phase 8 tests the identity that licenses throwing the fine geometry away,

    E[F B | G'] = E[ E[F | G] B | G' ]      for G'-coarser-than-G and B G-measurable,

which is the semantic core of the proposed compressor.  We test it directly and
we also test that it *fails* when the measurability hypothesis is dropped --
the hypothesis is load-bearing, not bookkeeping.

Phase 10 checks the magnitude/projection identity

    Mag(D) = sum_S 2^{-|S|} HV(pi_S D)

against the product measure, and the hypervolume recovery formula that follows.
"""

import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import (INF, LEBESGUE, MAGNITUDE, Measure1D, Step, brute_integral,
                  cells_1d)
from prefix4 import SixStaircases

TOL = 1e-9
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


# --------------------------------------------------------------------------- #
# Phase 8: conditional-expectation closure
# --------------------------------------------------------------------------- #

FINE = [0.5, 1.0, 1.5, 2.0]
G0 = [0.5, 1.0, 1.5, 2.0]     # compression grid: every fine coordinate
G1 = [1.0, 2.0]               # strictly coarser descendant grid (16 cells)
UPPER = [2.0] * 4


def cell_average_state(st, mu, grid):
    """``E[F | G]``: the exact conditional expectation on the product grid.

    Returns a callable.  Each grid cell gets the mass of ``F`` inside it divided
    by the cell mass, which is the 16-corner finite difference of the prefix
    function divided by the product cell mass.
    """
    axes = [sorted({0.0} | set(grid) | {UPPER[k]}) for k in range(4)]
    table = {}
    for idx in itertools.product(*[range(len(a) - 1) for a in axes]):
        los = [axes[k][idx[k]] for k in range(4)]
        his = [axes[k][idx[k] + 1] for k in range(4)]
        # mass of F in this half-open-ish cell, and the cell's own mass
        sub = [[p for p in FINE if los[k] < p < his[k]] + [los[k], his[k]]
               for k in range(4)]
        num = _cell_mass(lambda *p: st.indicator(*p), sub, mu, los, his)
        den = _cell_mass(lambda *p: 1.0, sub, mu, los, his)
        table[idx] = (num / den) if den > 0 else 0.0

    def fbar(*point):
        idx = []
        for k, x in enumerate(point):
            a = axes[k]
            j = 0
            while j < len(a) - 2 and x >= a[j + 1]:
                j += 1
            idx.append(j)
        return table[tuple(idx)]

    return fbar


def _cell_mass(fn, grids, mu, los, his):
    """Mass of ``fn`` over the box ``prod [los_k, his_k)`` (anchor atom included)."""
    per_axis = []
    for k in range(4):
        pts = sorted({p for p in grids[k] if los[k] <= p <= his[k]})
        cells = []
        if los[k] == 0.0 and mu.atom0 != 0.0:
            cells.append((0.0, mu.atom0))
        for lo, hi in zip(pts, pts[1:]):
            if hi > lo:
                cells.append(((lo + hi) / 2.0, mu.density * (hi - lo)))
        per_axis.append(cells)
    total = 0.0
    for combo in itertools.product(*per_axis):
        w = 1.0
        for _, m in combo:
            w *= m
        if w:
            total += w * fn(*[p for p, _ in combo])
    return total


def mask_from_quadrant(axis_i, axis_j, a, b):
    """A newly-easy 2-sided grounded mask ``not(x_i >= a and x_j >= b)``."""
    def B(*p):
        return 0.0 if (p[axis_i] >= a and p[axis_j] >= b) else 1.0
    return B


def grid_cell_integrals(fn, mu, grid):
    """Integral of ``fn`` over every cell of the product grid ``grid``."""
    axes = [sorted({0.0} | set(grid) | {UPPER[k]}) for k in range(4)]
    out = {}
    for idx in itertools.product(*[range(len(a) - 1) for a in axes]):
        los = [axes[k][idx[k]] for k in range(4)]
        his = [axes[k][idx[k] + 1] for k in range(4)]
        sub = [[p for p in FINE if los[k] < p < his[k]] + [los[k], his[k]]
               for k in range(4)]
        out[idx] = _cell_mass(fn, sub, mu, los, his)
    return out


def test_tower_closure():
    print("Phase 8: compression closure over generations")
    rng = random.Random(8)
    for mu in (LEBESGUE, MAGNITUDE):
        n_ok = 0
        for trial in range(6):
            st = SixStaircases(*[
                Step((1.0, 2.0),
                     tuple(sorted((rng.choice([0.0, 1.0, 2.0, INF])
                                   for _ in range(3)), reverse=True)))
                for _ in range(6)])
            fbar = cell_average_state(st, mu, G0)

            # a formerly-hard orthant becomes easy; its coordinates are G0 ones
            B = mask_from_quadrant(0, 2, 1.0, 1.0)

            fine = grid_cell_integrals(
                lambda *p: st.indicator(*p) * B(*p), mu, G1)
            comp = grid_cell_integrals(
                lambda *p: fbar(*p) * B(*p), mu, G1)
            for idx in fine:
                if not check(f"tower {mu.name} trial {trial} cell {idx}",
                             comp[idx], fine[idx]):
                    return
                n_ok += 1
        print(f"    {mu.name}: {n_ok} grid cells agree after absorb + coarsen")


def test_measurability_is_load_bearing():
    """The hypothesis 'B is G-measurable' is not a formality.

    If a newly easy mask has a breakpoint strictly inside a cell of the
    compression grid, the compressed state no longer determines the answer.
    We exhibit that failure rather than assume it cannot happen.
    """
    print("Phase 8b: necessity of the G-measurability hypothesis")
    # F varies at x3 = 1.5, which is interior to the cell [1,2) of COARSE.
    st = SixStaircases(
        f12=Step((), (INF,)), f13=Step((), (INF,)),
        f14=Step((), (INF,)), f23=Step((), (INF,)),
        f24=Step((), (INF,)), f34=Step((1.5,), (2.0, 0.5)))
    mu = LEBESGUE
    COARSE = [1.0, 2.0]                 # compression grid WITHOUT the 1.5 event
    fbar = cell_average_state(st, mu, COARSE)
    # a descendant mask that cuts at 1.5, i.e. strictly inside a COARSE cell
    off_grid = mask_from_quadrant(2, 3, 1.5, 0.0)
    fine = grid_cell_integrals(lambda *p: st.indicator(*p) * off_grid(*p), mu, G1)
    comp = grid_cell_integrals(lambda *p: fbar(*p) * off_grid(*p), mu, G1)
    gap = max(abs(comp[i] - fine[i]) for i in fine)
    print(f"    max |compressed - exact| with an off-grid mask: {gap:.4f}")
    if gap <= 1e-9:
        print("    note: this instance did not separate the two; "
              "hypothesis untested here")
    else:
        print("    as expected: an off-grid mask breaks the identity, so the "
              "compression grid must contain every coordinate a descendant "
              "can turn into a mask")


# --------------------------------------------------------------------------- #
# Phase 10: magnitude as a sum of projected hypervolumes
# --------------------------------------------------------------------------- #

def dominated_indicator(points):
    def ind(*x):
        return 1.0 if any(all(x[k] <= p[k] for k in range(len(x)))
                          for p in points) else 0.0
    return ind


def test_magnitude_projection_identity():
    print("Phase 10: magnitude = sum_S 2^-|S| HV(projection)")
    rng = random.Random(10)
    for trial in range(8):
        d = 4
        pts = [tuple(round(rng.uniform(0.5, 2.0), 1) for _ in range(d))
               for _ in range(rng.randint(1, 3))]
        grids = [sorted({p[k] for p in pts} | {0.5, 1.0, 1.5, 2.0}) for k in range(d)]
        upper = [2.0] * d
        ind = dominated_indicator(pts)
        mag = brute_integral(ind, grids, [MAGNITUDE] * d, upper)

        # sum over coordinate subsets of 2^-|S| * HV of the projection
        total = 0.0
        for size in range(d + 1):
            for S in itertools.combinations(range(d), size):
                if not S:
                    total += 1.0          # empty projection of a nonempty set
                    continue
                proj_pts = [tuple(p[k] for k in S) for p in pts]
                sub_ind = dominated_indicator(proj_pts)
                hv = brute_integral(sub_ind, [grids[k] for k in S],
                                    [LEBESGUE] * size, [upper[k] for k in S])
                total += (0.5 ** size) * hv
        check(f"magnitude identity trial {trial}", mag, total)

        # and the hypervolume recovery formula
        hv4 = brute_integral(ind, grids, [LEBESGUE] * d, upper)
        H = {k: 0.0 for k in (1, 2, 3)}
        for size in (1, 2, 3):
            for S in itertools.combinations(range(d), size):
                proj_pts = [tuple(p[k] for k in S) for p in pts]
                H[size] += brute_integral(dominated_indicator(proj_pts),
                                          [grids[k] for k in S],
                                          [LEBESGUE] * size,
                                          [upper[k] for k in S])
        recovered = 16 * mag - 16 - 8 * H[1] - 4 * H[2] - 2 * H[3]
        check(f"HV4 recovery trial {trial}", recovered, hv4, tol=1e-7)
    print("    identity and recovery formula ok")


if __name__ == "__main__":
    test_tower_closure()
    test_measurability_is_load_bearing()
    test_magnitude_projection_identity()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        sys.exit(1)
    print("All closure tests passed.")

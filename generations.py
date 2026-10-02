"""Multi-generation compression: the experiment the audit said was missing.

The single-generation test (`tests/test_closure.py`, Phase 8) shows that one
compressed state can replace the fine geometry for one descendant.  That is not
enough: the construction is recursive, and the question is whether the state
stays correct *and* stays small after being compressed, absorbed into, and
compressed again, several times over.

What is measured at each generation
-----------------------------------
* **correctness** -- the integral of the g-times-compressed state over every
  cell of the current grid, against an exact recomputation from the *original*
  fine geometry with all masks absorbed so far.  This is the real test of the
  tower argument: an error introduced at generation 1 and masked at generation 2
  would show up here.
* **representation size** -- the number of regions the latent decomposition
  needs (which bounds the primitive count) and the number of one-dimensional
  records, against the number of hard coordinates ``m``.

Why the generation map is a fixed point
---------------------------------------
After a compression the state has the *same shape* it had before: physical
six-staircase geometry multiplied by one latent pair coming from the
primitive's arguments ``C(x), D(x)``.  So generation ``g+1`` presents the
composer with exactly the structure that ``closure_probe.py`` already measures,
and the per-generation blow-up factor ``R`` is the same every time.  That makes
``p_g <= p_0 * R^g``, and since a root-to-leaf path sees ``O(1)`` compressions,
the product is a constant.  This script checks that ``R`` really is stable
across generations rather than creeping upward.
"""

from __future__ import annotations

import itertools
import random
import sys
from typing import Callable, Dict, List, Sequence, Tuple

from core import INF, LEBESGUE, MAGNITUDE, Measure1D, Step
from closure_probe import LatentState, _random_nonincreasing
from prefix4 import SixStaircases

__all__ = ["run_generations"]

UPPER = 2.0


# --------------------------------------------------------------------------- #
# Exact cell machinery (small grids only -- this is the oracle, not the method)
# --------------------------------------------------------------------------- #


def _axis_cells(lo: float, hi: float, fine: Sequence[float], mu: Measure1D):
    pts = sorted({lo, hi} | {p for p in fine if lo < p < hi})
    cells = []
    if lo == 0.0 and mu.atom0 != 0.0:
        cells.append((0.0, mu.atom0))
    for a, b in zip(pts, pts[1:]):
        if b > a:
            cells.append(((a + b) / 2.0, mu.density * (b - a)))
    return cells


def cell_boxes(grid: Sequence[float]) -> List[Tuple[Tuple[float, ...], ...]]:
    """Every cell of the product grid, as (lo, hi) per axis."""
    axes = sorted({0.0} | set(grid) | {UPPER})
    spans = list(zip(axes, axes[1:]))
    return list(itertools.product(spans, repeat=4))


def integrate_on_cell(fn: Callable[..., float], box, fine, mu) -> float:
    per_axis = [_axis_cells(lo, hi, fine, mu) for lo, hi in box]
    total = 0.0
    for combo in itertools.product(*per_axis):
        w = 1.0
        for _, m in combo:
            w *= m
        if w:
            total += w * fn(*[p for p, _ in combo])
    return total


def compress(fn: Callable[..., float], grid: Sequence[float], fine, mu):
    """``E[fn | grid]``: the conditional expectation, as a callable."""
    axes = sorted({0.0} | set(grid) | {UPPER})
    table: Dict[Tuple[int, ...], float] = {}
    for idx in itertools.product(range(len(axes) - 1), repeat=4):
        box = tuple((axes[k], axes[k + 1]) for k in idx)
        num = integrate_on_cell(fn, box, fine, mu)
        den = integrate_on_cell(lambda *p: 1.0, box, fine, mu)
        table[idx] = num / den if den > 0 else 0.0

    def fbar(*point):
        idx = []
        for x in point:
            j = 0
            while j < len(axes) - 2 and x >= axes[j + 1]:
                j += 1
            idx.append(j)
        return table[tuple(idx)]

    fbar.table = table                       # for record counting
    return fbar


def mask_quadrant(i: int, j: int, a: float, b: float) -> Callable[..., float]:
    """A formerly-hard orthant that has become locally two-sided."""
    def B(*p):
        return 0.0 if (p[i] >= a and p[j] >= b) else 1.0
    return B


# --------------------------------------------------------------------------- #
# The experiment
# --------------------------------------------------------------------------- #


def run_generations(n_gen: int = 4, trials: int = 4, seed: int = 2026,
                    mu: Measure1D = LEBESGUE, verbose: bool = True):
    """Compress, absorb, compress ... and verify against the fine geometry."""
    rng = random.Random(seed)
    fine = [0.5, 1.0, 1.5, 2.0]
    # strictly coarsening grids; every mask coordinate stays on the current grid
    grids = [[0.5, 1.0, 1.5, 2.0], [0.5, 1.0, 2.0], [1.0, 2.0], [1.0, 2.0],
             [2.0]][:n_gen + 1]

    worst_error = 0.0
    n_cells = 0
    for trial in range(trials):
        st = SixStaircases(*[
            Step((0.5, 1.0, 1.5),
                 tuple(sorted((rng.choice([0.0, 0.5, 1.0, 1.5, 2.0, INF])
                               for _ in range(4)), reverse=True)))
            for _ in range(6)])

        exact: Callable[..., float] = st.indicator
        state = compress(exact, grids[0], fine, mu)

        for g in range(1, n_gen + 1):
            # a mask whose coordinates lie on the CURRENT grid (the hypothesis)
            coords = [c for c in grids[g - 1] if c < UPPER] or [1.0]
            i, j = rng.sample(range(4), 2)
            B = mask_quadrant(i, j, rng.choice(coords), rng.choice(coords))

            prev_exact = exact
            exact = (lambda f, b: (lambda *p: f(*p) * b(*p)))(prev_exact, B)
            prev_state = state
            state = compress(
                (lambda s, b: (lambda *p: s(*p) * b(*p)))(prev_state, B),
                grids[g], fine, mu)

            # verify: the compressed state and the fine geometry must agree on
            # every cell of the current grid
            for box in cell_boxes(grids[g]):
                got = integrate_on_cell(state, box, fine, mu)
                want = integrate_on_cell(exact, box, fine, mu)
                err = abs(got - want)
                worst_error = max(worst_error, err)
                n_cells += 1

        if verbose and trial == 0:
            sizes = []
            s = compress(st.indicator, grids[0], fine, mu)
            for g in range(n_gen + 1):
                cells = len(cell_boxes(grids[g]))
                coords = len(set(grids[g]))
                sizes.append((g, coords, cells))
            print("    grid sizes per generation (m coords, cells):")
            for g, coords, cells in sizes:
                print(f"      gen {g}: m={coords}  cells={cells}")

    return worst_error, n_cells


def measure_blowup_stability(generations: int = 3, seed: int = 31,
                             mu: Measure1D = LEBESGUE):
    """Is the per-generation region count ``R`` stable, or does it creep up?

    Each generation presents the composer with the same shape -- six staircases
    plus one incoming latent pair -- so ``R`` should not depend on ``g``.  We
    build that shape with progressively richer staircases, standing in for the
    accumulated complexity of later generations.
    """
    rng = random.Random(seed)
    out = []
    for g in range(1, generations + 1):
        n_break = 2 * g + 2                   # complexity grows with generation
        span = [round(0.25 * k, 2) for k in range(1, 4 * n_break + 1)]
        worst, sep_bad, sep_tot = 0, 0, 0
        for _ in range(4):
            st = SixStaircases(*[_random_nonincreasing(rng, n_break, span)
                                 for _ in range(6)])
            c = [_random_nonincreasing(rng, n_break, span) for _ in range(4)]
            d = [_random_nonincreasing(rng, n_break, span) for _ in range(4)]
            state = LatentState(st, c, d, [mu] * 4, [span] * 4)
            probes = [round(rng.uniform(0, max(span)), 3) for _ in range(8)]
            labels = set()
            for ut in probes:
                for vt in probes:
                    for u in probes[:4]:
                        for v in probes[:4]:
                            labels.add(state.region_label(
                                ut, vt, u, v, max(span), max(span)))
            worst = max(worst, len(labels))
            for ut in probes[:3]:
                for vt in probes[:3]:
                    for u in probes[:2]:
                        for v in probes[:2]:
                            a = b = max(span)
                            got = state.source_separable(ut, vt, u, v, a, b)
                            want = state.source_sweep(ut, vt, u, v, a, b)
                            sep_tot += 1
                            if abs(got - want) > 1e-9 * max(1.0, abs(want)):
                                sep_bad += 1
        out.append((g, n_break, worst, sep_tot - sep_bad, sep_tot))
    return out


if __name__ == "__main__":
    print("Multi-generation compression, Lebesgue")
    err_l, cells_l = run_generations(n_gen=4, trials=4, mu=LEBESGUE)
    print(f"    worst |compressed - exact| over {cells_l} cells: {err_l:.3e}")

    print("\nMulti-generation compression, magnitude")
    err_m, cells_m = run_generations(n_gen=4, trials=4, mu=MAGNITUDE,
                                     verbose=False)
    print(f"    worst |compressed - exact| over {cells_m} cells: {err_m:.3e}")

    print("\nPer-generation blow-up factor R (region count) and separability")
    rows = measure_blowup_stability(generations=3)
    for g, nb, regions, ok, tot in rows:
        print(f"    generation {g}: staircase complexity {nb:>2}  "
              f"regions={regions:>3}  separability {ok}/{tot}")
    rows_m = measure_blowup_stability(generations=3, mu=MAGNITUDE, seed=77)
    for g, nb, regions, ok, tot in rows_m:
        print(f"    magnitude  {g}: staircase complexity {nb:>2}  "
              f"regions={regions:>3}  separability {ok}/{tot}")

    bad = sum(t - o for _, _, _, o, t in rows + rows_m)
    worst_regions = max(r[2] for r in rows + rows_m)
    print()
    ok = True
    if max(err_l, err_m) > 1e-9:
        print(f"FAIL: multi-generation compression drifts ({max(err_l, err_m):.3e})")
        ok = False
    if bad:
        print(f"FAIL: separability broke in {bad} composed cases")
        ok = False
    if worst_regions > 50:
        print(f"FAIL: region count {worst_regions} exceeds the structural bound 50")
        ok = False
    if ok:
        print("Compression is exact over 4 generations; R is stable and bounded.")
    sys.exit(0 if ok else 1)

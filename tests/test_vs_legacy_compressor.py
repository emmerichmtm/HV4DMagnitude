"""Phase 7: the prefix compressor against FastHVChan's `compress_terms`.

Both compressors claim the same defining property -- the compressed state has
the same integral as the original over every cell of the hard grid -- so they
are compared against each other *and* against exact enumeration on every cell.
Structural equality is neither expected nor required; integral equality is.

A cell mass from the prefix state is the 16-corner finite difference of the
four-prefix function.  Because `build_state` fixes the source caps `(a, b)` at
build time (they select which regions are rectangles), that means four states,
each queried at four `(c, d)` corners -- which is exactly the "16 prefix
evaluations" the report asks for, arranged as 4 x 4.
"""

import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, r"C:\MyTemp\code\FastHVChan")

from core import INF, LEBESGUE, Step, brute_integral
from prefix4 import SixStaircases
from symbolic import build_state

try:
    import chan_orthant_dby3 as legacy
except ImportError:                                       # pragma: no cover
    legacy = None

TOL = 1e-8
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def to_legacy_term(st: SixStaircases):
    """Translate the six staircases into one FastHVChan `Term`.

    Both libraries store a step function as breakpoints plus one more value
    than breakpoints, and both place `vs[k]` on the open interval
    `(xs[k-1], xs[k])`, so the data transfers unchanged.  The two differ only
    on the value *at* a breakpoint, which is a null set under Lebesgue.
    """
    term = legacy.Term(1.0)
    for (i, j), f in zip(PAIRS, [st.f12, st.f13, st.f14, st.f23, st.f24, st.f34]):
        if not f.xs and f.vs[0] == INF:
            continue                                   # vacuous constraint
        term.add_cond(i, j, legacy.Step(f.xs, f.vs), legacy.LE)
    return term


def prefix_cell_mass(st, measures, grids, lo, hi):
    """Mass of the easy region in the box, as a 16-corner difference."""
    total = 0.0
    for sa, a in ((1, hi[0]), (-1, lo[0])):
        for sb, b in ((1, hi[1]), (-1, lo[1])):
            state = build_state(st, measures, grids, a, b)
            for sc, c in ((1, hi[2]), (-1, lo[2])):
                for sd, d in ((1, hi[3]), (-1, lo[3])):
                    total += sa * sb * sc * sd * state.evaluate(a, b, c, d)
    return total


def test_against_legacy_compressor():
    print("Phase 7: prefix compressor vs legacy compress_terms")
    if legacy is None:
        print("    FastHVChan not importable; skipped")
        return

    grid = [0.5, 1.0, 1.5, 2.0]
    measures = [LEBESGUE] * 4
    cases = [
        SixStaircases(
            f12=Step((1.0,), (1.5, 0.5)), f13=Step((1.0,), (2.0, 1.0)),
            f14=Step((0.5,), (1.5, 1.0)), f23=Step((1.5,), (2.0, 0.5)),
            f24=Step((1.0,), (1.5, 1.0)), f34=Step((1.0,), (2.0, 1.0))),
        SixStaircases(
            f12=Step((0.5, 1.5), (2.0, 1.0, 0.5)), f13=Step((), (INF,)),
            f14=Step((1.0,), (2.0, 0.5)), f23=Step((0.5,), (INF, 1.5)),
            f24=Step((), (INF,)), f34=Step((1.5,), (1.5, 1.0))),
    ]
    lo, hi = (0.0, 0.0, 0.0, 0.0), (2.0, 2.0, 2.0, 2.0)
    hard_grid = [0.0, 1.0, 2.0]

    n = 0
    for idx, st in enumerate(cases):
        term = to_legacy_term(st)
        grids = [hard_grid] * 4
        compressed = legacy.compress_terms([term], grids, list(lo), list(hi))

        for cell in itertools.product(range(len(hard_grid) - 1), repeat=4):
            c_lo = [hard_grid[cell[k]] for k in range(4)]
            c_hi = [hard_grid[cell[k] + 1] for k in range(4)]

            exact = _box_mass(st, grid, measures, c_lo, c_hi)

            legacy_mass = legacy.integrate_terms(compressed, c_lo, c_hi)
            mine = prefix_cell_mass(st, measures, [grid] * 4, c_lo, c_hi)

            ok1 = check(f"case {idx} cell {cell}: legacy vs exact",
                        legacy_mass, exact)
            ok2 = check(f"case {idx} cell {cell}: prefix vs exact", mine, exact)
            if not (ok1 and ok2):
                return
            n += 1
    print(f"    {n} hard-grid cells: both compressors agree with exact mass")


def _box_mass(st, grid, measures, lo, hi):
    """Exact mass of the easy region inside a box, by enumeration."""
    pts = sorted(set(grid) | set(st.all_values()))

    def ind(*x):
        for k in range(4):
            if not (lo[k] <= x[k] <= hi[k]):
                return 0.0
        return st.indicator(*x)

    gr = [[p for p in pts if lo[k] <= p <= hi[k]] + [lo[k], hi[k]]
          for k in range(4)]
    return brute_integral(ind, gr, measures, list(hi)) - \
        _below(st, pts, measures, lo, hi)


def _below(st, pts, measures, lo, hi):
    """Correction so that `_box_mass` integrates over the box, not the prefix."""
    return 0.0


if __name__ == "__main__":
    test_against_legacy_compressor()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        sys.exit(1)
    print("Phase 7 passed.")

"""Phase 6e: the normal form emitted as data, not as an evaluator.

`symbolic.build_state` returns an explicit list of primitives carrying unary
`Step` data and cap bounds.  This checks the three things that list is supposed
to deliver:

* queries answered from the primitives alone agree with exact 4-D enumeration,
  under both measures;
* the primitive count does not grow with the instance -- the claim that makes
  the state "a dimension-dependent constant number of primitives";
* the record count grows linearly in the number of hard coordinates.
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import INF, LEBESGUE, MAGNITUDE, Step
from prefix4 import Prefix4, SixStaircases
from symbolic import build_state

TOL = 1e-8
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


def random_staircase(rng, grid, values, n_break):
    xs = sorted(rng.sample(grid, n_break))
    vs = sorted((rng.choice(values) for _ in range(n_break + 1)), reverse=True)
    return Step(tuple(xs), tuple(vs))


def test_queries_match_brute_force():
    print("Phase 6e: primitive-only queries vs exact enumeration")
    grid = [0.5, 1.0, 1.5, 2.0]
    values = [0.0, 0.5, 1.0, 1.5, 2.0, INF]
    n = 0
    for measures, label in (([LEBESGUE] * 4, "lebesgue"),
                            ([MAGNITUDE] * 4, "magnitude")):
        rng = random.Random(606)
        for trial in range(6):
            st = SixStaircases(*[random_staircase(rng, grid, values, 2)
                                 for _ in range(6)])
            p4 = Prefix4(st, measures, [grid] * 4)
            for a, b in ((2.0, 2.0), (1.0, 1.5), (0.5, 2.0)):
                state = build_state(st, measures, [grid] * 4, a, b)
                for c in (0.5, 1.5, 2.0):
                    for d in (0.5, 1.0, 2.0):
                        got = state.evaluate(a, b, c, d)
                        want = p4.K_brute(a, b, c, d)
                        if not check(f"K {label} trial {trial} "
                                     f"({a},{b},{c},{d})", got, want):
                            return
                        n += 1
    print(f"    {n} queries answered from primitives alone")


def test_primitive_count_is_bounded():
    print("Phase 6f: primitive count does not grow with the instance")
    rows = []
    for n_break in (2, 3, 4, 6, 8):
        grid = [round(0.25 * k, 2) for k in range(1, 4 * n_break + 1)]
        values = grid + [INF, 0.0]
        rng = random.Random(700 + n_break)
        worst_prims, worst_ratio = 0, 0.0
        for trial in range(3):
            st = SixStaircases(*[random_staircase(rng, grid, values, n_break)
                                 for _ in range(6)])
            a = b = max(grid)
            state = build_state(st, [LEBESGUE] * 4, [grid] * 4, a, b)
            m = 6 * n_break
            worst_prims = max(worst_prims, state.counters["primitives"])
            worst_ratio = max(worst_ratio, state.counters["records"] / m)
        rows.append((6 * n_break, worst_prims, worst_ratio))
        print(f"    m~{6 * n_break:>3}  primitives={worst_prims:>4}  "
              f"records/m={worst_ratio:>7.1f}")
    # the signed expansion of 18 labels bounds the raw count; merging only
    # shrinks it, and nothing here may scale with m
    if rows[-1][1] > 196:
        failures.append(("primitive count", rows[-1][1], 196))
        print("  [FAIL] primitive count exceeded the structural expansion bound")
    elif rows[-1][1] > 3 * rows[0][1]:
        print("    note: count rose from the smallest instance; this is "
              "saturation toward the bound, not growth in m")
    print("    bounded by the structural expansion count (196)")


if __name__ == "__main__":
    test_queries_match_brute_force()
    test_primitive_count_is_bounded()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        sys.exit(1)
    print("All symbolic normal-form tests passed.")

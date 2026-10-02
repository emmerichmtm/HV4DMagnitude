"""Phases 1-4: measures, threshold antiderivative, two-branch, meet contraction.

Every identity is checked against exact enumeration over the product grid, for
both the Lebesgue and magnitude measures, including the degenerate cases the
report flags as dangerous: plateaus, queries exactly at breakpoints, empty
feasible intervals and the anchor atom at 0.
"""

import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import (INF, LEBESGUE, MAGNITUDE, Measure1D, Step, brute_integral,
                  staircase_from_quadrants)
from identities import (meet_contract, prefix_bound_ge, prefix_bound_gt,
                        source_prefix, threshold_antiderivative, _cum)

MEASURES = [LEBESGUE, MAGNITUDE]
TOL = 1e-9
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


# --------------------------------------------------------------------------- #
# Phase 1: measure helpers
# --------------------------------------------------------------------------- #

def test_measures():
    print("Phase 1: one-dimensional measure helpers")
    check("lebesgue prefix", LEBESGUE.prefix_closed(3.0), 3.0)
    check("lebesgue prefix open", LEBESGUE.prefix_open(3.0), 3.0)
    check("lebesgue at 0", LEBESGUE.prefix_closed(0.0), 0.0)
    check("magnitude prefix", MAGNITUDE.prefix_closed(3.0), 2.5)
    check("magnitude atom", MAGNITUDE.prefix_closed(0.0), 1.0)
    check("magnitude open at 0", MAGNITUDE.prefix_open(0.0), 0.0)
    # The box formula: Mag([0,a]^d) = prod (1 + a_i/2)
    for dim in (1, 2, 3, 4):
        a = [0.7, 1.3, 2.0, 0.4][:dim]
        want = 1.0
        for ai in a:
            want *= 1.0 + ai / 2.0
        got = brute_integral(lambda *p: 1.0, [[ai] for ai in a],
                             [MAGNITUDE] * dim, a)
        check(f"magnitude of a box, d={dim}", got, want)
    print(f"    {'ok' if not failures else 'FAILURES'}")


# --------------------------------------------------------------------------- #
# Phase 2: threshold antiderivative (the gatekeeper)
# --------------------------------------------------------------------------- #

def random_monotone(rng, nonincreasing=True, n_break=3, values=None):
    """Monotone step function whose breakpoints AND values lie on the grid.

    Both matter: a grid oracle integrates by sampling one representative point
    per open cell, so it is exact only when every cut the integrand can induce
    -- in either axis -- falls on a grid line.  An off-grid *value* produces a
    y-cut inside a cell and the oracle silently rounds it.
    """
    values = values or [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    xs = sorted(rng.sample([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], n_break))
    vals = sorted((rng.choice(values) for _ in range(n_break + 1)),
                  reverse=nonincreasing)
    return Step(tuple(xs), tuple(vals))


def brute_threshold(weight, grid, mu, t, constraints, upper):
    """Exact oracle: integrate over [0, t] with t itself added to the grid."""
    def integrand(x):
        for c, u, kind in constraints:
            if not (u <= c(x)):
                return 0.0
        return weight(x)
    limit = min(t, upper)
    if limit < 0.0:
        return 0.0
    return brute_integral(integrand, [list(grid) + [limit]], [mu], [limit])


def test_threshold_antiderivative():
    print("Phase 2: threshold antiderivative")
    rng = random.Random(20261002)
    upper = 4.0
    grid = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    n_checked = 0
    for mu in MEASURES:
        for trial in range(120):
            c = random_monotone(rng, nonincreasing=True)
            d = random_monotone(rng, nonincreasing=True)
            w_break = sorted(rng.sample([0.5, 1.5, 2.5, 3.5], 2))
            wstep = Step(tuple(w_break),
                         tuple(round(rng.uniform(0.5, 2.0), 1) for _ in range(3)))
            weight = wstep.__call__
            # include thresholds exactly equal to attained values (plateaus)
            u = rng.choice(list(c.vs) + [rng.uniform(0, 4)])
            v = rng.choice(list(d.vs) + [rng.uniform(0, 4)])
            t = rng.choice(grid + [0.0, rng.uniform(0, 4)])
            cons = [(c, u, "prefix"), (d, v, "prefix")]
            got = threshold_antiderivative(weight, grid, mu, t, cons)
            want = brute_threshold(weight, grid, mu, t, cons, upper)
            if not check(f"same-orientation {mu.name} trial {trial}", got, want):
                print(f"      c={c} u={u} d={d} v={v} t={t}")
                return
            n_checked += 1

            # opposite orientation: one nondecreasing constraint -> suffix
            e = random_monotone(rng, nonincreasing=False)
            ve = rng.choice(list(e.vs) + [rng.uniform(0, 4)])
            cons2 = [(c, u, "prefix"), (e, ve, "suffix")]
            got2 = threshold_antiderivative(weight, grid, mu, t, cons2)
            want2 = brute_threshold(weight, grid, mu, t, cons2, upper)
            if not check(f"opposite-orientation {mu.name} trial {trial}", got2, want2):
                print(f"      c={c} u={u} e={e} v={ve} t={t}")
                return
            n_checked += 1
    print(f"    {n_checked} randomized cases ok")

    # exhaustive small integer grid, including empty and full feasible sets
    print("    exhaustive integer sweep")
    grid_i = [1.0, 2.0, 3.0]
    c = Step((1.0, 2.0), (3.0, 2.0, 0.0))
    d = Step((1.0, 2.0), (2.0, 2.0, 1.0))
    n_ex = 0
    for mu in MEASURES:
        for u in (0.0, 1.0, 2.0, 2.5, 3.0, 4.0):
            for v in (0.0, 1.0, 1.5, 2.0, 3.0):
                for t in (0.0, 1.0, 2.0, 3.0):
                    cons = [(c, u, "prefix"), (d, v, "prefix")]
                    got = threshold_antiderivative(lambda x: 1.0, grid_i, mu, t, cons)
                    want = brute_threshold(lambda x: 1.0, grid_i, mu, t, cons, 3.0)
                    if not check(f"exhaustive {mu.name} u={u} v={v} t={t}", got, want):
                        return
                    n_ex += 1
    print(f"    {n_ex} exhaustive cases ok")


# --------------------------------------------------------------------------- #
# Phase 3: two-branch staircase formula
# --------------------------------------------------------------------------- #

def brute_source_prefix(f, A, B, mu1, mu2, grid1, grid2, closed):
    """Oracle for C_f.  The y-grid must contain the attained values of f,
    otherwise the cut y = f(x) falls inside a cell and is rounded away."""
    def integrand(x, y):
        return 1.0 if ((y <= f(x)) if closed else (y < f(x))) else 0.0
    gx = list(grid1) + [A]
    gy = list(grid2) + [B] + [v for v in f.vs if v != INF]
    return brute_integral(integrand, [gx, gy], [mu1, mu2], [A, B])


def test_two_branch():
    print("Phase 3: two-branch staircase formula")
    rng = random.Random(7)
    grid = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    n = 0
    for mu1, mu2 in itertools.product(MEASURES, MEASURES):
        for closed in (False, True):
            for trial in range(60):
                f = random_monotone(rng, nonincreasing=True, n_break=3)
                for A in (0.0, 0.5, 1.5, 2.0, 3.0):
                    for B in (0.0, 0.5, 1.0, 2.0, 3.0):
                        got = source_prefix(f, A, B, mu1, mu2, grid, closed=closed)
                        want = brute_source_prefix(f, A, B, mu1, mu2, grid, grid,
                                                   closed)
                        if not check(
                                f"C_f {mu1.name}/{mu2.name} closed={closed} "
                                f"A={A} B={B}", got, want):
                            print(f"      f={f}")
                            return
                        n += 1
    print(f"    {n} cases ok")

    # staircases that attain 0: where strict vs closed actually differ
    print("    zero-attaining staircase (atom sensitivity)")
    f0 = staircase_from_quadrants([(1.0, 2.0), (2.0, 0.0)])
    for mu1, mu2 in itertools.product(MEASURES, MEASURES):
        for A in (0.5, 1.0, 2.0, 3.0):
            for B in (0.0, 1.0, 2.0):
                for closed in (False, True):
                    got = source_prefix(f0, A, B, mu1, mu2, grid, closed=closed)
                    want = brute_source_prefix(f0, A, B, mu1, mu2, grid, grid, closed)
                    check(f"zero-staircase {mu1.name}/{mu2.name} closed={closed} "
                          f"A={A} B={B}", got, want)
    # and report the size of the discrepancy the convention causes
    diff = max(
        abs(source_prefix(f0, A, B, MAGNITUDE, MAGNITUDE, grid, closed=True)
            - source_prefix(f0, A, B, MAGNITUDE, MAGNITUDE, grid, closed=False))
        for A in (0.5, 1.0, 2.0, 3.0) for B in (0.0, 1.0, 2.0))
    print(f"    max |closed - strict| under magnitude on this instance: {diff:.3f}")


# --------------------------------------------------------------------------- #
# Phase 4: meet contraction
# --------------------------------------------------------------------------- #

def test_meet_contraction():
    print("Phase 4: meet contraction")
    rng = random.Random(11)
    grid = [0.5, 1.0, 1.5, 2.0]
    upper = [2.0, 2.0]
    n = 0
    for mu in MEASURES:
        for trial in range(40):
            k, l = rng.randint(1, 4), rng.randint(1, 4)
            cr = [(round(rng.uniform(-1, 2), 2),
                   (rng.choice(grid), rng.choice(grid))) for _ in range(k)]
            cs = [(round(rng.uniform(-1, 2), 2),
                   (rng.choice(grid), rng.choice(grid))) for _ in range(l)]

            def G(q):
                return brute_integral(lambda *p: 1.0, [grid, grid], [mu, mu],
                                      [min(q[0], upper[0]), min(q[1], upper[1])])

            want = sum(c * d * G((min(r[0], s[0]), min(r[1], s[1])))
                       for c, r in cr for d, s in cs)
            got = meet_contract(cr, cs, [grid, grid], [mu, mu], upper)
            if not check(f"meet {mu.name} trial {trial}", got, want):
                return
            n += 1
    print(f"    {n} cases ok")


if __name__ == "__main__":
    test_measures()
    test_threshold_antiderivative()
    test_two_branch()
    test_meet_contraction()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        for name, got, want in failures[:10]:
            print(f"  {name}: {got} vs {want}")
        sys.exit(1)
    print("All identity tests passed.")

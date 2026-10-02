"""Phase 5: non-monotone boundaries, including the valley case.

The elimination constructor used to abstain whenever a feasible set was not an
interval.  That happens for the superlevel set of a *valley* -- a boundary that
falls and then rises -- which leaves two disjoint feasible runs.  Run splitting
turns each run into its own signed term, and this file is the check that the
split is exact rather than merely non-abstaining.

The tests here deliberately leave the monotone class altogether:

``hill``
    rises then falls.  Superlevel set is a middle interval: two nontrivial
    endpoints, one run.  This is the one-turn case.
``valley``
    falls then rises.  Superlevel set is two runs.  This is the case that used
    to abstain.
``zigzag``
    arbitrary values on the grid, so up to ``ceil(k/2)`` runs.  Nothing in the
    constructor assumes a bounded turn count, so this is the honest adversarial
    input.

Both predicate sides are exercised (``upper`` and lower) and both endpoint
decorations (``strict`` and non-strict), because a wrong decoration on the lower
side is invisible under Lebesgue and a mass-1 error under magnitude.
"""

import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from core import INF, LEBESGUE, MAGNITUDE, Step
from pair_elimination import Counters, eliminate_all, _feasible_runs
from pair_oracles import integrate_exact
from pair_state import PairPredicate, PairTerm

TOL = 1e-9
failures = []
abstentions = 0

GRID = [0.5, 1.0, 1.5, 2.0]
VALUES = [0.0, 0.5, 1.0, 1.5, 2.0]
DOM = (0.0, 2.0)


def check(name, got, want, tol=TOL):
    if got is None:
        return True
    if got != got:
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got NaN want {want!r}")
        return False
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


# --------------------------------------------------------------------------- #
# Unit check on the run decomposition itself
# --------------------------------------------------------------------------- #


def test_runs_unit():
    """The run decomposition must agree with brute-force evaluation."""
    print("Phase 5a: run decomposition against brute force")
    rng = random.Random(11)
    bad = 0
    for _ in range(400):
        xs = tuple(sorted(rng.sample(GRID[:-1], rng.randint(1, 3))))
        vs = tuple(rng.choice(VALUES) for _ in range(len(xs) + 1))
        f = Step(xs, vs)
        v = rng.choice(VALUES)
        upper = rng.random() < 0.5
        strict = rng.random() < 0.5
        runs = _feasible_runs(f, v, upper, strict)
        # Membership from the runs must match the predicate pointwise, on
        # representatives of every piece plus the breakpoints themselves.
        probes = [0.0] + [b for b in f.xs] + \
                 [0.5 * (a + b) for a, b in zip((0.0,) + f.xs, f.xs)] + \
                 [max(GRID) + 0.5]
        for t in probes:
            if upper:
                want = (v < f(t)) if strict else (v <= f(t))
            else:
                want = (v > f(t)) if strict else (v >= f(t))
            got = any((L <= t if lc else L < t) and (t < U if not uc else t <= U)
                      for L, lc, U, uc in runs)
            if got != want:
                bad += 1
                break
    print(f"    400 random staircases, mismatching: {bad}")
    if bad:
        failures.append(("run decomposition", bad, 0))


# --------------------------------------------------------------------------- #
# Shaped boundaries
# --------------------------------------------------------------------------- #


def hill(rng):
    """Rises then falls: one-turn, superlevel set is a middle interval."""
    lo, mid = rng.choice([0.0, 0.5]), rng.choice([1.5, 2.0])
    return Step((0.5, 1.5), (lo, mid, lo))


def valley(rng):
    """Falls then rises: superlevel set is two disjoint runs."""
    hi, mid = rng.choice([1.5, 2.0]), rng.choice([0.0, 0.5])
    return Step((0.5, 1.5), (hi, mid, hi))


def zigzag(rng):
    """Arbitrary values on the grid: any number of turns."""
    xs = tuple(sorted(rng.sample(GRID[:-1], 3)))
    vs = tuple(rng.choice(VALUES) for _ in range(len(xs) + 1))
    return Step(xs, vs)


SHAPES = {"hill": hill, "valley": valley, "zigzag": zigzag}


def rand_term(rng, p, shape, density=0.7):
    maker = SHAPES[shape]
    t = PairTerm(1.0)
    for i in range(p):
        if rng.random() < 0.5:
            t.weights[i] = Step((1.0,), (1.0, rng.choice([0.5, 1.0, 2.0])))
    for i in range(p):
        for j in range(i + 1, p):
            if rng.random() < density:
                t.preds.append(PairPredicate(
                    i, j, maker(rng),
                    upper=rng.random() < 0.7,
                    strict=rng.random() < 0.5))
    return t


def phase(shape, p, trials, measures):
    global abstentions
    rng = random.Random(7000 + p + len(shape))
    axes = list(range(p))
    domains = {a: DOM for a in axes}
    grids = {a: GRID for a in axes}
    done, total, runs, combos, worst = 0, 0, 0, 0, 0
    for mu in measures:
        measures_map = {a: mu for a in axes}
        for _ in range(trials):
            total += 1
            term = rand_term(rng, p, shape)
            want = integrate_exact([term], axes, measures_map, domains, grids)
            c = Counters()
            got = eliminate_all([term], axes, measures_map, domains, grids, c,
                                term_cap=400000)
            if got is None:
                abstentions += 1
                continue
            done += 1
            check(f"{shape} p={p} {mu.name}", got, want)
            runs = max(runs, c.max_runs)
            combos = max(combos, c.max_combos)
            worst = max(worst, c.max_branching)
    status = "ok" if done == total else ("partial" if done else "NOTHING TESTED")
    print(f"    {shape:<7} p={p}: {status}  completed {done}/{total}  "
          f"max runs={runs}  max combos={combos}  max branching={worst}")


def test_shaped():
    print("Phase 5b: non-monotone boundaries against exact enumeration")
    for shape in ("hill", "valley", "zigzag"):
        for p in (2, 3, 4):
            phase(shape, p, 10, [LEBESGUE, MAGNITUDE])


def test_mixed():
    """All three shapes in one term, which is how they arise in practice."""
    print("Phase 5c: mixed shapes in a single term")
    global abstentions
    rng = random.Random(999)
    for p in (3, 4, 5):
        axes = list(range(p))
        domains = {a: DOM for a in axes}
        grids = {a: GRID for a in axes}
        done, total, runs, combos = 0, 0, 0, 0
        for mu in (LEBESGUE, MAGNITUDE):
            mm = {a: mu for a in axes}
            for _ in range(8):
                total += 1
                t = PairTerm(1.0)
                for i in range(p):
                    if rng.random() < 0.5:
                        t.weights[i] = Step((1.0,), (1.0, rng.choice([0.5, 2.0])))
                for i in range(p):
                    for j in range(i + 1, p):
                        if rng.random() < 0.7:
                            maker = SHAPES[rng.choice(list(SHAPES))]
                            t.preds.append(PairPredicate(
                                i, j, maker(rng),
                                upper=rng.random() < 0.7,
                                strict=rng.random() < 0.5))
                want = integrate_exact([t], axes, mm, domains, grids)
                c = Counters()
                got = eliminate_all([t], axes, mm, domains, grids, c,
                                    term_cap=400000)
                if got is None:
                    abstentions += 1
                    continue
                done += 1
                check(f"mixed p={p} {mu.name}", got, want)
                runs = max(runs, c.max_runs)
                combos = max(combos, c.max_combos)
        status = "ok" if done == total else "partial"
        print(f"    mixed   p={p}: {status}  completed {done}/{total}  "
              f"max runs={runs}  max combos={combos}")


if __name__ == "__main__":
    test_runs_unit()
    test_shaped()
    test_mixed()
    print()
    print(f"abstentions (outside implemented grammar): {abstentions}")
    if failures:
        print(f"{len(failures)} FAILURES")
        for n, g, w in failures[:8]:
            print(f"  {n}: {g} vs {w}")
        sys.exit(1)
    print("No mismatches.")

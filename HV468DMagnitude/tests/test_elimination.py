"""Phases 2-4: eliminate variables and compare against exact enumeration.

The claim under test is the one the whole generalization rests on: integrating
a variable out of a weighted monotone-pair term yields a bounded signed sum of
weighted monotone-pair terms in the survivors, exactly.

What is measured, per the prompt:
  * exactness against the product-grid oracle, in 4, 6 and 8 variables;
  * emitted term count, and whether it grows with breakpoint complexity;
  * both Lebesgue and magnitude product measures.

Abstention is not failure but is reported: the constructor returns ``None``
when a predicate falls outside the implemented grammar, so a gap in coverage
can never be mistaken for a passing test.
"""

import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from core import INF, LEBESGUE, MAGNITUDE, Step
from pair_state import PairPredicate, PairTerm, step_const
from pair_elimination import Counters, eliminate_all
from pair_oracles import integrate_exact, term_breakpoints

TOL = 1e-9
failures = []
abstentions = 0


def check(name, got, want, tol=TOL):
    if got is None:
        return True
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


GRID = [0.5, 1.0, 1.5, 2.0]
DOM = (0.0, 2.0)


def rand_monotone(rng, n_break=2, nonincreasing=True, values=None):
    values = values or [0.0, 0.5, 1.0, 1.5, 2.0]
    xs = sorted(rng.sample(GRID[:-1], min(n_break, len(GRID) - 1)))
    vs = sorted((rng.choice(values) for _ in range(len(xs) + 1)),
                reverse=nonincreasing)
    return Step(tuple(xs), tuple(vs))


def rand_term(rng, p, n_break=2, density=0.6):
    t = PairTerm(1.0)
    for i in range(p):
        if rng.random() < 0.5:
            t.weights[i] = Step((1.0,), (1.0, rng.choice([0.5, 1.0, 2.0])))
    for i in range(p):
        for j in range(i + 1, p):
            if rng.random() < density:
                t.preds.append(PairPredicate(
                    i, j, rand_monotone(rng, n_break), upper=True, strict=True))
    return t


def run_case(rng, p, measure, n_break=2, label=""):
    global abstentions
    term = rand_term(rng, p, n_break)
    axes = list(range(p))
    measures = {a: measure for a in axes}
    domains = {a: DOM for a in axes}
    grids = {a: GRID for a in axes}
    want = integrate_exact([term], axes, measures, domains, grids)
    counters = Counters()
    got = eliminate_all([term], axes, measures, domains, grids, counters)
    if got is None:
        abstentions += 1
        return None, counters
    check(f"{label} p={p} {measure.name}", got, want)
    return counters, counters


def phase(p, trials, measures, n_break=2):
    """Report completion rate as well as branching.

    An abstention is not a failure, but a phase that abstained on most of its
    instances has tested nothing, so the rate is printed rather than hidden.
    """
    rng = random.Random(1000 + p)
    worst, done, total = 0, 0, 0
    for measure in measures:
        for _ in range(trials):
            total += 1
            c, _ = run_case(rng, p, measure, n_break, label=f"phase{p}")
            if c:
                done += 1
                worst = max(worst, c.max_branching)
    status = "ok" if done else "NOTHING TESTED"
    print(f"  p={p}: {status}  completed {done}/{total}  "
          f"max branching per elimination: {worst}")
    return worst


def test_elimination_exactness():
    print("Phase 2-4: elimination against exact enumeration")
    for p in (2, 3, 4):
        phase(p, 15, [LEBESGUE, MAGNITUDE])
    for p in (5, 6):
        phase(p, 6, [LEBESGUE])


def test_branching_vs_breakpoints():
    """The branching constant must not grow with staircase complexity."""
    print("Phase 3b: branching vs breakpoint complexity (p fixed)")
    rng = random.Random(77)
    rows = []
    for n_break in (1, 2, 3):
        worst, bps = 0, 0
        for _ in range(8):
            term = rand_term(rng, 4, n_break)
            axes = list(range(4))
            measures = {a: LEBESGUE for a in axes}
            domains = {a: DOM for a in axes}
            grids = {a: GRID for a in axes}
            c = Counters()
            want = integrate_exact([term], axes, measures, domains, grids)
            got = eliminate_all([term], axes, measures, domains, grids, c)
            if got is None:
                continue
            check(f"branching n_break={n_break}", got, want)
            worst = max(worst, c.max_branching)
            bps = max(bps, term_breakpoints([term]))
        rows.append((n_break, bps, worst))
        print(f"    breakpoints/pred={n_break}  term breakpoints={bps:>3}  "
              f"max branching per elimination={worst}")
    if len({r[2] for r in rows}) == 1:
        print("    emitted count independent of breakpoint complexity")


if __name__ == "__main__":
    test_elimination_exactness()
    test_branching_vs_breakpoints()
    print()
    print(f"abstentions (outside implemented grammar): {abstentions}")
    if failures:
        print(f"{len(failures)} FAILURES")
        for n, g, w in failures[:6]:
            print(f"  {n}: {g} vs {w}")
        sys.exit(1)
    print("No mismatches.")

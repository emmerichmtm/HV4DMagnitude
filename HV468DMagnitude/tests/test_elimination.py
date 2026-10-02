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
    if got != got:                                  # NaN is never a pass
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got NaN want {want!r}")
        return False
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


GRID = [0.5, 1.0, 1.5, 2.0]
DOM = (0.0, 2.0)
#: The oracle enumerates the full product grid, so its cost is
#: (cells per axis)^p.  Beyond six variables the grid has to shrink or the
#: oracle, not the constructor, becomes the bottleneck.
COARSE = [1.0, 2.0]

#: The oracle evaluates one representative per cell, which is exact only if
#: every function in the term is constant on each cell.  So a generated
#: staircase must put *both* its breakpoints and its values on the same grid the
#: oracle cells are cut from -- a value off the grid becomes a region boundary
#: inside a cell, and the oracle silently returns the wrong number.  Mixing the
#: fine staircases with the coarse grid is exactly the artefact that made the
#: first p=8 run look like a constructor failure.
VALUES = {
    4: [0.0, 0.5, 1.0, 1.5, 2.0],
    2: [0.0, 1.0, 2.0],
}


def grid_for(p):
    return GRID if p <= 6 else COARSE


def values_for(p):
    return VALUES[4] if p <= 6 else VALUES[2]


def rand_monotone(rng, grid, values, n_break=2, nonincreasing=True):
    """A monotone staircase with breakpoints and values on ``grid``."""
    interior = [g for g in grid if g < max(grid)]
    xs = sorted(rng.sample(interior, min(n_break, len(interior))))
    vs = sorted((rng.choice(values) for _ in range(len(xs) + 1)),
                reverse=nonincreasing)
    return Step(tuple(xs), tuple(vs))


def rand_term(rng, p, n_break=2, density=0.6):
    grid, values = grid_for(p), values_for(p)
    t = PairTerm(1.0)
    for i in range(p):
        if rng.random() < 0.5:
            t.weights[i] = Step((1.0,), (1.0, rng.choice([0.5, 1.0, 2.0])))
    for i in range(p):
        for j in range(i + 1, p):
            if rng.random() < density:
                t.preds.append(PairPredicate(
                    i, j, rand_monotone(rng, grid, values, n_break),
                    upper=True, strict=True))
    return t


def run_case(rng, p, measure, n_break=2, label="", cap=None):
    global abstentions
    term = rand_term(rng, p, n_break)
    axes = list(range(p))
    grid = grid_for(p)
    measures = {a: measure for a in axes}
    domains = {a: DOM for a in axes}
    grids = {a: grid for a in axes}
    want = integrate_exact([term], axes, measures, domains, grids)
    counters = Counters()
    got = eliminate_all([term], axes, measures, domains, grids, counters,
                        term_cap=cap)
    if got is None:
        abstentions += 1
        return None, counters
    check(f"{label} p={p} {measure.name}", got, want)
    return counters, counters


def phase(p, trials, measures, n_break=2, cap=None):
    """Report completion rate as well as branching.

    An abstention is not a failure, but a phase that abstained on most of its
    instances has tested nothing, so the rate is printed rather than hidden.
    """
    rng = random.Random(1000 + p)
    worst, done, total = 0, 0, 0
    for measure in measures:
        for _ in range(trials):
            total += 1
            c, _ = run_case(rng, p, measure, n_break, label=f"phase{p}",
                            cap=cap)
            if c:
                done += 1
                worst = max(worst, c.max_branching)
    status = "ok" if done == total else ("partial" if done else "NOTHING TESTED")
    print(f"  p={p}: {status}  completed {done}/{total}  "
          f"max branching per elimination: {worst}")
    return worst


def test_elimination_exactness():
    print("Phase 2-4: elimination against exact enumeration")
    for p in (2, 3, 4):
        phase(p, 15, [LEBESGUE, MAGNITUDE])
    for p in (5, 6):
        phase(p, 6, [LEBESGUE, MAGNITUDE])
    for p in (8,):
        phase(p, 4, [LEBESGUE, MAGNITUDE], cap=200000)


def test_branching_vs_breakpoints():
    """Obligation 7: the branching constant must not grow with staircase
    complexity.

    The graph of predicates has to be held fixed while only the staircase
    resolution varies, or the comparison measures two things at once.  Each
    row therefore re-seeds the generator, so every row draws the *same*
    sequence of predicate graphs and differs only in ``n_break``.
    """
    print("Phase 3b: branching vs breakpoints (same graphs, p=4 fixed)")
    rows = []
    for n_break in (1, 2, 3):
        rng = random.Random(77)                 # identical graphs per row
        worst, bps, done = 0, 0, 0
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
            done += 1
            check(f"branching n_break={n_break}", got, want)
            worst = max(worst, c.max_branching)
            bps = max(bps, term_breakpoints([term]))
        rows.append((n_break, bps, worst))
        print(f"    breaks/pred={n_break}  term breakpoints={bps:>3}  "
              f"completed={done}/8  max branching={worst}")
    if len({r[2] for r in rows}) == 1:
        print("    branching independent of breakpoint complexity")
    else:
        lo, hi = rows[0][2], rows[-1][2]
        print(f"    branching moved {lo} -> {hi} as breakpoints "
              f"{rows[0][1]} -> {rows[-1][1]}")


def test_breakpoints_along_chain():
    """Obligation 6: breakpoints per surviving term must stay bounded.

    If each elimination doubled the breakpoints of the terms it emits, the
    per-term cost would grow with chain depth even at bounded branching, and
    the fixed-dimensional claim would be empty.
    """
    print("Phase 3c: breakpoints per term along the elimination chain")
    rng = random.Random(4242)
    for p in (6, 8):
        term = rand_term(rng, p, 2)
        axes = list(range(p))
        grid = grid_for(p)
        measures = {a: LEBESGUE for a in axes}
        domains = {a: DOM for a in axes}
        grids = {a: grid for a in axes}
        from pair_elimination import eliminate_variable
        state, trace = [term], []
        for x in axes:
            nxt = []
            for t in state:
                got = eliminate_variable(t, x, measures[x], domains[x],
                                         grids[x])
                if got is None:
                    nxt = None
                    break
                nxt.extend(got)
            if nxt is None:
                trace.append("abstain")
                break
            state = nxt
            if not state:
                break
            trace.append(round(term_breakpoints(state) / len(state), 2))
        print(f"    p={p}: breakpoints/term after each elimination: {trace}")


if __name__ == "__main__":
    test_elimination_exactness()
    test_branching_vs_breakpoints()
    test_breakpoints_along_chain()
    print()
    print(f"abstentions (outside implemented grammar): {abstentions}")
    if failures:
        print(f"{len(failures)} FAILURES")
        for n, g, w in failures[:6]:
            print(f"  {n}: {g} vs {w}")
        sys.exit(1)
    print("No mismatches.")

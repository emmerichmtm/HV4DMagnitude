"""Measure the growth of the pair-elimination constructor.

The report's complexity claim rests on a branching constant ``B_p`` that
depends on the number of variables ``p`` but *not* on the instance size ``N``
(the number of staircase breakpoints).  That separation is what makes the
method fixed-parameter in ``p`` and linear-ish in ``N``; if ``B_p`` grew with
``N`` the stated bound would be wrong.

Two sweeps, deliberately kept apart:

``size``
    ``p`` and the predicate graph are held fixed while only the number of
    breakpoints per staircase grows.  Each row re-seeds the generator, so every
    row draws the *same* sequence of graphs and differs in resolution alone.
    Without that, the sweep measures graph density and resolution at once and
    tells you nothing about either.

``dimension``
    breakpoints held fixed while ``p`` grows.

Exactness is not checked here.  The oracle needs its cells aligned to the
staircases, which caps the resolution far below what is needed to see a trend
in ``N``; exactness lives in ``tests/test_elimination.py`` and the two
measurements are kept separate so neither constrains the other.
"""

from __future__ import annotations

import os
import random
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from core import LEBESGUE, MAGNITUDE, Step
from pair_elimination import Counters, eliminate_all
from pair_oracles import term_breakpoints
from pair_state import PairPredicate, PairTerm

DOM = (0.0, 2.0)
TERM_CAP = 300000


def ladder(k: int):
    """``k`` interior grid points on ``(0, 2]``, as both breaks and values."""
    return [2.0 * i / k for i in range(1, k + 1)]


def rand_monotone(rng, grid, n_break, nonincreasing=True):
    interior = [g for g in grid if g < max(grid)]
    xs = sorted(rng.sample(interior, min(n_break, len(interior))))
    vals = [0.0] + grid
    vs = sorted((rng.choice(vals) for _ in range(len(xs) + 1)),
                reverse=nonincreasing)
    return Step(tuple(xs), tuple(vs))


def rand_term(rng, p, grid, n_break, density=0.6):
    t = PairTerm(1.0)
    for i in range(p):
        if rng.random() < 0.5:
            t.weights[i] = Step((1.0,), (1.0, rng.choice([0.5, 1.0, 2.0])))
    for i in range(p):
        for j in range(i + 1, p):
            if rng.random() < density:
                t.preds.append(PairPredicate(
                    i, j, rand_monotone(rng, grid, n_break),
                    upper=True, strict=True))
    return t


def measure(p, n_break, k, trials, seed, measure_obj=LEBESGUE):
    """Return (breakpoints, max branching, mean branching, peak state, done)."""
    rng = random.Random(seed)
    grid = ladder(k)
    axes = list(range(p))
    measures = {a: measure_obj for a in axes}
    domains = {a: DOM for a in axes}
    grids = {a: grid for a in axes}
    brs, peaks, bps, done, abst = [], [], [], 0, 0
    for _ in range(trials):
        term = rand_term(rng, p, grid, n_break)
        c = Counters()
        got = eliminate_all([term], axes, measures, domains, grids, c,
                            term_cap=TERM_CAP)
        bps.append(term_breakpoints([term]))
        if got is None:
            abst += 1
            continue
        done += 1
        brs.append(c.max_branching)
        peaks.append(c.peak_state)
    if not brs:
        return (max(bps) if bps else 0, 0, 0.0, 0, done, abst)
    return (max(bps), max(brs), statistics.mean(brs), max(peaks), done, abst)


def sweep_size(p, ks, trials=10):
    print(f"\n  instance size at p={p} (graph held fixed, resolution varies)")
    print("    breaks/pred  term breakpoints  max B  mean B  peak state  done")
    for k in ks:
        n_break = k - 1
        bps, mx, mean, peak, done, abst = measure(
            p, n_break, k, trials, seed=90210)
        print(f"    {n_break:>11}  {bps:>16}  {mx:>5}  {mean:>6.1f}  "
              f"{peak:>10}  {done:>2}/{trials}"
              + (f"  ({abst} abstained)" if abst else ""))


def sweep_dimension(ps, k=4, trials=10):
    print(f"\n  dimension at fixed resolution (k={k} grid points)")
    print("    p   term breakpoints  max B  mean B  peak state  done")
    for p in ps:
        bps, mx, mean, peak, done, abst = measure(
            p, 2, k, trials, seed=31337)
        print(f"    {p:>2}  {bps:>16}  {mx:>5}  {mean:>6.1f}  "
              f"{peak:>10}  {done:>2}/{trials}"
              + (f"  ({abst} abstained)" if abst else ""))


if __name__ == "__main__":
    print("Pair-elimination growth probe")
    print("  B = terms emitted by one elimination from one incoming term")
    sweep_size(4, [2, 3, 4, 6, 8, 12, 16])
    sweep_size(6, [2, 3, 4, 6, 8])
    sweep_dimension([2, 3, 4, 5, 6, 7, 8])

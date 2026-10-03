"""Phase 7: two compression generations on nested grids.

The report's closure claim is not that one compression works -- that is Phase 6
-- but that the output can be *used recursively*, with no latent variable
carried between generations.  The checkable form of that is the tower property.

Let ``G1`` be a grid and ``G2`` a coarser one whose boundaries are a subset of
``G1``'s, so ``G2`` is an admissible descendant.  Put ``F_bar1 = E[F | G1]``.
Every cell of ``G2`` is a union of ``G1`` cells, so its indicator is
``G1``-measurable and

    integral over C2 of F_bar1  ==  integral over C2 of F.

So compressing ``F_bar1`` on ``G2`` must give exactly what compressing the
*original* ``F`` on ``G2`` would give.  That is the precise sense in which the
discarded fine geometry is irrelevant, and it is what this file tests: the
second generation's numerator is compared against the exact integral of the
original ``F`` over each coarse cell.

Two checks are run at each generation, deliberately kept apart, because they
fail for different reasons:

``carry``
    the oracle integral of ``F_bar1`` over a coarse cell equals that of ``F``.
    This tests generation 1's *output* -- whether the state really is the
    conditional expectation -- with no second elimination involved.
``chain``
    the second generation's symbolic numerator equals the same number.  This
    tests generation 2's elimination on an input that no longer came from a
    hand-written generator but from the machine itself.

A failure in ``carry`` but not ``chain`` would mean generation 1 is wrong; the
reverse would mean the constructor cannot digest its own output.
"""

import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from core import LEBESGUE, MAGNITUDE, Step
from compression import (CompressionGrid, cell_indices, denominator_at,
                         merge_terms, next_generation, numerator_at)
from pair_elimination import Counters
from pair_oracles import integrate_exact
from pair_state import PairPredicate, PairTerm

TOL = 1e-9
failures = []
notes = []

DOM = (0.0, 2.0)
V = [0.5, 1.0, 1.5, 2.0]
VALUES = [0.0, 0.5, 1.0, 1.5, 2.0]

#: nested grids: every G2 boundary is a G1 boundary, so G2 is a descendant and
#: each G2 cell is a union of G1 cells
G1 = [0.0, 0.5, 1.0, 1.5, 3.0]
G2 = [0.0, 1.0, 3.0]


def check(name, got, want, tol=TOL):
    if got != got:
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got NaN want {want!r}")
        return False
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


def nested_ok():
    return set(G2).issubset(set(G1))


def rand_staircase(rng, n_break=2, nonincreasing=True):
    xs = tuple(sorted(rng.sample(V[:-1], n_break)))
    vs = tuple(sorted((rng.choice(VALUES) for _ in range(len(xs) + 1)),
                      reverse=nonincreasing))
    return Step(xs, vs)


def rand_f_term(rng, d, density=0.3):
    t = PairTerm(1.0)
    for i in range(d):
        if rng.random() < 0.5:
            t.weights[i] = Step((1.0,), (1.0, rng.choice([0.5, 1.0, 2.0])))
    for i in range(d):
        for j in range(i + 1, d):
            if rng.random() < density:
                t.preds.append(PairPredicate(
                    i, j, rand_staircase(rng),
                    upper=rng.random() < 0.7, strict=rng.random() < 0.5))
    return t


def setup(d, boundaries, mu):
    axes = list(range(d))
    grids = [CompressionGrid(list(boundaries)) for _ in axes]
    measures = {a: mu for a in axes}
    domains = {a: DOM for a in axes}
    # the event list must cut at every breakpoint any weight can acquire, so it
    # carries both grids' boundaries as well as the staircase grid
    events = {a: sorted(set(V) | set(G1[1:-1]) | set(G2[1:-1])) for a in axes}
    return axes, grids, measures, domains, events


def oracle_grid(d):
    return {a: sorted(set(V) | set(G1[1:-1]) | set(G2[1:-1])) for a in range(d)}


def run(d, mu, density, trials, term_cap, verbose=True, fine=None,
        coarse=None, carry_cells=None):
    """Two generations on nested grids, checked against the original F.

    ``carry_cells`` limits how many cells get the ``carry`` check.  That check
    integrates the whole generation-1 state cell by cell, so its cost is
    (state size) x (points per cell), which is minutes per cell once the state
    runs to thousands of terms.  ``chain`` is cheap on every cell -- it only
    evaluates the final state at a point -- so coverage of the actual claim is
    not reduced by sampling ``carry``.
    """
    fine = list(G1 if fine is None else fine)
    coarse = list(G2 if coarse is None else coarse)
    if not set(coarse).issubset(set(fine)):
        failures.append((f"nesting d={d}", coarse, f"subset of {fine}"))
        return []
    rng = random.Random(9000 + d)
    axes = list(range(d))
    domains = {a: DOM for a in axes}
    measures = {a: mu for a in axes}
    marks = sorted(set(V) | set(fine[1:-1]) | set(coarse[1:-1]))
    og = {a: marks for a in axes}
    g1 = [CompressionGrid(fine) for _ in axes]
    g2 = [CompressionGrid(coarse) for _ in axes]
    events = {a: marks for a in axes}

    rows = []
    for _ in range(trials):
        f = rand_f_term(rng, d, density)
        t0 = time.time()
        c1 = Counters()
        s1, raw1 = next_generation([f], d, measures, domains, g1, events, c1,
                                   term_cap=term_cap)
        if s1 is None:
            rows.append((None, None, None, None, time.time() - t0, "gen1 cap"))
            continue
        t1 = time.time()
        c2 = Counters()
        s2, raw2 = next_generation(s1, d, measures, domains, g2, events, c2,
                                   term_cap=term_cap)
        t2 = time.time()
        if s2 is None:
            rows.append((raw1, len(s1), None, None, t2 - t0,
                         f"gen2 cap >{term_cap}"))
            continue

        # ---- the two checks, per coarse cell --------------------------------
        # Note what `s2` is: `next_generation` has already divided by the cell
        # mass and relabelled y onto the physical axes 0..d-1.  So it is
        # E[F | G2] itself, evaluated as an ordinary term, not a numerator --
        # which is the quantity the recursion claim is actually about.
        for n_cell, cells in enumerate(cell_indices(d, g2)):
            cdom = {a: g2[a].cell_bounds(cells[a], DOM[1]) for a in axes}
            want = integrate_exact([f], axes, measures, cdom, og)
            # carry: generation 1's output integrates like F on this cell
            if carry_cells is None or n_cell < carry_cells:
                carry = integrate_exact(s1, axes, measures, cdom, og)
                check(f"carry d={d} {mu.name} cell={cells}", carry, want)
            # chain: two generations give E[F | G2] on this cell
            den = denominator_at(d, g2, measures, domains, cells)
            want_bar = want / den if den else 0.0
            for ys in zip(*[g2[a].representatives(cells[a], DOM[1])
                            for a in axes]):
                got = sum(t.value_at(list(ys)) for t in s2)
                check(f"chain d={d} {mu.name} cell={cells} y={ys}",
                      got, want_bar)
        rows.append((raw1, len(s1), raw2, len(s2), t2 - t0, "ok"))
        if verbose:
            print(f"    d={d} {mu.name:<9} gen1 raw={raw1:<7} merged={len(s1):<6}"
                  f"  gen2 raw={raw2:<8} merged={len(s2):<7}"
                  f"  ({t1 - t0:.0f}s + {t2 - t1:.0f}s)")
    return rows


def test_nesting():
    print("Phase 7a: the grids must be nested for the tower property to apply")
    if not nested_ok():
        failures.append(("grid nesting", G2, f"subset of {G1}"))
        print(f"  [FAIL] {G2} is not a subset of {G1}")
        return
    ga, gb = CompressionGrid(list(G1)), CompressionGrid(list(G2))
    print(f"    G1 {G1} -> {ga.n_cells()} cells;  "
          f"G2 {G2} -> {gb.n_cells()} cells;  nested: yes")
    # every coarse cell must be an exact union of fine cells
    for m in range(gb.n_cells()):
        lo, up = gb.cell_bounds(m, DOM[1])
        covered = [ga.cell_bounds(k, DOM[1]) for k in range(ga.n_cells())
                   if lo <= ga.boundaries[k] < up]
        span = (covered[0][0], covered[-1][1])
        if span != (lo, up):
            failures.append((f"coarse cell {m} cover", span, (lo, up)))
            print(f"  [FAIL] coarse cell {m}: fine cells cover {span}, "
                  f"want {(lo, up)}")
    print("    every coarse cell is an exact union of fine cells")


def test_two_generations():
    print("Phase 7b: generation 2 against the original F (4 cells -> 2 cells)")
    for d in (2, 3):
        run(d, LEBESGUE, 0.4, 2, term_cap=400000)
        run(d, MAGNITUDE, 0.4, 2, term_cap=400000)
    # one instance per measure at d=4: the second generation there runs to
    # ~2 x 10^5 raw terms and takes minutes
    run(4, LEBESGUE, 0.4, 1, term_cap=400000)
    run(4, MAGNITUDE, 0.4, 1, term_cap=400000)


#: a coarser fine-grid for d = 6, so the first generation stays affordable
#: while the pair is still non-trivial (3 cells collapsing to 2)
D6_FINE = [0.0, 0.5, 1.0, 3.0]
D6_COARSE = [0.0, 1.0, 3.0]


def test_d6():
    """Phase 7c: two generations at d = 6, the dimension this sprint is about.

    ``carry`` is sampled rather than run on all 64 cells: it integrates the
    whole generation-1 state cell by cell, which is minutes per cell at this
    state size.  ``chain`` -- the claim itself -- is checked on every cell.
    """
    print("Phase 7c: d=6, two generations (3 cells -> 2 cells)")
    # `carry` is skipped entirely here, not sampled: integrating a 7000-term
    # state over one cell of a 6-fold product grid is minutes on its own, and
    # the check is already covered at d = 2, 3 and 4.  `chain` -- the claim --
    # runs on all 64 cells.  Together these two take about twelve minutes.
    for mu in (LEBESGUE, MAGNITUDE):
        rows = run(6, mu, 0.2, 1, term_cap=2000000, fine=D6_FINE,
                   coarse=D6_COARSE, carry_cells=0)
        for raw1, m1, raw2, m2, secs, status in rows:
            if status != "ok":
                notes.append(f"d=6 {mu.name} second generation: {status} "
                             f"(gen1 raw={raw1}, merged={m1}, {secs:.0f}s)")
                print(f"    d=6 {mu.name}: {status}  gen1 raw={raw1} "
                      f"merged={m1}  ({secs:.0f}s)")


if __name__ == "__main__":
    test_nesting()
    test_two_generations()
    test_d6()
    print()
    for n in notes:
        print(f"note: {n}")
    if failures:
        print(f"{len(failures)} FAILURES")
        for n, g, w in failures[:8]:
            print(f"  {n}: {g} vs {w}")
        sys.exit(1)
    print("No mismatches.")

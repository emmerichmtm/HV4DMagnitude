"""Phase 6: compression in ``2d`` variables, for d = 4 and d = 6.

This is the construction the general theorem actually needs, and the one that
distinguishes "validated at d = 8" from "validated at d = 4": the static
elimination chain works in ``p = d`` variables, but compression builds a term in
``p = 2d``.

The property under test is the definition of the conditional expectation, not a
weaker proxy.  For every cell ``C`` of the product grid,

    N_F(y) == integral of F over C,     for y in C,

where the left side is produced by building the ``2d``-variable integrand and
eliminating the ``d`` old variables symbolically, and the right side by exact
enumeration of ``F`` restricted to ``C``.  If those agree on every cell then
``F_bar = N_F / mass(C)`` is the conditional expectation, and the tower-property
argument of the report applies to it.

Each cell is probed at its floor *and* its midpoint.  The floor of the first
cell is the anchor ``0``, which is the single point where the cell convention
can be wrong without Lebesgue noticing.

An aggregate check is included as well: summing the numerators over all cells
must reproduce the integral of ``F`` over the whole ground set.  It is implied by
the per-cell checks, but it fails loudly if cells are double-counted or missed,
which a per-cell check on its own would not catch.
"""

import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from core import LEBESGUE, MAGNITUDE, Step
from compression import (CompressionGrid, cell_indices, compress,
                         denominator_at, divide_by_cell_mass, merge_terms,
                         numerator_at, reindex_to_physical)
from pair_elimination import Counters
from pair_oracles import integrate_exact
from pair_state import PairPredicate, PairTerm

TOL = 1e-9
failures = []
abstentions = 0

DOM = (0.0, 2.0)
#: staircase breakpoints and values, kept on one grid so the oracle's cells are
#: aligned with every function in the term (see test_elimination for why)
V = [0.5, 1.0, 1.5, 2.0]
VALUES = [0.0, 0.5, 1.0, 1.5, 2.0]
#: cell boundaries; the last one is deliberately above the ceiling so every cell
#: is [g_m, g_{m+1}) with the same decoration
CELLS = [0.0, 1.0, 3.0]


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


def rand_staircase(rng, n_break=2, nonincreasing=True):
    xs = tuple(sorted(rng.sample(V[:-1], n_break)))
    vs = tuple(sorted((rng.choice(VALUES) for _ in range(len(xs) + 1)),
                      reverse=nonincreasing))
    return Step(xs, vs)


def rand_f_term(rng, d, density=0.3):
    """A weighted monotone-pair term in the ``d`` physical variables."""
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


def run_instance(rng, d, mu, density, label, term_cap=400000):
    """One F: compress it, then check every cell against the oracle."""
    global abstentions
    axes = list(range(d))
    grids = [CompressionGrid(list(CELLS)) for _ in axes]
    measures = {a: mu for a in axes}
    domains = {a: DOM for a in axes}
    events = {a: sorted(set(V) | set(grids[a].interior)) for a in axes}
    oracle_grid = {a: V for a in axes}

    f_term = rand_f_term(rng, d, density)
    counters = Counters()
    terms = compress(f_term, d, measures, domains, grids, events, counters,
                     term_cap=term_cap)
    if terms is None:
        abstentions += 1
        return None, 0

    # The conditional expectation, as terms rather than as a quotient of
    # numbers.  If the report's "only changes unary step densities" is right,
    # these evaluate to numerator/denominator everywhere.
    bar = divide_by_cell_mass(terms, d, grids, measures, domains)
    # And the closure claim: the output must mention no x variable at all, so
    # relabelling y back to the physical axes is a pure rename.
    try:
        physical = reindex_to_physical(bar, d)
    except ValueError as exc:
        failures.append((f"{label} d={d} {mu.name} closure", str(exc), "clean"))
        print(f"  [FAIL] {label} d={d} {mu.name} closure: {exc}")
        physical = None

    total_num = 0.0
    for cells in cell_indices(d, grids):
        cell_dom = {a: grids[a].cell_bounds(cells[a], DOM[1]) for a in axes}
        want = integrate_exact([f_term], axes, measures, cell_dom, oracle_grid)
        den = denominator_at(d, grids, measures, domains, cells)
        for ys in zip(*[grids[a].representatives(cells[a], DOM[1])
                        for a in axes]):
            got = numerator_at(terms, d, ys)
            check(f"{label} d={d} {mu.name} cell={cells} y={ys}", got, want)
            # F_bar as a term must equal the quotient
            got_bar = numerator_at(bar, d, ys)
            check(f"{label} d={d} {mu.name} Fbar cell={cells} y={ys}",
                  got_bar, want / den if den else 0.0)
            if physical is not None:
                got_phys = sum(t.value_at(list(ys)) for t in physical)
                check(f"{label} d={d} {mu.name} reindexed cell={cells}",
                      got_phys, got_bar)
        # the aggregate uses one representative per cell
        floors = [grids[a].cell_bounds(cells[a], DOM[1])[0] for a in axes]
        total_num += numerator_at(terms, d, floors)

    whole = integrate_exact([f_term], axes, measures, domains, oracle_grid)
    check(f"{label} d={d} {mu.name} SUM over cells", total_num, whole)
    return counters, len(terms)


def phase(d, trials, density, measures, term_cap=400000):
    """``final`` is the report's C_d: output terms from ONE incoming term."""
    rng = random.Random(4000 + d)
    done, total, peak, runs, comb, final = 0, 0, 0, 0, 0, 0
    t0 = time.time()
    for mu in measures:
        for _ in range(trials):
            total += 1
            c, n_out = run_instance(rng, d, mu, density, "phase6", term_cap)
            if c is None:
                continue
            done += 1
            peak = max(peak, c.peak_state)
            runs = max(runs, c.max_runs)
            comb = max(comb, c.max_combos)
            final = max(final, n_out)
    status = "ok" if done == total else ("partial" if done else "NOTHING TESTED")
    print(f"  d={d} (p={2 * d}): {status}  completed {done}/{total}  "
          f"output terms (C_d)={final}  peak state={peak}  "
          f"max runs={runs}  max combos={comb}  ({time.time() - t0:.0f}s)")


def test_grid_convention():
    """The cell maps must be right-continuous and partition the ground set."""
    print("Phase 6a: cell convention")
    g = CompressionGrid(list(CELLS))
    pred, succ = g.pred_step(), g.succ_step()
    for y, want_lo, want_up in ((0.0, 0.0, 1.0), (0.5, 0.0, 1.0),
                                (1.0, 1.0, 3.0), (2.0, 1.0, 3.0)):
        ok = pred(y) == want_lo and succ(y) == want_up
        if not ok:
            failures.append((f"cell map at y={y}",
                             (pred(y), succ(y)), (want_lo, want_up)))
            print(f"  [FAIL] cell map at y={y}: "
                  f"({pred(y)}, {succ(y)}) want ({want_lo}, {want_up})")
    # masses must sum to the measure of the whole axis
    for mu in (LEBESGUE, MAGNITUDE):
        got = sum(g.cell_mass(m, mu, DOM[1]) for m in range(g.n_cells()))
        want = mu.atom0 + mu.density * (DOM[1] - DOM[0])
        check(f"cell masses sum ({mu.name})", got, want)
    print(f"    boundaries {CELLS} -> {g.n_cells()} cells per axis, "
          f"masses sum correctly under both measures")


def test_compression():
    print("Phase 6b: N_F(y) against the exact integral over each cell")
    phase(4, 4, 0.3, [LEBESGUE, MAGNITUDE])
    phase(6, 2, 0.3, [LEBESGUE, MAGNITUDE])


def test_merge():
    """Phase 6c: how far does canonical merging collapse C_d?

    `\\ref{prop:complexity}` concedes that nothing proves a signed sum of
    distinct primitives can be merged back to a constant number.  This measures
    what the only geometry-free merge actually recovers -- adding the
    coefficients of terms that are identical, and dropping what cancels -- and
    checks that the merged state still represents the same function.

    Reported, not asserted: it is a cost measurement, and the only pass/fail
    here is that merging preserves the value on every cell.
    """
    print("Phase 6c: canonical merging of the output state")
    print("    d  measure     raw C_d   merged  ratio   exact after merge")
    for d in (4, 6):
        axes = list(range(d))
        grids = [CompressionGrid(list(CELLS)) for _ in axes]
        domains = {a: DOM for a in axes}
        events = {a: sorted(set(V) | set(grids[a].interior)) for a in axes}
        oracle_grid = {a: V for a in axes}
        rng = random.Random(4000 + d)
        for mu in (LEBESGUE, MAGNITUDE):
            measures = {a: mu for a in axes}
            f_term = rand_f_term(rng, d, 0.3)
            terms = compress(f_term, d, measures, domains, grids, events,
                             Counters(), term_cap=400000)
            if terms is None:
                print(f"    {d:>2}  {mu.name:<10} abstained")
                continue
            merged = merge_terms(terms)
            bad = 0
            for cells in cell_indices(d, grids):
                cell_dom = {a: grids[a].cell_bounds(cells[a], DOM[1])
                            for a in axes}
                want = integrate_exact([f_term], axes, measures, cell_dom,
                                       oracle_grid)
                ys = [grids[a].cell_bounds(cells[a], DOM[1])[0] for a in axes]
                if abs(numerator_at(merged, d, ys) - want) > TOL:
                    bad += 1
            if bad:
                failures.append((f"merge d={d} {mu.name}", f"{bad} cells", "0"))
                print(f"  [FAIL] merge d={d} {mu.name}: {bad} cells changed")
            ratio = len(terms) / max(len(merged), 1)
            print(f"    {d:>2}  {mu.name:<10} {len(terms):>7}  "
                  f"{len(merged):>7}  {ratio:>4.1f}x   "
                  f"{'yes' if not bad else 'NO'}")
    print("    (an empty row means that F integrated to zero everywhere)")


if __name__ == "__main__":
    test_grid_convention()
    test_compression()
    test_merge()
    print()
    print(f"abstentions (outside implemented grammar): {abstentions}")
    if failures:
        print(f"{len(failures)} FAILURES")
        for n, g, w in failures[:8]:
            print(f"  {n}: {g} vs {w}")
        sys.exit(1)
    print("No mismatches.")

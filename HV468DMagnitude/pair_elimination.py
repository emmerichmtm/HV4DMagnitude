"""Variable elimination for weighted monotone-pair terms.

This implements the constructor the generalized report's Lemma 2 (one-variable
elimination) describes, and exposes the pair version as two applications of it.

A structural observation, recorded in the audit
-----------------------------------------------
The report proves Lemma 2 (one variable) and Lemma 3 (a pair) separately.  For
*closure* purposes the pair lemma is not needed: the mutual ``x``--``y``
staircase is, from ``x``'s point of view, just one more bound candidate that
happens to depend on ``y``, so eliminating ``x`` and then ``y`` by Lemma 2
twice already yields the claim.  The pair version buys a better constant and a
cleaner closed form (the two-branch staircase formula), not extra generality.
We therefore test Lemma 2 as the load-bearing statement and obtain the pair
result by iteration.

What one elimination does
-------------------------
Every predicate incident to the eliminated variable ``x`` becomes, after
decorated inversion, a *finite union of decorated intervals* of feasible ``x``,
whose endpoints are unary functions of exactly one surviving variable.  With
``O_p(1)`` candidates per side, the surviving space is partitioned by which
candidate is active; on each label the integral is ``Q(U) - Q(L^-)`` for the
unary cumulative ``Q`` of ``x``'s own density, and every label condition
compares two unary maps, which is again a pair predicate (or a unary filter when
both sides involve the same variable).

Inversion is a union of runs, not a single interval
---------------------------------------------------
Inverting ``xi_j <| f(x)`` into a constraint on ``x`` does not in general give a
single-sided bound, and does not even give a single interval:

* when ``f`` is monotone the feasible ``x``-set is a prefix or a suffix, but the
  *degenerate* pieces (all ``x`` feasible, no ``x`` feasible) belong to the
  opposite family, so forcing one side misrepresents them;
* when ``f`` has one turn -- which this class provably contains, see
  ``pair_state`` -- the superlevel set of a *hill* is a middle interval with two
  nontrivial endpoints, while that of a *valley* is the union of two disjoint
  runs.

The integral over a disjoint union is the sum of the integrals, so a feasible
set with ``r`` runs is handled by splitting the term into ``r`` signed terms,
the ``k``-th carrying the ``k``-th run's endpoints.  Two things make this less
obvious than it sounds.  First, the number of runs varies *per piece* of the
surviving variable, so the split is by run *index*: pieces with fewer runs pad
their missing indices with the degenerate ``[0, 0)``, which carries exactly zero
mass under either measure.  Second, several predicates incident to ``x``
intersect, and the intersection of two unions of runs is the union over pairs of
runs -- so the branches multiply across incident predicates.  The product is
still bounded by a function of ``p`` alone, which is what the lemma needs, and
``Counters.max_runs`` records the worst run count actually seen so the claim can
be checked rather than assumed.

Abstention must be total
------------------------
Every "cannot express this" path returns ``None`` all the way out.  An earlier
version let one of them -- a comparison against a non-monotone bound -- fall
through as *unconstrained*, which silently deletes a constraint and inflates the
integral.  Dropping a constraint is not a conservative approximation; it is a
wrong answer that looks like a passing test until the term class is wide enough
to reach it.  ``_impose`` therefore has three outcomes: a (possibly split) list
of pieces, the empty list for "identically false", and ``None`` for abstain.
With run splitting in place no path abstains any more, but the distinction is
kept: it is what makes a future gap visible instead of silent.
"""

from __future__ import annotations

import itertools
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core import INF, Measure1D, Step
from identities import _cum
from pair_state import (DecoratedBound, PairPredicate, PairTerm, step_const,
                        merge_same_orientation)

__all__ = ["Candidate", "eliminate_variable", "eliminate_all", "Counters"]


class Counters:
    """Instrumentation required by the prompt.

    ``max_branching`` is the quantity the lemma actually bounds: the largest
    number of terms produced from *one* incoming term by *one* elimination.
    A running total would conflate that with the depth of the chain.

    ``max_runs`` and ``max_combos`` exist to check the run-splitting story.
    ``max_runs`` is the largest number of disjoint feasible runs seen in a
    single inversion, and ``max_combos`` the largest product of run counts over
    the predicates incident to one eliminated variable.  Both must stay bounded
    in the instance size or the construction is not what the lemma claims.
    """

    def __init__(self):
        self.labels = 0
        self.terms_in = 0
        self.terms_out = 0
        self.eliminations = 0
        self.max_branching = 0
        self.max_labels_per_elim = 0
        self.peak_state = 0
        self.max_runs = 0
        self.max_combos = 0

    def record(self, produced: int, labels: int) -> None:
        self.terms_in += 1
        self.terms_out += produced
        self.max_branching = max(self.max_branching, produced)
        self.max_labels_per_elim = max(self.max_labels_per_elim, labels)

    def note_runs(self, runs: int) -> None:
        self.max_runs = max(self.max_runs, runs)

    def note_combos(self, combos: int) -> None:
        self.max_combos = max(self.max_combos, combos)

    def __repr__(self):
        return (f"Counters(elims={self.eliminations}, "
                f"max_branching={self.max_branching}, "
                f"max_labels={self.max_labels_per_elim}, "
                f"max_runs={self.max_runs}, max_combos={self.max_combos})")


# --------------------------------------------------------------------------- #
# Candidate bounds
# --------------------------------------------------------------------------- #


class Candidate:
    """A bound on the eliminated variable, as a unary map of one survivor.

    ``var is None`` marks a constant bound (the ground set's own limits).
    ``upper`` says which side; ``closed`` whether the endpoint is attained.
    """

    __slots__ = ("var", "f", "upper", "closed", "tag")

    def __init__(self, var: Optional[int], f: Step, upper: bool, closed: bool,
                 tag: str = ""):
        self.var, self.f, self.upper, self.closed, self.tag = (
            var, f, upper, closed, tag)

    def at(self, point: Sequence[float]) -> float:
        return self.f(0.0) if self.var is None else self.f(point[self.var])

    def is_vacuous(self) -> bool:
        """True when this bound is implied by the ground set and can be dropped."""
        vs = self.f.vs
        if self.upper:
            return all(v == INF for v in vs)
        return all(v == 0.0 for v in vs) and self.closed

    def __repr__(self):
        side = "U" if self.upper else "L"
        return f"<{side} {self.tag} var={self.var}>"


# --------------------------------------------------------------------------- #
# Decorated inversion
# --------------------------------------------------------------------------- #


def _feasible_runs(f: Step, v: float, upper: bool,
                   strict: bool) -> List[Tuple[float, bool, float, bool]]:
    """Maximal feasible runs of ``{x : v <| f(x)}``, left to right.

    Each run is ``(L, L_closed, U, U_closed)``; the list is empty when nothing
    is feasible.  ``Step`` is right-continuous: piece ``k >= 1`` of ``f`` spans
    ``[xs[k-1], xs[k])`` and piece 0 spans ``[0, xs[0])``.  Hence a run starting
    at piece ``k`` is closed at ``xs[k-1]`` and one ending at piece ``k`` is open
    at ``xs[k]``.  The decorations are therefore uniform by construction --
    lower closed, upper open -- which is what lets a missing run be padded with
    ``[0, 0)`` and carry exactly zero mass under either measure.

    A monotone ``f`` yields one run, a hill one, a valley two; nothing here
    assumes any of those, so a boundary that has accumulated several turns is
    handled on the same footing.
    """
    if upper:
        ok = (lambda t: t > v) if strict else (lambda t: t >= v)
    else:
        ok = (lambda t: t < v) if strict else (lambda t: t <= v)
    mask = [ok(val) for val in f.vs]
    runs: List[Tuple[float, bool, float, bool]] = []
    n = len(mask)
    k = 0
    while k < n:
        if not mask[k]:
            k += 1
            continue
        start = k
        while k < n and mask[k]:
            k += 1
        last = k - 1
        L = 0.0 if start == 0 else f.xs[start - 1]
        U = INF if last == n - 1 else f.xs[last]
        runs.append((L, True, U, False))
    return runs


def _cum_atom(weight, atom_val: float, grid: Sequence[float], mu: Measure1D,
              t: float, closed: bool) -> float:
    """``_cum`` with the anchor's weight supplied separately.

    Identical to ``identities._cum`` except that the ``delta_0`` component is
    scaled by ``atom_val`` rather than by ``weight(0)``.  The two differ exactly
    when a term owns or disowns the anchor point -- see ``PairTerm`` -- and
    under Lebesgue (``atom0 = 0``) they never differ at all.
    """
    if t <= 0.0:
        return mu.atom0 * atom_val if (closed and t == 0.0) else 0.0
    total = mu.atom0 * atom_val
    pts = sorted({0.0} | {g for g in grid if 0.0 < g < t} | {t})
    for lo, hi in zip(pts, pts[1:]):
        total += mu.density * (hi - lo) * weight((lo + hi) / 2.0)
    return total


def _samples_for(breaks: Sequence[float],
                 interior: bool = False) -> List[float]:
    """One representative per piece of a step function with these breakpoints.

    ``interior=True`` picks a point strictly inside the first piece instead of
    its left endpoint ``0``.  That matters when the sampled relation is
    left-continuous: then the value *at* ``0`` differs from the value on
    ``(0, breaks[0])``, and sampling at ``0`` would spread the single-point
    value across the whole piece -- an error of full measure, not a
    measure-zero one.

    The degenerate case needs the same treatment and is easy to miss: with no
    breakpoints at all the single piece is ``[0, inf)``, whose interior is
    ``(0, inf)``.  Any positive representative does, because the function is
    constant there; returning ``0`` instead would describe the anchor and claim
    it holds everywhere.  That is exactly how a constant bound ``f == 0`` --
    the constraint ``x_j <= 0`` -- came to be read as "every ``x`` is feasible"
    for every ``x_j``, inventing mass on cells where the term is identically
    zero.
    """
    if not breaks:
        return [1.0 if interior else 0.0]
    out = [0.5 * breaks[0] if interior else 0.0]      # the piece [0, breaks[0])
    for a, b in zip(breaks, breaks[1:]):
        out.append(0.5 * (a + b))
    out.append(breaks[-1] + 1.0)
    return out


def _right_continuous(upper: bool, strict: bool) -> bool:
    """Is the inverted feasible set right-continuous in the surviving value?

    Inverting ``v <| f(x)`` gives a feasible ``x``-set that moves as ``v``
    crosses the values ``f`` attains.  For ``v < f(x)`` and ``v >= f(x)`` the
    set is constant on ``[a, b)`` -- right-continuous, matching ``Step``.  For
    ``v <= f(x)`` and ``v > f(x)`` it is constant on ``(a, b]`` instead, and a
    right-continuous step can then be correct everywhere except at the jumps.
    """
    return upper == strict


#: One branch of an inverted predicate: the bounds it contributes, plus weight
#: modifications to apply, as ``(axis, atom_multiplier, continuous_step)``.
#: A separate atom multiplier is what lets a branch own the anchor point
#: exclusively, or disown it, which a single step cannot express.
Branch = Tuple[List[Candidate], List[Tuple[int, float, Optional[Step]]]]


def _invert_to_branches(p: PairPredicate, x: int,
                        counters: Optional[Counters] = None
                        ) -> Optional[List[Branch]]:
    """Turn a predicate incident to ``x`` into alternative bound sets.

    Returns one branch per run index, or ``None`` to abstain.  An empty list
    means the predicate is unsatisfiable, so the whole term is zero -- the
    caller gets that for free, since a product over an empty branch list has no
    combinations.

    For a left-continuous relation an extra branch is emitted that owns the
    anchor ``xi_j = 0`` alone.  There the predicate's truth does not depend on
    ``x`` at all: ``0 <= f(x)`` holds for every ``x`` because every value in
    this class is non-negative, and ``0 > f(x)`` holds for none.  So the anchor
    branch either drops the constraint entirely or does not exist, and in both
    cases the run branches disown the anchor.  The two families therefore
    partition ``xi_j``, which is what keeps the sum exact rather than
    double-counting the atom.
    """
    if p.j == x:
        # x is the constrained side: the predicate already reads x <| f(xi_i),
        # so f is the bound itself and needs no inversion.  f may have turns
        # here; nothing in the cumulative step requires a bound to be monotone.
        return [([Candidate(p.i, p.f, p.upper, not p.strict, "direct")], [])]

    assert p.i == x
    f = p.f
    right_cont = _right_continuous(p.upper, p.strict)
    # The threshold on x is a step function of the value v = xi_j, and it can
    # only change where v crosses an attained value of f.  Breakpoints must be
    # strictly positive (the argument xi_j lives on [0, inf)); the anchor v = 0
    # is handled by its own branch below rather than by a breakpoint, which a
    # right-continuous step could not place.
    breaks = sorted({v for v in f.vs if v not in (INF, -INF) and v > 0.0})
    samples = _samples_for(breaks, interior=not right_cont)
    runs_per_sample = [_feasible_runs(f, v, p.upper, p.strict) for v in samples]
    rmax = max((len(r) for r in runs_per_sample), default=0)
    if counters:
        counters.note_runs(rmax)

    bx = tuple(breaks)
    branches: List[Branch] = []
    for k in range(rmax):
        lows: List[float] = []
        ups: List[float] = []
        feas: List[float] = []
        for runs in runs_per_sample:
            if k < len(runs):
                L, _lc, U, _uc = runs[k]
                lows.append(L)
                ups.append(U)
                feas.append(1.0)
            else:
                # This piece has fewer runs than the branch index.  [0, 0) --
                # lower closed, upper open at the same point -- carries mass
                # Q(0, open) - Q(0, open) = 0 under Lebesgue and magnitude
                # alike, so the anchor atom is neither kept nor double-counted.
                lows.append(0.0)
                ups.append(0.0)
                feas.append(0.0)
        cands = [Candidate(p.j, Step(bx, tuple(lows)), False, True, "inv-lo"),
                 Candidate(p.j, Step(bx, tuple(ups)), True, False, "inv-up")]
        cands = [c for c in cands if not c.is_vacuous()]
        step = None if all(v == 1.0 for v in feas) else Step(bx, tuple(feas))
        # atom multiplier: the anchor belongs to the anchor branch when the
        # relation is left-continuous, and otherwise follows the step.
        atom_mult = 0.0 if not right_cont else (1.0 if step is None
                                                else step(0.0))
        mods = [(p.j, atom_mult, step)]
        branches.append((cands, mods))

    if not right_cont and p.upper:
        # `0 <= f(x)` holds for every x, so at the anchor the predicate imposes
        # nothing: no candidates, and the continuous part is handed to the run
        # branches above.
        branches.append(([], [(p.j, 1.0, step_const(0.0))]))
    return branches


# --------------------------------------------------------------------------- #
# One elimination
# --------------------------------------------------------------------------- #


def eliminate_variable(term: PairTerm, x: int, mu: Measure1D,
                       domain: Tuple[float, float],
                       grid: Sequence[float],
                       counters: Optional[Counters] = None
                       ) -> Optional[List[PairTerm]]:
    """Integrate ``x`` out of one term.

    Returns a signed list of terms in the remaining variables, or ``None`` if
    some predicate falls outside the implemented grammar (the caller then knows
    the constructor abstained rather than guessed).
    """
    lo, hi = domain
    rest: List[PairPredicate] = []
    branch_sets: List[List[Branch]] = []
    for p in term.preds:
        if x not in (p.i, p.j):
            rest.append(p)
            continue
        got = _invert_to_branches(p, x, counters)
        if got is None:
            return None
        branch_sets.append(got)

    if counters:
        counters.eliminations += 1
        combos = 1
        for bs in branch_sets:
            combos *= max(len(bs), 1)
        counters.note_combos(combos)

    out: List[PairTerm] = []
    n_labels = 0
    # The feasible set of x is the intersection of the incident predicates'
    # feasible sets, and intersecting unions of runs gives the union over
    # choices of one run from each -- hence the product.  Each combination is a
    # single intersection of intervals, which the active-bound labelling below
    # handles exactly as before.
    for combo in itertools.product(*branch_sets):
        base = term.clone()
        base.preds = rest
        lowers: List[Candidate] = [
            Candidate(None, step_const(lo), False, True, "dom")]
        uppers: List[Candidate] = [
            Candidate(None, step_const(hi), True, True, "dom")]
        for cands, mods in combo:
            for axis, atom_mult, step in mods:
                base.mul_split(axis, atom_mult, step)
            for c in cands:
                (uppers if c.upper else lowers).append(c)
        if base.is_dead():
            continue

        qx = base.weights.get(x, step_const(1.0))
        atom_x = base.atom_weight(x)
        events = sorted(set(grid) | set(qx.xs))

        def Q(t: float, closed: bool, _qx=qx, _a=atom_x, _ev=events) -> float:
            return _cum_atom(_qx.__call__, _a, _ev, mu, t, closed=closed)

        for lo_c, up_c in itertools.product(lowers, uppers):
            n_labels += 1
            # Q(U) - Q(L^-) on the region where these two are the active bounds
            for cand, sign in ((up_c, 1.0), (lo_c, -1.0)):
                seed = base.clone()
                seed.coeff *= sign
                seed.weights.pop(x, None)
                seed.atoms.pop(x, None)        # x is consumed by Q, not stored
                seed.preds = list(rest)
                pieces = _apply_label(seed, lo_c, up_c, lowers, uppers)
                if pieces is None:
                    return None                # abstention propagates outward
                for piece in pieces:
                    if not _absorb_cumulative(piece, cand, Q, hi):
                        continue
                    piece.preds = merge_same_orientation(piece.preds)
                    if not piece.is_dead():
                        out.append(piece)
    if counters:
        counters.record(len(out), n_labels)
    return out


def _apply_label(piece: PairTerm, lo_c: Candidate, up_c: Candidate,
                 lowers: Sequence[Candidate],
                 uppers: Sequence[Candidate]) -> Optional[List[PairTerm]]:
    """Impose "these two candidates are the active bounds", plus feasibility.

    The active lower bound is the pointwise maximum of the lower candidates and
    the active upper the pointwise minimum, so ties must be broken or the same
    point is labelled twice.  We use the *first winner*: among candidates
    attaining the extremum, the one with the smallest index owns the region.
    Hence ``lo_c`` at index ``k`` requires a strict inequality against every
    earlier candidate and a non-strict one against every later candidate.
    Without this, two numerically equal bounds -- which inversion produces
    routinely, since distinct predicates often invert to the same domain limit
    -- each claim the region and the integral comes out 2x or 4x too large.

    A single condition can split the region it constrains (see
    ``_comparison_predicate``), so this threads a *list* of pieces through the
    conditions, fanning out as needed.  An empty list means the label is
    identically false.
    """
    conditions: List[Tuple[Candidate, Candidate, bool]] = []
    k_lo = next(i for i, c in enumerate(lowers) if c is lo_c)
    for i, other in enumerate(lowers):
        if other is not lo_c:
            conditions.append((other, lo_c, i < k_lo))
    k_up = next(i for i, c in enumerate(uppers) if c is up_c)
    for i, other in enumerate(uppers):
        if other is not up_c:
            conditions.append((up_c, other, i < k_up))
    conditions.append((lo_c, up_c, False))         # feasibility L <= U

    pieces = [piece]
    for a, b, strict in conditions:
        nxt: List[PairTerm] = []
        for pc in pieces:
            got = _impose(pc, a, b, strict)
            if got is None:
                return None
            nxt.extend(got)
        pieces = nxt
        if not pieces:
            return []
    return pieces


def _impose(piece: PairTerm, a: Candidate, b: Candidate,
            strict: bool) -> Optional[List[PairTerm]]:
    """Impose ``a(.) < b(.)`` (strict) or ``a(.) <= b(.)``.

    Returns the resulting pieces (one, or several when the condition splits),
    the empty list when it is identically false, and ``None`` when it cannot be
    expressed in the term class -- which the caller must propagate as an
    abstention, never swallow.
    """
    if a.var is None and b.var is None:
        held = (a.f(0.0) < b.f(0.0)) if strict else (a.f(0.0) <= b.f(0.0))
        return [piece] if held else []
    if a.var is None:                                  # const <(=) b(xi_q)
        piece.mul_weight(b.var, _indicator_const(b.f, a.f(0.0), True, strict))
        return [] if piece.is_dead() else [piece]
    if b.var is None:                                  # a(xi_p) <(=) const
        piece.mul_weight(a.var, _indicator_const(a.f, b.f(0.0), False, strict))
        return [] if piece.is_dead() else [piece]
    if a.var == b.var:                                 # same variable: filter
        piece.mul_weight(a.var, _indicator_pointwise(a.f, b.f, strict))
        return [] if piece.is_dead() else [piece]
    return _comparison_predicate(piece, a, b, strict)


def _indicator_const(f: Step, c: float, const_is_lower: bool,
                     strict: bool) -> Step:
    """Indicator of ``c <(=) f(.)`` or of ``f(.) <(=) c``."""
    if const_is_lower:
        ok = (lambda v: c < v) if strict else (lambda v: c <= v)
    else:
        ok = (lambda v: v < c) if strict else (lambda v: v <= c)
    return Step(f.xs, tuple(1.0 if ok(v) else 0.0 for v in f.vs))


def _indicator_pointwise(a: Step, b: Step, strict: bool) -> Step:
    ok = (lambda u, v: u < v) if strict else (lambda u, v: u <= v)
    xs = tuple(sorted(set(a.xs) | set(b.xs)))
    vs = [1.0 if ok(a(0.0), b(0.0)) else 0.0]
    for t in xs:
        vs.append(1.0 if ok(a(t), b(t)) else 0.0)
    return Step(xs, tuple(vs))


def _comparison_predicate(piece: PairTerm, a: Candidate, b: Candidate,
                          strict: bool) -> Optional[List[PairTerm]]:
    """Impose ``a(xi_p) <(=) b(xi_q)`` as pair predicates on ``(p, q)``.

    For each piece of ``a`` the admissible ``xi_q`` set is the superlevel set of
    ``b`` at the level ``a`` takes there.  That set is a union of runs, so the
    condition becomes a signed sum over run indices, each expressed as up to two
    pair predicates with a ``0/1`` filter on ``xi_p`` for the pieces where that
    run is absent.  The runs are disjoint, so the pieces returned here partition
    the incoming one and no mass is counted twice.
    """
    cuts = [c for c in sorted(set(a.f.xs)) if c > 0.0]
    reps = [0.0] + cuts                      # right-continuity: xs[k-1] is a rep
    runs_per = [_feasible_runs(b.f, a.f(t), True, strict) for t in reps]
    rmax = max((len(r) for r in runs_per), default=0)
    if rmax == 0:
        return []                            # identically false

    cx = tuple(cuts)
    out: List[PairTerm] = []
    for k in range(rmax):
        lows: List[float] = []
        ups: List[float] = []
        feas: List[float] = []
        for runs in runs_per:
            if k < len(runs):
                L, _lc, U, _uc = runs[k]
                lows.append(L)
                ups.append(U)
                feas.append(1.0)
            else:
                # Absent run: make it unsatisfiable both ways.  The filter
                # already zeroes the piece; the infinite lower bound is belt and
                # braces, and never reaches a weight, only an indicator.
                lows.append(INF)
                ups.append(0.0)
                feas.append(0.0)
        pc = piece.clone() if rmax > 1 else piece
        if not all(v == 1.0 for v in feas):
            pc.mul_weight(a.var, Step(cx, tuple(feas)))
        if not all(v == 0.0 for v in lows):
            pc.preds.append(PairPredicate(a.var, b.var, Step(cx, tuple(lows)),
                                          upper=False, strict=False))
        if not all(v == INF for v in ups):
            pc.preds.append(PairPredicate(a.var, b.var, Step(cx, tuple(ups)),
                                          upper=True, strict=True))
        if not pc.is_dead():
            out.append(pc)
    return out


def _absorb_cumulative(piece: PairTerm, cand: Candidate, Q, hi: float) -> bool:
    """Multiply in the cumulative evaluated at this candidate.

    The mass of the feasible interval is ``Q(U) - Q(L^-)``.  So the endpoint
    flag has to be *inverted* on the lower side: if ``L`` is attained, what must
    be subtracted is the mass strictly below it, and if ``L`` is excluded, the
    mass up to and including it.  Under Lebesgue this makes no difference; under
    the magnitude measure at ``L = 0`` it decides whether the anchor atom is
    kept or deleted, and getting it backwards silently loses mass 1.

    Endpoints are clamped to the ground set's ceiling ``hi`` before the measure
    is taken.  An inverted bound carries ``+inf`` on the pieces where its
    predicate constrains nothing, and the cumulative of an unbounded endpoint
    diverges: ``_cum`` integrates ``density * (inf - g) * weight(inf)``, which
    is ``inf`` when that last value is positive and ``0 * inf = NaN`` when it is
    zero.  Those pieces are already excluded by the label indicators, but a NaN
    is not annihilated by multiplying with zero, so it has to be prevented
    rather than cancelled.  Clamping is exact here because ``x`` is confined to
    ``[lo, hi]``: "unconstrained above" *is* ``hi``.

    This clamp is at the measure step only.  Clamping the candidate itself,
    before the label comparisons, is wrong -- there ``+inf`` correctly means
    "weaker than every other bound", and collapsing it onto ``hi`` makes a
    vacuous bound tie with the ground set and claim regions it does not own.
    An earlier attempt clamped both sides everywhere and broke p=4 and p=5.
    """
    closed = cand.closed if cand.upper else (not cand.closed)
    if cand.var is None:
        piece.coeff *= Q(min(cand.f(0.0), hi), closed)
        return piece.coeff != 0.0
    xs = cand.f.xs
    vs = tuple(Q(min(v, hi), closed) for v in cand.f.vs)
    piece.mul_weight(cand.var, Step(xs, vs))
    return True


# --------------------------------------------------------------------------- #
# Full elimination
# --------------------------------------------------------------------------- #


def eliminate_all(terms: Sequence[PairTerm], order: Sequence[int],
                  measures: Dict[int, Measure1D],
                  domains: Dict[int, Tuple[float, float]],
                  grids: Dict[int, Sequence[float]],
                  counters: Optional[Counters] = None,
                  term_cap: Optional[int] = None
                  ) -> Optional[float]:
    """Integrate every variable in ``order`` out; returns the scalar value.

    ``term_cap`` aborts (as an abstention) when the state exceeds that many
    terms, so a blow-up is reported as "not measured" rather than running the
    machine out of memory.
    """
    state = [t.clone() for t in terms]
    for x in order:
        nxt: List[PairTerm] = []
        for t in state:
            got = eliminate_variable(t, x, measures[x], domains[x], grids[x],
                                     counters)
            if got is None:
                return None
            nxt.extend(got)
        state = nxt
        if counters:
            counters.peak_state = max(counters.peak_state, len(state))
        if term_cap is not None and len(state) > term_cap:
            return None
        if not state:
            return 0.0
    return sum(t.coeff for t in state)

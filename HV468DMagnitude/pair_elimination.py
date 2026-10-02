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
decorated inversion, a *decorated interval* of feasible ``x``, whose two
endpoints are unary functions of exactly one surviving variable.  With ``O(p)``
candidates per side, the surviving space is partitioned by which candidate is
active; on each label the integral is ``Q(U) - Q(L^-)`` for the unary cumulative
``Q`` of ``x``'s own density, and every label condition compares two unary maps,
which is again a pair predicate (or a unary filter when both sides involve the
same variable).

Interval-valued inversion, and why one side is not enough
---------------------------------------------------------
Inverting ``xi_j <| f(x)`` into a constraint on ``x`` cannot in general produce
a single-sided bound:

* when ``f`` is monotone the feasible ``x``-set is a prefix or a suffix, but the
  *degenerate* pieces (all ``x`` feasible, no ``x`` feasible) belong to the
  opposite family, so forcing one side misrepresents them;
* when ``f`` has one turn -- which this class provably contains, see
  ``pair_state`` -- the superlevel set of a hill is a genuine middle interval
  with two nontrivial endpoints.

So inversion returns a lower *and* an upper candidate, plus a 0/1 feasibility
weight marking the pieces where the interval is empty.  A set of ``x`` that is
not an interval at all (the superlevel set of a valley) is outside the grammar
and the constructor abstains.

Abstention must be total
------------------------
Every "cannot express this" path returns ``None`` all the way out.  An earlier
version let one of them -- a comparison against a non-monotone bound -- fall
through as *unconstrained*, which silently deletes a constraint and inflates the
integral.  Dropping a constraint is not a conservative approximation; it is a
wrong answer that looks like a passing test until the term class is wide enough
to reach it.  ``_impose`` therefore has three outcomes: imposed, identically
false, and abstain.
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

#: ``_feasible_interval`` returns this when the feasible set is empty.  It is a
#: distinct object rather than ``None`` so that "empty" can never be confused
#: with "outside the grammar", and never with the numeric value 0.
EMPTY_IV = "empty"


class Counters:
    """Instrumentation required by the prompt.

    ``max_branching`` is the quantity the lemma actually bounds: the largest
    number of terms produced from *one* incoming term by *one* elimination.
    A running total would conflate that with the depth of the chain.
    """

    def __init__(self):
        self.labels = 0
        self.terms_in = 0
        self.terms_out = 0
        self.eliminations = 0
        self.max_branching = 0
        self.max_labels_per_elim = 0

    def record(self, produced: int, labels: int) -> None:
        self.terms_in += 1
        self.terms_out += produced
        self.max_branching = max(self.max_branching, produced)
        self.max_labels_per_elim = max(self.max_labels_per_elim, labels)

    def __repr__(self):
        return (f"Counters(elims={self.eliminations}, "
                f"max_branching={self.max_branching}, "
                f"max_labels={self.max_labels_per_elim})")


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


def _feasible_interval(f: Step, v: float, upper: bool, strict: bool):
    """Feasible ``x``-set of ``v <| f(x)`` as a decorated interval.

    Returns ``(L, L_closed, U, U_closed)``, or ``EMPTY_IV``, or ``None`` when
    the set is not an interval (so the caller must abstain).

    ``Step`` is right-continuous: piece ``k >= 1`` of ``f`` spans
    ``[xs[k-1], xs[k])`` and piece 0 spans ``[0, xs[0])``.  Hence a feasible run
    starting at piece ``k`` is closed at ``xs[k-1]``, and one ending at piece
    ``k`` is open at ``xs[k]``.  The decorations are therefore uniform by
    construction -- lower closed, upper open -- which is what lets the empty
    pieces be encoded as ``[0, 0)`` with exactly zero mass under either measure.
    """
    if upper:
        ok = (lambda t: t > v) if strict else (lambda t: t >= v)
    else:
        ok = (lambda t: t < v) if strict else (lambda t: t <= v)
    mask = [ok(val) for val in f.vs]
    if not any(mask):
        return EMPTY_IV
    first = mask.index(True)
    last = len(mask) - 1 - mask[::-1].index(True)
    if not all(mask[first:last + 1]):
        return None                       # two or more runs: not an interval
    L = 0.0 if first == 0 else f.xs[first - 1]
    U = INF if last == len(mask) - 1 else f.xs[last]
    return (L, True, U, False)


def _samples_for(breaks: Sequence[float]) -> List[float]:
    """One representative per piece of a step function with these breakpoints."""
    if not breaks:
        return [0.0]
    out = [0.0]                                       # the piece [0, breaks[0])
    for a, b in zip(breaks, breaks[1:]):
        out.append(0.5 * (a + b))
    out.append(breaks[-1] + 1.0)
    return out


def _invert_to_bounds(p: PairPredicate, x: int):
    """Turn a predicate incident to ``x`` into bounds on ``x``.

    Returns ``(candidates, feasibility_weight, feasibility_axis)`` or ``None``
    to abstain.  ``candidates`` holds at most one lower and one upper bound;
    the feasibility weight, when present, is a 0/1 step on the surviving
    variable that is 0 exactly where the feasible interval is empty.
    """
    if p.j == x:
        # x is the constrained side: the predicate already reads x <| f(xi_i),
        # so f is the bound itself and needs no inversion.  f may have a turn
        # here; nothing in the cumulative step requires a bound to be monotone.
        return ([Candidate(p.i, p.f, p.upper, not p.strict, "direct")],
                None, None)

    assert p.i == x
    f = p.f
    # The threshold on x is a step function of the value v = xi_j, and it can
    # only change where v crosses an attained value of f.  Breakpoints must be
    # strictly positive (the argument xi_j lives on [0, inf)), so a value of 0
    # contributes no breakpoint -- the piece below it is unreachable.
    breaks = sorted({v for v in f.vs if v not in (INF, -INF) and v > 0.0})
    samples = _samples_for(breaks)

    lows: List[float] = []
    ups: List[float] = []
    feas: List[float] = []
    for v in samples:
        iv = _feasible_interval(f, v, p.upper, p.strict)
        if iv is None:
            return None                   # not an interval: outside the grammar
        if iv == EMPTY_IV:
            # [0, 0) -- lower closed at 0, upper open at 0 -- carries mass
            # Q(0, open) - Q(0, open) = 0 under Lebesgue and under magnitude
            # alike, so the anchor atom is neither kept nor double-counted.
            lows.append(0.0)
            ups.append(0.0)
            feas.append(0.0)
            continue
        L, _lc, U, _uc = iv
        lows.append(L)
        ups.append(U)
        feas.append(1.0)

    bx = tuple(breaks)
    cands = [Candidate(p.j, Step(bx, tuple(lows)), False, True, "inv-lo"),
             Candidate(p.j, Step(bx, tuple(ups)), True, False, "inv-up")]
    cands = [c for c in cands if not c.is_vacuous()]
    weight = None if all(v == 1.0 for v in feas) else Step(bx, tuple(feas))
    return (cands, weight, p.j)


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
    lowers: List[Candidate] = [Candidate(None, step_const(lo), False, True, "dom")]
    uppers: List[Candidate] = [Candidate(None, step_const(hi), True, True, "dom")]

    base = term.clone()
    rest: List[PairPredicate] = []
    for p in term.preds:
        if x not in (p.i, p.j):
            rest.append(p)
            continue
        got = _invert_to_bounds(p, x)
        if got is None:
            return None
        cands, weight, axis = got
        if weight is not None:
            base.mul_weight(axis, weight)
        for c in cands:
            (uppers if c.upper else lowers).append(c)
    base.preds = rest

    if counters:
        counters.eliminations += 1
    if base.is_dead():
        if counters:
            counters.record(0, 0)
        return []

    qx = base.weights.get(x, step_const(1.0))
    events = sorted(set(grid) | set(qx.xs))

    def Q(t: float, closed: bool) -> float:
        return _cum(qx.__call__, events, mu, t, closed=closed)

    out: List[PairTerm] = []
    n_labels = 0
    for lo_c, up_c in itertools.product(lowers, uppers):
        n_labels += 1
        # Q(U) - Q(L^-) on the region where these two are the active bounds
        for cand, sign in ((up_c, 1.0), (lo_c, -1.0)):
            piece = base.clone()
            piece.coeff *= sign
            piece.weights.pop(x, None)
            piece.preds = list(rest)
            ok = _apply_label(piece, lo_c, up_c, lowers, uppers)
            if ok is None:
                return None                    # abstention propagates outward
            if not ok:
                continue
            if not _absorb_cumulative(piece, cand, Q):
                continue
            piece.preds = merge_same_orientation(piece.preds)
            if not piece.is_dead():
                out.append(piece)
    if counters:
        counters.record(len(out), n_labels)
    return out


def _apply_label(piece: PairTerm, lo_c: Candidate, up_c: Candidate,
                 lowers: Sequence[Candidate],
                 uppers: Sequence[Candidate]) -> Optional[bool]:
    """Impose "these two candidates are the active bounds", plus feasibility.

    The active lower bound is the pointwise maximum of the lower candidates and
    the active upper the pointwise minimum, so ties must be broken or the same
    point is labelled twice.  We use the *first winner*: among candidates
    attaining the extremum, the one with the smallest index owns the region.
    Hence ``lo_c`` at index ``k`` requires a strict inequality against every
    earlier candidate and a non-strict one against every later candidate.
    Without this, two numerically equal bounds -- which interval inversion
    produces routinely, since distinct predicates often invert to the same
    domain limit -- each claim the region and the integral comes out 2x or 4x
    too large.
    """
    k_lo = next(i for i, c in enumerate(lowers) if c is lo_c)
    for i, other in enumerate(lowers):
        if other is lo_c:
            continue
        got = _impose(piece, other, lo_c, strict=(i < k_lo))
        if got is None or not got:
            return got
    k_up = next(i for i, c in enumerate(uppers) if c is up_c)
    for i, other in enumerate(uppers):
        if other is up_c:
            continue
        got = _impose(piece, up_c, other, strict=(i < k_up))
        if got is None or not got:
            return got
    return _impose(piece, lo_c, up_c, strict=False)    # feasibility L <= U


def _impose(piece: PairTerm, a: Candidate, b: Candidate,
            strict: bool) -> Optional[bool]:
    """Impose ``a(.) < b(.)`` (strict) or ``a(.) <= b(.)``.

    Returns ``True`` when imposed, ``False`` when identically false, and
    ``None`` when the condition cannot be expressed in the term class -- which
    the caller must propagate as an abstention, never swallow.
    """
    if a.var is None and b.var is None:
        return (a.f(0.0) < b.f(0.0)) if strict else (a.f(0.0) <= b.f(0.0))
    if a.var is None:                                  # const <(=) b(xi_q)
        piece.mul_weight(b.var, _indicator_const(b.f, a.f(0.0), True, strict))
        return True
    if b.var is None:                                  # a(xi_p) <(=) const
        piece.mul_weight(a.var, _indicator_const(a.f, b.f(0.0), False, strict))
        return True
    if a.var == b.var:                                 # same variable: filter
        piece.mul_weight(a.var, _indicator_pointwise(a.f, b.f, strict))
        return True
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
                          strict: bool) -> Optional[bool]:
    """Impose ``a(xi_p) <(=) b(xi_q)`` as pair predicates on ``(p, q)``.

    For each piece of ``a`` the admissible ``xi_q`` set is the superlevel set of
    ``b`` at the level ``a`` takes there.  That set is an interval whenever
    ``b`` is monotone or has a single hill, so it is expressed as up to two pair
    predicates; where it is empty, a 0/1 filter on ``xi_p`` records it.  If it
    is not an interval the condition is inexpressible and we abstain -- we do
    *not* treat it as unconstrained.
    """
    cuts = [c for c in sorted(set(a.f.xs)) if c > 0.0]
    reps = [0.0] + cuts                      # right-continuity: xs[k-1] is a rep
    lows: List[float] = []
    ups: List[float] = []
    feas: List[float] = []
    for t in reps:
        iv = _feasible_interval(b.f, a.f(t), True, strict)
        if iv is None:
            return None
        if iv == EMPTY_IV:
            lows.append(INF)
            ups.append(0.0)
            feas.append(0.0)
            continue
        L, _lc, U, _uc = iv
        lows.append(L)
        ups.append(U)
        feas.append(1.0)

    cx = tuple(cuts)
    if not all(v == 1.0 for v in feas):
        piece.mul_weight(a.var, Step(cx, tuple(feas)))
    if not any(v == 1.0 for v in feas):
        return True                          # the filter already zeroes it all
    if not all(v == 0.0 for v in lows):
        piece.preds.append(PairPredicate(a.var, b.var, Step(cx, tuple(lows)),
                                         upper=False, strict=False))
    if not all(v == INF for v in ups):
        piece.preds.append(PairPredicate(a.var, b.var, Step(cx, tuple(ups)),
                                         upper=True, strict=True))
    return True


def _absorb_cumulative(piece: PairTerm, cand: Candidate, Q) -> bool:
    """Multiply in the cumulative evaluated at this candidate.

    The mass of the feasible interval is ``Q(U) - Q(L^-)``.  So the endpoint
    flag has to be *inverted* on the lower side: if ``L`` is attained, what must
    be subtracted is the mass strictly below it, and if ``L`` is excluded, the
    mass up to and including it.  Under Lebesgue this makes no difference; under
    the magnitude measure at ``L = 0`` it decides whether the anchor atom is
    kept or deleted, and getting it backwards silently loses mass 1.
    """
    closed = cand.closed if cand.upper else (not cand.closed)
    if cand.var is None:
        piece.coeff *= Q(cand.f(0.0), closed)
        return piece.coeff != 0.0
    xs = cand.f.xs
    vs = tuple(Q(v, closed) for v in cand.f.vs)
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
        if term_cap is not None and len(state) > term_cap:
            return None
        if not state:
            return 0.0
    return sum(t.coeff for t in state)

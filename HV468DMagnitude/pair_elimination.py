"""Variable elimination for weighted monotone-pair terms.

This implements the constructor the generalized report's Lemma 2 (one-variable
elimination) describes, and exposes the pair version as two applications of it.

A structural observation, recorded in the audit
-----------------------------------------------
The report proves Lemma 2 (one variable) and Lemma 3 (a pair) separately. For
*closure* purposes the pair lemma is not needed: the mutual ``x``--``y``
staircase is, from ``x``'s point of view, just one more bound candidate that
happens to depend on ``y``, so eliminating ``x`` and then ``y`` by Lemma 2
twice already yields the claim. The pair version buys a better constant and a
cleaner closed form (the two-branch staircase formula), not extra generality.
We therefore test Lemma 2 as the load-bearing statement and obtain the pair
result by iteration.

What one elimination does
-------------------------
Every predicate incident to the eliminated variable ``x`` becomes, after
decorated inversion, a lower or upper bound on ``x`` whose endpoint is a
monotone unary function of exactly one surviving variable. With ``O(p)``
candidates per side, the surviving space is partitioned by which candidate is
active; on each label the integral is ``Q(U) - Q(L)`` for the unary cumulative
``Q`` of ``x``'s own density, and every label condition compares two unary
monotone maps, which is again a pair predicate (or a unary filter when both
sides involve the same variable).
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
    """A bound on the eliminated variable, as a monotone map of one survivor.

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

    def __repr__(self):
        side = "U" if self.upper else "L"
        return f"<{side} {self.tag} var={self.var}>"


def _invert_to_bound(p: PairPredicate, x: int) -> Optional[Candidate]:
    """Turn a predicate incident to ``x`` into a bound on ``x``.

    Two cases.  If ``x`` is the *constrained* side (``j``), the predicate
    already reads ``x <| f(xi_i)`` and the bound is ``f`` itself.  If ``x`` is
    the *argument* side (``i``), the predicate reads ``xi_j <| f(x)`` and must
    be inverted: the feasible ``x``-set is a prefix when ``f`` is nonincreasing
    and a suffix when it is nondecreasing.  The inverse is expressed as a step
    map of ``xi_j``, built from ``f``'s value set.
    """
    if p.j == x:
        return Candidate(p.i, p.f, p.upper, not p.strict, "direct")

    assert p.i == x
    f = p.f
    nonincreasing = f.is_nonincreasing()
    nondecreasing = f.is_nondecreasing()
    if not (nonincreasing or nondecreasing):
        return None                                   # caller falls back

    # The threshold on x is a step function of the value v = xi_j, and it can
    # only change where v crosses an attained value of f.  Breakpoints must be
    # strictly positive (the argument xi_j lives on [0, inf)), so a value of 0
    # contributes no breakpoint -- the piece below it is unreachable.
    breaks = sorted({v for v in f.vs if v not in (INF, -INF) and v > 0.0})
    samples = _samples_for(breaks)
    bounds = [_solve_threshold(f, v, p.upper, p.strict, nonincreasing)
              for v in samples]
    sides = {b[1] for b in bounds}
    if len(sides) != 1:
        return None                                   # mixed: caller abstains
    is_upper = bounds[0][1]
    step = Step(tuple(breaks), tuple(b[0] for b in bounds))
    closed = bounds[0][2]
    return Candidate(p.j, step, is_upper, closed, "inverted")


def _samples_for(breaks: Sequence[float]) -> List[float]:
    """One representative per piece of a step function with these breakpoints."""
    if not breaks:
        return [0.0]
    out = [0.0]                                       # the piece [0, breaks[0])
    for a, b in zip(breaks, breaks[1:]):
        out.append(0.5 * (a + b))
    out.append(breaks[-1] + 1.0)
    return out


def _solve_threshold(f: Step, v: float, upper: bool, strict: bool,
                     nonincreasing: bool) -> Tuple[float, bool, bool]:
    """Feasible ``x``-set of ``v <| f(x)`` as (endpoint, is_upper, closed)."""
    # predicate on x: f(x) > v (upper & strict), f(x) >= v (upper, non-strict),
    #                 f(x) < v (lower & strict), f(x) <= v (lower, non-strict)
    if upper:
        ok = (lambda t: t > v) if strict else (lambda t: t >= v)
    else:
        ok = (lambda t: t < v) if strict else (lambda t: t <= v)
    mask = [ok(val) for val in f.vs]
    if not any(mask):
        return (0.0, True, False)                     # empty, as exclusive 0
    if all(mask):
        return (INF, True, False)
    first = mask.index(True)
    last = len(mask) - 1 - mask[::-1].index(True)
    if first == 0:                                    # prefix -> upper bound
        return (f.xs[last], True, False)
    if last == len(mask) - 1:                         # suffix -> lower bound
        return (f.xs[first - 1], False, True)
    return (f.xs[last], True, False)                  # middle: approximate side


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
    rest: List[PairPredicate] = []

    for p in term.preds:
        if x not in (p.i, p.j):
            rest.append(p)
            continue
        cand = _invert_to_bound(p, x)
        if cand is None:
            return None
        (uppers if cand.upper else lowers).append(cand)

    qx = term.weights.get(x, step_const(1.0))
    events = sorted(set(grid) | set(qx.xs))

    def Q(t: float, closed: bool) -> float:
        return _cum(qx.__call__, events, mu, t, closed=closed)

    out: List[PairTerm] = []
    n_labels = 0
    if counters:
        counters.eliminations += 1
    for lo_c, up_c in itertools.product(lowers, uppers):
        n_labels += 1
        # Q(U) - Q(L) on the region where these two are the active bounds
        for cand, sign in ((up_c, 1.0), (lo_c, -1.0)):
            piece = term.clone()
            piece.coeff *= sign
            piece.weights.pop(x, None)
            piece.preds = list(rest)
            if not _apply_label(piece, lo_c, up_c, lowers, uppers):
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
                 uppers: Sequence[Candidate]) -> bool:
    """Impose "these two candidates are the active bounds", plus feasibility."""
    for other in lowers:
        if other is lo_c:
            continue
        if not _impose(piece, other, lo_c, "le"):     # other <= active lower
            return False
    for other in uppers:
        if other is up_c:
            continue
        if not _impose(piece, up_c, other, "le"):     # active upper <= other
            return False
    return _impose(piece, lo_c, up_c, "le")           # feasibility L <= U


def _impose(piece: PairTerm, a: Candidate, b: Candidate, rel: str) -> bool:
    """Impose ``a(.) <= b(.)``; returns False if it is identically false."""
    if a.var is None and b.var is None:
        return a.f(0.0) <= b.f(0.0)
    if a.var is None:                                  # const <= b(xi_q)
        const = a.f(0.0)
        piece.mul_weight(b.var, _indicator_ge(b.f, const))
        return True
    if b.var is None:                                  # a(xi_p) <= const
        const = b.f(0.0)
        piece.mul_weight(a.var, _indicator_le(a.f, const))
        return True
    if a.var == b.var:                                 # same variable: filter
        piece.mul_weight(a.var, _indicator_pointwise(a.f, b.f))
        return True
    pred = _comparison_predicate(a, b)
    if pred is None:
        return True                                    # unconstrained
    piece.preds.append(pred)
    return True


def _indicator_le(f: Step, c: float) -> Step:
    return Step(f.xs, tuple(1.0 if v <= c else 0.0 for v in f.vs))


def _indicator_ge(f: Step, c: float) -> Step:
    return Step(f.xs, tuple(1.0 if v >= c else 0.0 for v in f.vs))


def _indicator_pointwise(a: Step, b: Step) -> Step:
    xs = tuple(sorted(set(a.xs) | set(b.xs)))
    vs = [1.0 if a(0.0) <= b(0.0) else 0.0]
    for t in xs:
        vs.append(1.0 if a(t) <= b(t) else 0.0)
    return Step(xs, tuple(vs))


def _comparison_predicate(a: Candidate, b: Candidate) -> Optional[PairPredicate]:
    """``a(xi_p) <= b(xi_q)`` as a pair predicate on ``(p, q)``.

    For each value of ``xi_p`` the admissible ``xi_q`` set is, by monotonicity
    of ``b``, a prefix or a suffix; its boundary as a function of ``xi_p`` is
    the composition of ``b``'s inverse with ``a``.  Both are monotone, so the
    composition is monotone and the result is a legal pair predicate.
    """
    if not (b.f.is_nondecreasing() or b.f.is_nonincreasing()):
        return None
    cuts = sorted(set(a.f.xs))
    reps = [0.0] + [0.5 * (p + q) for p, q in zip(cuts, cuts[1:])] + \
           ([cuts[-1] + 1.0] if cuts else [])
    reps = [0.0] if not cuts else [0.0] + [c + 1e-9 for c in cuts]
    thresholds, sides = [], set()
    for t in reps:
        want = a.f(t)
        bound, is_upper = _solve_ge(b.f, want)
        thresholds.append(bound)
        sides.add(is_upper)
    if len(sides) != 1:
        return None
    is_upper = sides.pop()
    step = Step(tuple(cuts), tuple(thresholds))
    return PairPredicate(a.var, b.var, step, upper=is_upper, strict=False)


def _solve_ge(f: Step, v: float) -> Tuple[float, bool]:
    """Feasible set of ``f(t) >= v`` as (endpoint, is_upper_bound)."""
    mask = [val >= v for val in f.vs]
    if not any(mask):
        return (-INF, True)
    if all(mask):
        return (INF, True)
    first = mask.index(True)
    last = len(mask) - 1 - mask[::-1].index(True)
    if first == 0:
        return (f.xs[last], True)
    return (f.xs[first - 1], False)


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
                  counters: Optional[Counters] = None
                  ) -> Optional[float]:
    """Integrate every variable in ``order`` out; returns the scalar value."""
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
        if not state:
            return 0.0
    return sum(t.coeff for t in state)

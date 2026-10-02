"""Weighted monotone-pair terms: the state class of the generalized report.

A term in ``p`` variables is

    Phi(xi) = coeff * prod_i q_i(xi_i) * prod_e [ xi_{j_e} <| f_e(xi_{i_e}) ]

with each ``q_i`` a unary step density, each ``f_e`` a monotone step map, and
``<|`` a decorated strict/non-strict comparison.  A *state* is a signed sum of
such terms.

Deviations from Definition 1 of the report, recorded deliberately
------------------------------------------------------------------
The report's definition keeps "at most one lower and one upper monotone
predicate per unordered pair, since same-orientation predicates are merged by
min/max".  That merge is only available when the two maps have the *same*
orientation.  Elimination can produce, for one pair, an upper predicate whose
map is nondecreasing *and* one whose map is nonincreasing; their conjunction is

    xi_j <| min( h_up(xi_i), h_down(xi_i) )

which has one turning point and is **not monotone**.  The four-dimensional
treatment allows exactly this (its "latent one-turn normal form"); the general
definition as written does not.

We therefore let a pair carry a *list* of predicates rather than one per side.
The list length stays bounded by a function of ``p`` -- each elimination adds
only O_p(1) of them -- so nothing in the asymptotics changes.  This is a repair
of the class definition, not a counterexample; see ``notes/audit.md``.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import INF, Measure1D, Step          # the audited 4-D primitives

__all__ = ["DecoratedBound", "PairPredicate", "PairTerm", "EMPTY", "FULL",
           "step_const", "merge_same_orientation"]


# --------------------------------------------------------------------------- #
# Decorated bounds
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DecoratedBound:
    """An endpoint carrying its own semantics.

    ``kind`` is ``"empty"``, ``"finite"`` or ``"full"``.  The numeric value 0 is
    *never* used to mean "empty": under the magnitude measure the singleton
    ``{0}`` has mass 1 while the empty set has mass 0, and conflating them was a
    real source of error in the 4-D audit.
    """

    kind: str
    value: float = 0.0
    closed: bool = False

    @staticmethod
    def empty() -> "DecoratedBound":
        return DecoratedBound("empty")

    @staticmethod
    def full() -> "DecoratedBound":
        return DecoratedBound("full")

    @staticmethod
    def finite(value: float, closed: bool) -> "DecoratedBound":
        return DecoratedBound("finite", value, closed)

    def is_empty(self) -> bool:
        return self.kind == "empty"

    def is_full(self) -> bool:
        return self.kind == "full"

    def as_upper(self) -> Tuple[float, bool]:
        """Numeric upper limit and whether it is attained."""
        if self.kind == "empty":
            return (-INF, False)
        if self.kind == "full":
            return (INF, False)
        return (self.value, self.closed)

    def as_lower(self) -> Tuple[float, bool]:
        if self.kind == "empty":
            return (INF, False)
        if self.kind == "full":
            return (0.0, True)              # the ground set starts at the anchor
        return (self.value, self.closed)


EMPTY = DecoratedBound.empty()
FULL = DecoratedBound.full()


def step_const(value: float) -> Step:
    return Step((), (value,))


# --------------------------------------------------------------------------- #
# Predicates and terms
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PairPredicate:
    """``xi_j  <|  f(xi_i)`` (upper) or ``xi_j  |>  f(xi_i)`` (lower).

    ``strict`` records whether the comparison is strict.  Under Lebesgue the
    flag is irrelevant; under the magnitude measure it decides whether the
    anchor atom is inside the region, so it is carried explicitly.
    """

    i: int
    j: int
    f: Step
    upper: bool = True
    strict: bool = True

    def holds(self, xi_i: float, xi_j: float) -> bool:
        bound = self.f(xi_i)
        if self.upper:
            return xi_j < bound if self.strict else xi_j <= bound
        return xi_j > bound if self.strict else xi_j >= bound

    def key(self):
        return (self.i, self.j, self.f.xs, self.f.vs, self.upper, self.strict)


@dataclass
class PairTerm:
    """One weighted monotone-pair term.

    The weight on an axis has **two** components, because the magnitude measure
    does: ``mu_M = delta_0 + (1/2) lambda``.  ``weights[axis]`` is the density
    against the continuous part, and ``atoms[axis]`` the weight at the anchor
    point ``0``, stored absolutely rather than as a multiplier.  When no atom
    weight is recorded the default is ``1.0``, so a term that never touches the
    anchor behaves exactly as before.

    Why they have to be independent
    -------------------------------
    A right-continuous ``Step`` cannot distinguish the single point ``0`` from
    the interval ``[0, x_0)``, so with one component the indicators ``[xi = 0]``
    and ``[xi > 0]`` are both inexpressible.  Those two are precisely what
    decorated inversion needs: inverting a non-strict upper or a strict lower
    predicate gives a bound that is *left*-continuous in the surviving
    variable, which a right-continuous step can represent everywhere except at
    its jumps -- and the jump at ``0`` is the one that carries mass, since the
    magnitude atom there has mass 1 while every other single point has mass 0.
    Splitting the weight lets that one point be carried exactly instead of
    approximated.

    Under Lebesgue ``atom0 = 0`` and the whole mechanism is inert, which is why
    this was invisible until the magnitude oracle was run against mixed
    strict/non-strict predicates.
    """

    coeff: float = 1.0
    weights: Dict[int, Step] = field(default_factory=dict)
    preds: List[PairPredicate] = field(default_factory=list)
    atoms: Dict[int, float] = field(default_factory=dict)

    def clone(self) -> "PairTerm":
        return PairTerm(self.coeff, dict(self.weights), list(self.preds),
                        dict(self.atoms))

    def atom_weight(self, axis: int) -> float:
        return self.atoms.get(axis, 1.0)

    def value_at(self, point: Sequence[float]) -> float:
        """Evaluate the integrand at a point (used only by the oracles).

        At ``point[axis] == 0`` the atom weight is used in place of the step,
        which is what makes the oracle agree with the split representation.
        The oracle gives the anchor its own cell with representative exactly
        ``0.0``, so this branch is hit precisely on the atom.
        """
        v = self.coeff
        if v == 0.0:
            return 0.0
        for axis in set(self.weights) | set(self.atoms):
            t = point[axis]
            if t == 0.0:
                v *= self.atom_weight(axis)
            else:
                w = self.weights.get(axis)
                if w is not None:
                    v *= w(t)
            if v == 0.0:
                return 0.0
        for p in self.preds:
            if not p.holds(point[p.i], point[p.j]):
                return 0.0
        return v

    def mul_weight(self, axis: int, s: Step) -> None:
        """Multiply in a weight that is an ordinary function of the axis.

        A ``Step`` is right-continuous, so ``s(0.0)`` *is* its value at the
        anchor; both components are scaled accordingly.
        """
        cur = self.weights.get(axis)
        self.weights[axis] = s if cur is None else _mul(cur, s)
        self.atoms[axis] = self.atom_weight(axis) * s(0.0)

    def mul_split(self, axis: int, atom_mult: float,
                  step: Optional[Step] = None) -> None:
        """Multiply the atom and the continuous part by different factors.

        This is the operation an ordinary weight cannot perform: ``[xi > 0]`` is
        ``mul_split(axis, 0.0, step_const(1.0))`` and ``[xi == 0]`` is
        ``mul_split(axis, 1.0, step_const(0.0))``.
        """
        self.atoms[axis] = self.atom_weight(axis) * atom_mult
        if step is not None:
            cur = self.weights.get(axis)
            self.weights[axis] = step if cur is None else _mul(cur, step)

    def variables(self) -> List[int]:
        out = set(self.weights) | set(self.atoms)
        for p in self.preds:
            out.add(p.i)
            out.add(p.j)
        return sorted(out)

    def breakpoints(self) -> int:
        n = sum(len(w.xs) for w in self.weights.values())
        n += sum(len(p.f.xs) for p in self.preds)
        return n

    def is_dead(self) -> bool:
        """True only when the term is identically zero.

        An axis whose continuous density vanishes is *not* enough: the anchor
        atom may still carry weight, and under magnitude that is mass 1.  The
        axis is dead only when both components vanish.
        """
        if self.coeff == 0.0:
            return True
        for axis in set(self.weights) | set(self.atoms):
            w = self.weights.get(axis)
            step_dead = w is not None and all(v == 0.0 for v in w.vs)
            if step_dead and self.atom_weight(axis) == 0.0:
                return True
        return False


def _mul(a: Step, b: Step) -> Step:
    xs = tuple(sorted(set(a.xs) | set(b.xs)))
    vs = [a(0.0) * b(0.0)] + [a(x) * b(x) for x in xs]
    return Step(xs, tuple(vs))


def merge_same_orientation(preds: Sequence[PairPredicate]
                           ) -> List[PairPredicate]:
    """Merge predicates of the same pair, side and orientation by min/max.

    Opposite orientations are *not* merged -- their conjunction has a turning
    point and leaves the monotone class.  They are kept as separate predicates,
    which is the repair described in this module's docstring.
    """
    buckets: Dict[Tuple, PairPredicate] = {}
    out: List[PairPredicate] = []
    for p in preds:
        direction = (0 if p.f.is_nondecreasing() else
                     1 if p.f.is_nonincreasing() else 2)
        if direction == 2:                       # already one-turn: keep as is
            out.append(p)
            continue
        key = (p.i, p.j, p.upper, p.strict, direction)
        if key in buckets:
            old = buckets[key]
            merged = old.f.minimum(p.f) if p.upper else old.f.maximum(p.f)
            buckets[key] = PairPredicate(p.i, p.j, merged, p.upper, p.strict)
        else:
            buckets[key] = p
    out.extend(buckets.values())
    return out

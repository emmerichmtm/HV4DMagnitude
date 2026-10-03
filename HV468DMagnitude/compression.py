"""Compression in ``2d`` variables: the report's general compression integrand.

This builds `\\eqref{eq:generalcompression}`,

    N_F(y) = int F(x) prod_i 1[ pi_i^-(y_i) < x_i <= pi_i^+(y_i) ] dmu(x),

as a *single* weighted monotone-pair term in the ``2d`` variables ``(x, y)``,
eliminates the ``d`` old variables with the pair-elimination constructor, and
divides by the cell masses.  The claim being tested is the one the general
theorem rests on: that the compressed function lands back in the same class, so
no latent variable has to persist between compression generations.

Why the term really is in the class
-----------------------------------
The old geometry supplies pair factors among the ``x_i``.  Each grid interval
supplies two more, between ``x_i`` and its own ``y_i``: the predecessor and
successor boundaries are step functions of ``y_i``, so
``pi_i^-(y_i) < x_i <= pi_i^+(y_i)`` is a pair of ordinary pair predicates with
``x_i`` on the constrained side.  Nothing else is needed, and in particular no
constraint ever links two different ``y`` variables.

The cell convention, and why it is not the report's
----------------------------------------------------
The report writes the cell as ``(pi^-, pi^+]``, left-open and right-closed.
That is the natural way to write it in prose and the wrong way to represent it
here: a family of ``(a, b]`` cells is left-continuous in ``y``, and a
right-continuous ``Step`` cannot express which cell ``y`` falls in.  The
mismatch is harmless for ``y > 0`` -- it misplaces only the boundary points,
which carry no mass -- but at ``y = 0`` the magnitude atom has mass 1 and the
answer would simply be wrong.

So cells here are ``[g_m, g_{m+1})``, right-open, with the top boundary placed
*beyond* the domain ceiling so that every cell has the same decoration and the
last one still closes at ``hi`` by the ground set rather than by a special case.
The anchor then sits inside the first cell, where it belongs, and both boundary
maps are right-continuous step functions of ``y``.

This is a representation choice, not a change of content: any partition into
intervals gives a valid conditioning, and this one is the partition that the
state class can actually carry.
"""

from __future__ import annotations

import itertools
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core import INF, Measure1D, Step
from pair_elimination import Counters, eliminate_variable
from pair_state import PairPredicate, PairTerm

__all__ = ["CompressionGrid", "build_compression_term", "compress",
           "numerator_at", "denominator_at", "conditional_at", "cell_indices",
           "reciprocal_mass_step", "divide_by_cell_mass", "merge_terms",
           "reindex_to_physical"]


@dataclass
class CompressionGrid:
    """A future-complete grid on one axis.

    ``boundaries`` is ``[0 = g_0, g_1, ..., g_k]`` with ``g_k`` strictly above
    the domain ceiling, so the cells ``[g_m, g_{m+1})`` cover ``[0, hi]`` with a
    single decoration throughout.
    """

    boundaries: List[float]

    def __post_init__(self):
        b = self.boundaries
        if b[0] != 0.0:
            raise ValueError("the first boundary must be the anchor 0")
        if any(q <= p for p, q in zip(b, b[1:])):
            raise ValueError("boundaries must be strictly increasing")

    @property
    def interior(self) -> Tuple[float, ...]:
        return tuple(self.boundaries[1:-1])

    def pred_step(self) -> Step:
        """``pi^-(y)``: the lower boundary of the cell containing ``y``."""
        return Step(self.interior, tuple(self.boundaries[:-1]))

    def succ_step(self) -> Step:
        """``pi^+(y)``: the upper boundary of the cell containing ``y``."""
        return Step(self.interior, tuple(self.boundaries[1:]))

    def n_cells(self) -> int:
        return len(self.boundaries) - 1

    def cell_of(self, y: float) -> int:
        for m in range(self.n_cells()):
            if self.boundaries[m] <= y < self.boundaries[m + 1]:
                return m
        raise ValueError(f"{y} is outside the grid")

    def cell_bounds(self, m: int, hi: float) -> Tuple[float, float]:
        """The cell clipped to the ground set, as the oracle wants it."""
        return (self.boundaries[m], min(self.boundaries[m + 1], hi))

    def representatives(self, m: int, hi: float) -> List[float]:
        """Points of cell ``m`` to evaluate at: its floor and its midpoint.

        The floor matters: for ``m = 0`` it is the anchor ``0``, which is the one
        point where getting the cell convention wrong changes the answer.
        """
        lo, up = self.cell_bounds(m, hi)
        out = [lo]
        mid = 0.5 * (lo + up)
        if mid > lo:
            out.append(mid)
        return out

    def cell_mass(self, m: int, mu: Measure1D, hi: float) -> float:
        lo, up = self.cell_bounds(m, hi)
        mass = mu.density * (up - lo)
        if lo == 0.0:
            mass += mu.atom0          # the anchor lives in the first cell
        return mass


def build_compression_term(f_term: PairTerm, d: int,
                           grids: Sequence[CompressionGrid]) -> PairTerm:
    """The integrand of `\\eqref{eq:generalcompression}` as one term in 2d vars.

    ``x`` occupies axes ``0 .. d-1`` and ``y`` axes ``d .. 2d-1``.  Each axis
    contributes two predicates with ``x_i`` on the *constrained* side, so no
    inversion is needed to read them as bounds on ``x_i``:

    * ``x_i >= pi_i^-(y_i)`` -- lower, non-strict;
    * ``x_i <  pi_i^+(y_i)`` -- upper, strict.

    Both are the right-continuous decorations, which is what the cell convention
    above was chosen to achieve.
    """
    term = f_term.clone()
    for i in range(d):
        y = d + i
        g = grids[i]
        term.preds.append(PairPredicate(y, i, g.pred_step(),
                                        upper=False, strict=False))
        term.preds.append(PairPredicate(y, i, g.succ_step(),
                                        upper=True, strict=True))
    return term


def compress(f_term: PairTerm, d: int,
             measures: Dict[int, Measure1D],
             domains: Dict[int, Tuple[float, float]],
             grids: Sequence[CompressionGrid],
             event_grids: Dict[int, Sequence[float]],
             counters: Optional[Counters] = None,
             term_cap: Optional[int] = None) -> Optional[List[PairTerm]]:
    """Eliminate the ``d`` old variables, leaving a signed sum of terms in ``y``.

    Returns ``None`` if the constructor abstained or the state exceeded
    ``term_cap``, and ``[]`` when the numerator is identically zero.
    """
    terms = [build_compression_term(f_term, d, grids)]
    for x in range(d):
        nxt: List[PairTerm] = []
        for t in terms:
            got = eliminate_variable(t, x, measures[x], domains[x],
                                     event_grids[x], counters)
            if got is None:
                return None
            nxt.extend(got)
        terms = nxt
        if counters:
            counters.peak_state = max(counters.peak_state, len(terms))
        if term_cap is not None and len(terms) > term_cap:
            return None
        if not terms:
            return []
    return terms


def reciprocal_mass_step(g: CompressionGrid, mu: Measure1D,
                         hi: float) -> Step:
    """``1 / mu(cell(y))`` as a unary step density in ``y``.

    The report's claim that dividing by the cell mass "only changes unary step
    densities" holds exactly, and this is why: the mass is constant on each
    cell, so as a function of ``y`` it is a step function with breakpoints at
    the interior boundaries -- and with the right-open cell convention it is
    right-continuous, hence representable.  Note the anchor is inside the first
    cell, so under magnitude the first cell's mass carries the atom and every
    ``y`` in that cell, including ``y = 0``, divides by the same number.
    """
    vals = []
    for m in range(g.n_cells()):
        mass = g.cell_mass(m, mu, hi)
        vals.append(0.0 if mass == 0.0 else 1.0 / mass)
    return Step(g.interior, tuple(vals))


def divide_by_cell_mass(terms: Sequence[PairTerm], d: int,
                        grids: Sequence[CompressionGrid],
                        measures: Dict[int, Measure1D],
                        domains: Dict[int, Tuple[float, float]]
                        ) -> List[PairTerm]:
    """Turn the numerator into the conditional expectation, as terms.

    This is the step that makes the output usable recursively: the result is
    again a weighted monotone-pair state, in the new physical variables, with no
    latent variable carried over -- which is the substance of the report's
    closure claim.  Only unary densities change; no predicate is touched.
    """
    out = []
    for t in terms:
        s = t.clone()
        for i in range(d):
            s.mul_weight(d + i, reciprocal_mass_step(grids[i], measures[i],
                                                     domains[i][1]))
        out.append(s)
    return out


def merge_terms(terms: Sequence[PairTerm]) -> List[PairTerm]:
    """Collapse terms that are identical apart from their coefficient.

    This is the *canonical* merge, the only one available without reasoning
    about the geometry: terms whose weights and predicates agree exactly have
    their coefficients added, and terms whose coefficients cancel to zero are
    dropped.  `\\ref{prop:complexity}` is explicit that nothing proves a signed
    sum of distinct primitives collapses to a constant number, so measuring how
    much this recovers is the honest way to find out what the representation
    actually costs.
    """
    buckets: Dict[Tuple, PairTerm] = {}
    order: List[Tuple] = []
    for t in terms:
        key = (tuple(sorted((a, w.xs, w.vs) for a, w in t.weights.items())),
               tuple(sorted((a, v) for a, v in t.atoms.items())),
               tuple(sorted(p.key() for p in t.preds)))
        if key in buckets:
            buckets[key].coeff += t.coeff
        else:
            s = t.clone()
            buckets[key] = s
            order.append(key)
    return [buckets[k] for k in order if buckets[k].coeff != 0.0]


def reindex_to_physical(terms: Sequence[PairTerm], d: int) -> List[PairTerm]:
    """Relabel the ``y`` axes ``d..2d-1`` back to the physical axes ``0..d-1``.

    Needed to feed one generation's output into the next.  That this is a pure
    relabelling -- no residual ``x`` axis, no auxiliary variable -- is exactly
    the property the report needs for repeated compression.
    """
    out = []
    for t in terms:
        s = PairTerm(t.coeff, {}, [], {})
        for axis, w in t.weights.items():
            if axis < d:
                raise ValueError(f"axis {axis} is not a y variable; "
                                 "the x variables were not all eliminated")
            s.weights[axis - d] = w
        for axis, a in t.atoms.items():
            s.atoms[axis - d] = a
        for p in t.preds:
            if p.i < d or p.j < d:
                raise ValueError("a predicate still mentions an x variable")
            s.preds.append(PairPredicate(p.i - d, p.j - d, p.f, p.upper,
                                         p.strict))
        out.append(s)
    return out


def cell_indices(d: int, grids: Sequence[CompressionGrid]):
    """Every cell of the product grid, as a tuple of per-axis cell indices."""
    return itertools.product(*[range(g.n_cells()) for g in grids])


def numerator_at(terms: Sequence[PairTerm], d: int,
                 y_values: Sequence[float]) -> float:
    """Evaluate the compressed numerator at one ``y``."""
    point = [0.0] * (2 * d)
    for i, v in enumerate(y_values):
        point[d + i] = v
    return sum(t.value_at(point) for t in terms)


def denominator_at(d: int, grids: Sequence[CompressionGrid],
                   measures: Dict[int, Measure1D],
                   domains: Dict[int, Tuple[float, float]],
                   cells: Sequence[int]) -> float:
    out = 1.0
    for i, m in enumerate(cells):
        out *= grids[i].cell_mass(m, measures[i], domains[i][1])
    return out


def conditional_at(terms: Sequence[PairTerm], d: int,
                   grids: Sequence[CompressionGrid],
                   measures: Dict[int, Measure1D],
                   domains: Dict[int, Tuple[float, float]],
                   cells: Sequence[int],
                   y_values: Sequence[float]) -> float:
    """``E[F | grid]`` on this cell: the numerator over the cell mass."""
    den = denominator_at(d, grids, measures, domains, cells)
    if den == 0.0:
        return 0.0
    return numerator_at(terms, d, y_values) / den

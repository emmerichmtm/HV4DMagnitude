"""Core objects for the product-measure / magnitude 4-D prefix experiment.

Everything here is deliberately simple and slow where that buys clarity: this
module exists so that the *claims* in the accompanying report can be checked
against exact brute force, not so that anything runs fast.

Conventions, stated once and relied on everywhere
-------------------------------------------------

**Ground set.** All coordinates live in ``[0, inf)``.  The anchor ``0`` is a
real point of the space, not a formal lower limit: under the magnitude measure
it carries an atom, so "the prefix below 0" and "the prefix up to and including
0" are different sets with different masses.  Nothing here may silently
translate a coordinate.

**One-dimensional measure.**  Every measure used is of the form

    mu = atom0 * delta_0 + density * lambda

on ``[0, inf)``.  ``Lebesgue`` is ``atom0=0, density=1`` and ``Magnitude`` is
``atom0=1, density=1/2``; the latter is the one-dimensional factor of
``(delta_0 + lambda/2)^{⊗d}``, whose product over d coordinates gives the
ell_1 magnitude of an anchored downward-closed set.  Points ``t > 0`` have
measure zero, so the *only* place where a closed/open endpoint distinction can
change a mass is at ``0``.

**Step functions.**  ``Step`` is piecewise constant and **right-continuous**:
with breakpoints ``x_0 < x_1 < ... < x_{k-1}`` (all strictly positive) and
values ``v_0, ..., v_k``,

    f(t) = v_0           for 0 <= t < x_0
    f(t) = v_{j+1}       for x_j <= t < x_{j+1}
    f(t) = v_k           for t >= x_{k-1}

Right-continuity is forced by the geometry, not chosen for convenience: the
staircase bounding a union of grounded quadrants ``{x_i >= a_t, x_j >= b_t}``
is ``f(x) = min{b_t : a_t <= x}``, which drops exactly *at* ``a_t``.

**Staircase constraints are strict.**  The region left uncovered by that union
is ``{x_j < f(x_i)}``, with a strict inequality.  The report writes
``[x_j <= f_ij(x_i)]``.  Under Lebesgue measure the two agree (they differ on a
null set), but under the magnitude measure they differ whenever ``f`` can take
the value 0, because the atom at ``x_j = 0`` then has positive mass.  We
implement the strict form and provide ``closed=True`` to reproduce the
report's form for comparison.

**Generalized inverse.**  For nonincreasing right-continuous ``f`` the set
``{x : f(x) >= u}`` is a prefix of the form ``[0, theta)`` -- half open, with
``theta`` possibly ``+inf`` and possibly ``0`` (meaning empty).  We return
``theta`` together with the fact that it is an *exclusive* bound, because
treating it as inclusive silently adds the atom at 0 when the set is empty.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from dataclasses import dataclass
from typing import Callable, List, Sequence, Tuple

INF = math.inf

__all__ = [
    "Measure1D", "LEBESGUE", "MAGNITUDE", "Step",
    "inv_ge", "inv_le", "brute_integral", "cells_1d",
]


# --------------------------------------------------------------------------- #
# One-dimensional product-measure factors
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Measure1D:
    """``atom0 * delta_0 + density * lambda`` on ``[0, inf)``."""

    atom0: float
    density: float
    name: str = ""

    def prefix_closed(self, t: float) -> float:
        """mu([0, t]).  Zero for t < 0."""
        if t < 0.0:
            return 0.0
        return self.atom0 + self.density * t

    def prefix_open(self, t: float) -> float:
        """mu([0, t)).  Note this is 0 at t = 0, not ``atom0``."""
        if t <= 0.0:
            return 0.0
        return self.atom0 + self.density * t

    def interval(self, lo: float, hi: float, lo_closed=True, hi_closed=True) -> float:
        """Mass of the interval between ``lo`` and ``hi``."""
        if hi < lo:
            return 0.0
        mass = self.density * (hi - lo)
        if lo == 0.0 and lo_closed:
            mass += self.atom0
        if lo == hi:
            return self.atom0 if (lo == 0.0 and lo_closed and hi_closed) else 0.0
        return mass


LEBESGUE = Measure1D(atom0=0.0, density=1.0, name="lebesgue")
MAGNITUDE = Measure1D(atom0=1.0, density=0.5, name="magnitude")


# --------------------------------------------------------------------------- #
# Right-continuous step functions on [0, inf)
# --------------------------------------------------------------------------- #


class Step:
    """Right-continuous piecewise-constant function on ``[0, inf)``.

    ``xs`` are strictly increasing positive breakpoints; ``vs`` has one more
    entry than ``xs``.  ``vs[0]`` is the value on ``[0, xs[0])`` and ``vs[j+1]``
    the value on ``[xs[j], xs[j+1])``.
    """

    __slots__ = ("xs", "vs")

    def __init__(self, xs: Sequence[float], vs: Sequence[float]) -> None:
        if len(vs) != len(xs) + 1:
            raise ValueError("len(vs) must be len(xs) + 1")
        if any(b <= 0.0 for b in xs):
            raise ValueError("breakpoints must be strictly positive")
        if any(b >= c for b, c in zip(xs, xs[1:])):
            raise ValueError("breakpoints must be strictly increasing")
        # drop breakpoints that do not change the value
        cxs: List[float] = []
        cvs: List[float] = [vs[0]]
        for b, v in zip(xs, vs[1:]):
            if v != cvs[-1]:
                cxs.append(b)
                cvs.append(v)
        self.xs = tuple(cxs)
        self.vs = tuple(cvs)

    @staticmethod
    def const(value: float) -> "Step":
        return Step((), (value,))

    def __call__(self, t: float) -> float:
        """Value at ``t``; right-continuous, so f(x_j) is the value to the right."""
        return self.vs[bisect_right(self.xs, t)]

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Step(xs={self.xs}, vs={self.vs})"

    def breakpoints(self) -> Tuple[float, ...]:
        return self.xs

    def is_nonincreasing(self) -> bool:
        return all(a >= b for a, b in zip(self.vs, self.vs[1:]))

    def is_nondecreasing(self) -> bool:
        return all(a <= b for a, b in zip(self.vs, self.vs[1:]))

    def minimum(self, other: "Step") -> "Step":
        return self._combine(other, min)

    def maximum(self, other: "Step") -> "Step":
        return self._combine(other, max)

    def _combine(self, other: "Step", op) -> "Step":
        xs = tuple(sorted(set(self.xs) | set(other.xs)))
        vs = [op(self(0.0), other(0.0))]
        for b in xs:
            vs.append(op(self(b), other(b)))
        return Step(xs, vs)


def staircase_from_quadrants(points: Sequence[Tuple[float, float]]) -> Step:
    """Boundary of the region left uncovered by grounded quadrants.

    Each ``(a, b)`` is the vertex of ``{x_i >= a, x_j >= b}``.  The uncovered
    region is ``{x_j < f(x_i)}`` with ``f(x) = min{b_t : a_t <= x}``, which is
    nonincreasing and right-continuous.  ``f = +inf`` where no quadrant applies.
    """
    best: dict = {}
    for a, b in points:
        best[a] = min(best.get(a, INF), b)
    xs = sorted(k for k in best if k > 0.0)
    running = min([best[k] for k in best if k <= 0.0], default=INF)
    vs = [running]
    for b in xs:
        running = min(running, best[b])
        vs.append(running)
    return Step(tuple(xs), tuple(vs))


# --------------------------------------------------------------------------- #
# Generalized inverses
# --------------------------------------------------------------------------- #


def inv_ge(f: Step, u: float) -> float:
    """``theta`` with ``{x >= 0 : f(x) >= u} == [0, theta)``, for nonincreasing f.

    The bound is **exclusive**.  ``theta == 0`` means the set is empty (note
    that this is *not* the same as the set ``{0}``), and ``theta == inf`` means
    the constraint is vacuous.  Callers must therefore use ``prefix_open`` and
    never ``prefix_closed`` on this value, or the magnitude atom at 0 is added
    to an empty set.
    """
    if not f.is_nonincreasing():
        raise ValueError("inv_ge expects a nonincreasing step function")
    if f.vs[0] < u:
        return 0.0                      # empty: even at x = 0 the value is too small
    # last index whose value is >= u
    last = 0
    for j, v in enumerate(f.vs):
        if v >= u:
            last = j
        else:
            break
    if last == len(f.vs) - 1:
        return INF
    return f.xs[last]                   # value drops below u exactly at this point


def inv_le(f: Step, u: float) -> float:
    """``theta`` with ``{x >= 0 : f(x) <= u} == [0, theta)``, for nondecreasing f.

    Same exclusive-bound convention as :func:`inv_ge`.
    """
    if not f.is_nondecreasing():
        raise ValueError("inv_le expects a nondecreasing step function")
    if f.vs[0] > u:
        return 0.0
    last = 0
    for j, v in enumerate(f.vs):
        if v <= u:
            last = j
        else:
            break
    if last == len(f.vs) - 1:
        return INF
    return f.xs[last]


# --------------------------------------------------------------------------- #
# Exact brute force over a finite grid
# --------------------------------------------------------------------------- #


def cells_1d(grid: Sequence[float], mu: Measure1D,
             upper: float) -> List[Tuple[float, float]]:
    """Decompose ``[0, upper]`` into (representative point, mass) pairs.

    The decomposition is exact for any function that is constant on the open
    intervals between consecutive grid points: the atom at 0 becomes its own
    cell, and every open interval contributes its representative midpoint.
    Points ``t > 0`` carry no mass and are therefore not enumerated -- which is
    exactly why a step function's value *at* a breakpoint can only matter at 0.
    """
    pts = sorted({0.0} | {g for g in grid if 0.0 < g <= upper} | {upper})
    pts = [p for p in pts if p <= upper]
    cells: List[Tuple[float, float]] = []
    if mu.atom0 != 0.0:
        cells.append((0.0, mu.atom0))
    for lo, hi in zip(pts, pts[1:]):
        if hi > lo:
            cells.append(((lo + hi) / 2.0, mu.density * (hi - lo)))
    return cells


def brute_integral(integrand: Callable[..., float],
                   grids: Sequence[Sequence[float]],
                   measures: Sequence[Measure1D],
                   uppers: Sequence[float]) -> float:
    """Exact integral of a grid-aligned ``integrand`` over a prefix box.

    ``integrand`` takes one argument per coordinate.  It must be constant on
    each open cell of the product grid; its value at the anchor 0 is sampled
    separately, so atoms are handled correctly.
    """
    per_axis = [cells_1d(g, m, u) for g, m, u in zip(grids, measures, uppers)]
    total = 0.0

    def walk(axis: int, point: List[float], weight: float) -> float:
        if weight == 0.0:
            return 0.0
        if axis == len(per_axis):
            value = integrand(*point)
            return weight * value if value else 0.0
        acc = 0.0
        for rep, mass in per_axis[axis]:
            point.append(rep)
            acc += walk(axis + 1, point, weight * mass)
            point.pop()
        return acc

    total = walk(0, [], 1.0)
    return total

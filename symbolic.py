"""The four-prefix normal form emitted as *data* rather than as an evaluator.

`prefix4.Prefix4` computes K(a,b,c,d) by walking the target plane.  That is
enough to check the mathematics but it is not the object the report claims
exists: a dimension-bounded list of primitives

    J(C, D) = int_0^C alpha(u) * Gamma(min(D, h_up(u), h_down(u))) d nu(u),
    Gamma(t) = int_0^t beta(s) d omega(s),

carrying explicit unary data, from which any prefix query is answered by
arithmetic on precomputed cumulative records.  This module builds that list.

How the boundaries are obtained
-------------------------------
The region boundaries are monotone step functions built from generalized
inverses, their compositions and pointwise minima.  Rather than deriving each
one in closed form -- the composition rules are where sign and orientation
errors breed -- each boundary is computed *exactly* by evaluating its defining
predicate at the finitely many candidate breakpoints.  This is exact, not
fitted: every function in sight is piecewise constant with breakpoints in that
candidate set, so a predicate that is monotone in the latent plane is pinned
down by its values there.  The construction then asserts the monotonicity it
expects, so an orientation error fails loudly instead of silently.

Region algebra
--------------
A region is a conjunction of atoms, each possibly negated:

    ('u_le', s)   u <= s            ('v_le', t)   v <= t
    ('inc', phi)  v <= phi(u), phi nondecreasing
    ('dec', phi)  v <= phi(u), phi nonincreasing

Negations expand by ``[not A] = 1 - [A]``, so a region becomes a signed sum of
pure conjunctions; a pure conjunction intersects to
``{u <= C} x {v <= min(D, h_up(u), h_down(u))}``, which is exactly the shape of
a primitive.  Two conditions of the same orientation merge by pointwise
minimum, which is the closure property tested in Phase 5.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from core import INF, Measure1D, Step
from identities import _cum, prefix_bound_gt
from prefix4 import Bound, Prefix4, SixStaircases, cap, mass


def min_bound(x: Bound, y: Bound) -> Bound:
    """Tighter of two prefix bounds; an exclusive bound wins a tie.

    Needed because the construction mixes sources with different endpoint
    semantics: query arguments are inclusive, while every bound derived from a
    generalized inverse is exclusive.  Minimising the raw numbers would put the
    anchor atom on the wrong side whenever a bound is 0.
    """
    if x[0] < y[0]:
        return x
    if y[0] < x[0]:
        return y
    return (x[0], x[1] and y[1])

__all__ = ["PrefixPrimitive", "PrefixState4D", "build_state"]


# --------------------------------------------------------------------------- #
# Primitive
# --------------------------------------------------------------------------- #


@dataclass
class PrefixPrimitive:
    """One term of the normal form, as data.

    ``c_of`` and ``d_of`` map the prefix arguments ``(a,b,c,d)`` to the limits
    ``C`` and ``D``.  They are minima of the arguments with monotone unary
    images of them, so they introduce no new coordinates -- the property the
    shared-record bound relies on.
    """

    sign: float
    alpha: Step                      # weight in the latent u
    beta: Step                       # weight in the latent v
    h_up: Optional[Step]             # nondecreasing bound on v given u
    h_down: Optional[Step]           # nonincreasing bound on v given u
    u_bound: Bound                   # cap on the latent u, from the region
    v_bound: Bound                   # cap on the latent v, from the region
    label: Tuple = ()                # provenance, for instrumentation

    def limits(self, a: float, b: float, c: float, d: float) -> Tuple[Bound, Bound]:
        """``(C, D)`` for a query: the region caps met with the query box.

        Both are minima of query arguments with monotone unary images of them,
        so no new coordinate is introduced -- the property the shared-record
        bound needs.
        """
        return (min_bound((c, True), self.u_bound),
                min_bound((d, True), self.v_bound))

    def bound_v(self, u: float, D: Bound) -> Bound:
        """The v-prefix at this ``u``, with its endpoint semantics.

        Boundaries produced by ``boundary_step`` are exclusive (the set is
        ``[0, psi(u))``) while the query argument ``D`` is inclusive, so the
        two cannot simply be minimised together: at a bound of 0 the magnitude
        atom belongs to one and not the other.
        """
        t_open = INF
        if self.h_up is not None:
            t_open = min(t_open, self.h_up(u))
        if self.h_down is not None:
            t_open = min(t_open, self.h_down(u))
        return min_bound(D, (t_open, False))


class PrefixState4D:
    """A list of primitives plus the records needed to evaluate them."""

    def __init__(self, prims: Sequence[PrefixPrimitive],
                 mu_u: Measure1D, mu_v: Measure1D,
                 grid_u: Sequence[float], grid_v: Sequence[float]):
        self.prims = list(prims)
        self.mu_u, self.mu_v = mu_u, mu_v
        self.grid_u = sorted(set(grid_u))
        self.grid_v = sorted(set(grid_v))
        self.counters: Dict[str, int] = {
            "primitives": len(self.prims),
            "records": 0,
            "queries": 0,
        }
        self._gamma_cache: Dict[int, Dict[float, float]] = {}
        self.prune()
        self._register_records()

    def prune(self) -> None:
        """Drop structurally empty primitives and merge identical ones.

        The signed expansion of the region algebra emits many terms that cannot
        contribute: a weight that is identically zero, or an empty prefix.
        Merging the survivors by structure keeps the reported primitive count
        honest -- the raw expansion count overstates the state.
        """
        merged: Dict[Tuple, PrefixPrimitive] = {}
        for prim in self.prims:
            if all(v == 0.0 for v in prim.alpha.vs):
                continue
            if all(v == 0.0 for v in prim.beta.vs):
                continue
            if prim.u_bound == (0.0, False) or prim.v_bound == (0.0, False):
                continue
            key = ((prim.alpha.xs, prim.alpha.vs),
                   (prim.beta.xs, prim.beta.vs),
                   None if prim.h_up is None else (prim.h_up.xs, prim.h_up.vs),
                   None if prim.h_down is None else (prim.h_down.xs, prim.h_down.vs),
                   prim.u_bound, prim.v_bound)
            if key in merged:
                merged[key].sign += prim.sign
            else:
                merged[key] = PrefixPrimitive(
                    prim.sign, prim.alpha, prim.beta, prim.h_up, prim.h_down,
                    prim.u_bound, prim.v_bound, prim.label)
        self.prims = [q for q in merged.values() if q.sign != 0.0]
        self.counters["primitives"] = len(self.prims)

    def _register_records(self) -> None:
        """Count the one-dimensional records the state stores.

        One record per primitive per grid coordinate on each axis: the step
        value together with the cumulative at the left endpoint.  This is the
        quantity the report bounds by O(m).
        """
        per_prim = len(self.grid_u) + len(self.grid_v)
        self.counters["records"] = len(self.prims) * per_prim

    def Gamma(self, prim: PrefixPrimitive, bound: Bound) -> float:
        """``int beta d omega`` over the v-prefix, honouring its endpoint."""
        t, closed = bound
        if t < 0.0:
            return 0.0
        return _cum(prim.beta.__call__, self._v_events(prim), self.mu_v, t,
                    closed=closed)

    def _v_events(self, prim: PrefixPrimitive) -> List[float]:
        return sorted(set(self.grid_v) | set(prim.beta.xs))

    def evaluate(self, a: float, b: float, c: float, d: float) -> float:
        """``K(a,b,c,d)`` from the primitives alone."""
        self.counters["queries"] += 1
        total = 0.0
        for prim in self.prims:
            C, D = prim.limits(a, b, c, d)
            if C[0] < 0.0 or D[0] < 0.0:
                continue
            total += prim.sign * self._integrate(prim, C, D)
        return total

    def _integrate(self, prim: PrefixPrimitive, C: Bound, D: Bound) -> float:
        """``int_{[0,C]} alpha(u) Gamma(min(D, h(u))) d nu(u)``."""
        def weight(u: float) -> float:
            av = prim.alpha(u)
            if av == 0.0:
                return 0.0
            return av * self.Gamma(prim, prim.bound_v(u, D))
        return _cum(weight, self._u_events(prim), self.mu_u, C[0],
                    closed=C[1])

    def _u_events(self, prim: PrefixPrimitive) -> List[float]:
        events = set(self.grid_u) | set(prim.alpha.xs)
        if prim.h_up is not None:
            events |= set(prim.h_up.xs)
        if prim.h_down is not None:
            events |= set(prim.h_down.xs)
        return sorted(events)


# --------------------------------------------------------------------------- #
# Exact monotone boundaries from their defining predicates
# --------------------------------------------------------------------------- #


def boundary_step(pred: Callable[[float, float], bool],
                  u_points: Sequence[float], v_points: Sequence[float],
                  expect: str) -> Step:
    """Exact ``psi`` with ``{v : pred(u, v)} = [0, psi(u))`` on the candidates.

    ``pred`` must be downward closed in ``v`` for each ``u``; ``expect`` is the
    orientation asserted of the result (``'inc'`` or ``'dec'``), so a wrong
    orientation surfaces here rather than as a wrong number later.
    """
    us = sorted(set(u_points))
    vs = sorted(set(v for v in v_points if v > 0.0))
    values: List[float] = []
    for u in us:
        bound = INF
        for v in vs:
            if not pred(u, v):
                bound = v
                break
        if not pred(u, 0.0):
            bound = 0.0
        values.append(bound)
    xs, vls = [], [values[0]]
    for u, val in zip(us[1:], values[1:]):
        if val != vls[-1] and u > 0.0:
            xs.append(u)
            vls.append(val)
    step = Step(tuple(xs), tuple(vls))
    if expect == "inc" and not step.is_nondecreasing():
        raise AssertionError(f"expected a nondecreasing boundary, got {step}")
    if expect == "dec" and not step.is_nonincreasing():
        raise AssertionError(f"expected a nonincreasing boundary, got {step}")
    return step


def cap_u(pred: Callable[[float], bool], u_points: Sequence[float]) -> Bound:
    """``{u : pred(u)} = [0, s)`` for a prefix predicate, as an EXCLUSIVE bound."""
    us = sorted(set(u for u in u_points if u > 0.0))
    if not pred(0.0):
        return (0.0, False)
    for u in us:
        if not pred(u):
            return (u, False)
    return (INF, False)


# --------------------------------------------------------------------------- #
# Region algebra
# --------------------------------------------------------------------------- #

Atom = Tuple[str, object]


def expand(atoms: Sequence[Tuple[Atom, bool]]
           ) -> List[Tuple[float, List[Atom]]]:
    """Expand negations into a signed sum of pure conjunctions."""
    out: List[Tuple[float, List[Atom]]] = [(1.0, [])]
    for atom, negated in atoms:
        nxt: List[Tuple[float, List[Atom]]] = []
        for sign, conj in out:
            if negated:
                nxt.append((sign, list(conj)))             # the "1" of 1 - [A]
                nxt.append((-sign, list(conj) + [atom]))   # minus [A]
            else:
                nxt.append((sign, list(conj) + [atom]))
        out = nxt
    return out


def conjunction_shape(conj: Sequence[Atom]
                      ) -> Tuple[Bound, Bound, Optional[Step], Optional[Step]]:
    """Intersect a pure conjunction into ``(u_bound, v_bound, h_up, h_down)``."""
    u_cap: Bound = (INF, False)
    v_cap: Bound = (INF, False)
    h_up: Optional[Step] = None
    h_down: Optional[Step] = None
    for kind, payload in conj:
        if kind == "u_le":
            u_cap = min_bound(u_cap, payload)
        elif kind == "v_le":
            v_cap = min_bound(v_cap, payload)
        elif kind == "inc":
            h_up = payload if h_up is None else h_up.minimum(payload)
        elif kind == "dec":
            h_down = payload if h_down is None else h_down.minimum(payload)
        else:
            raise ValueError(f"unknown atom {kind!r}")
    return u_cap, v_cap, h_up, h_down


# --------------------------------------------------------------------------- #
# The constructor
# --------------------------------------------------------------------------- #


def build_state(st: SixStaircases, measures: Sequence[Measure1D],
                grids: Sequence[Sequence[float]],
                a: float, b: float) -> PrefixState4D:
    """Emit the primitives for the six-staircase state at source caps ``(a, b)``.

    ``a`` and ``b`` are fixed at build time because they select which regions
    are rectangles; ``c`` and ``d`` remain free query arguments.  This matches
    the report, where the source caps enter as anchored target rectangles.
    """
    p4 = Prefix4(st, measures, grids)
    mu3, mu4 = measures[2], measures[3]

    # Candidate latent coordinates: the values of the staircases that define
    # rho/sigma, plus the grid.  Every boundary is constant between these.
    u_pts = sorted({0.0} | {v for v in st.f13.vs if v != INF}
                   | {v for v in st.f23.vs if v != INF}
                   | set(grids[2]) | set(st.f34.xs))
    v_pts = sorted({0.0} | {v for v in st.f14.vs if v != INF}
                   | {v for v in st.f24.vs if v != INF}
                   | set(grids[3]) | {v for v in st.f34.vs if v != INF})

    def A_of(u, v):
        return cap(a, min(p4.rho3(u), p4.rho4(v)))

    def B_of(u, v):
        return cap(b, min(p4.sig3(u), p4.sig4(v)))

    # --- fixed boundary curves ------------------------------------------- #
    # which cap is active for the x1 prefix
    sel_a_u = cap_u(lambda u: p4.rho3(u) >= a, u_pts)      # u-side of the rectangle
    sel_a_v = cap_u(lambda v: p4.rho4(v) >= a, v_pts)
    sel_b_u = cap_u(lambda u: p4.sig3(u) >= b, u_pts)
    sel_b_v = cap_u(lambda v: p4.sig4(v) >= b, v_pts)
    # rho3(u) <= rho4(v): nondecreasing frontier in v
    phi_a = boundary_step(lambda u, v: p4.rho3(u) <= p4.rho4(v),
                          u_pts, v_pts, expect="inc")
    phi_b = boundary_step(lambda u, v: p4.sig3(u) <= p4.sig4(v),
                          u_pts, v_pts, expect="inc")
    # branch 2 (the cap bites) is downward closed in the target plane
    psi = boundary_step(lambda u, v: not p4._branch1(A_of(u, v), B_of(u, v)),
                        u_pts, v_pts, expect="dec")
    # the target staircase g = f34 : v < g(u)
    g_bound = st.f34

    prims: List[PrefixPrimitive] = []

    for sA in ("a", "u", "v"):
        for sB in ("b", "u", "v"):
            for branch in (1, 2):
                atoms: List[Tuple[Atom, bool]] = []
                # target box and the g staircase are common to every region
                atoms.append((("dec", g_bound), False))

                # --- selector for the x1 cap --------------------------------
                if sA == "a":
                    atoms.append((("u_le", sel_a_u), False))
                    atoms.append((("v_le", sel_a_v), False))
                elif sA == "u":
                    atoms.append((("u_le", sel_a_u), True))     # rho3 < a
                    atoms.append((("inc", phi_a), False))       # rho3 <= rho4
                else:
                    atoms.append((("v_le", sel_a_v), True))     # rho4 < a
                    atoms.append((("inc", phi_a), True))        # rho4 <  rho3

                # --- selector for the x2 cap --------------------------------
                if sB == "b":
                    atoms.append((("u_le", sel_b_u), False))
                    atoms.append((("v_le", sel_b_v), False))
                elif sB == "u":
                    atoms.append((("u_le", sel_b_u), True))
                    atoms.append((("inc", phi_b), False))
                else:
                    atoms.append((("v_le", sel_b_v), True))
                    atoms.append((("inc", phi_b), True))

                # --- branch -------------------------------------------------
                atoms.append((("dec", psi), branch == 1))

                for alpha, beta in _integrand_factors(p4, sA, sB, branch, a, b,
                                                      u_pts, v_pts):
                    for sign, conj in expand(atoms):
                        u_cap, v_cap, h_up, h_down = conjunction_shape(conj)
                        prims.append(PrefixPrimitive(
                            sign=sign, alpha=alpha, beta=beta,
                            h_up=h_up, h_down=h_down,
                            u_bound=u_cap, v_bound=v_cap,
                            label=(sA, sB, branch)))

    return PrefixState4D(prims, mu3, mu4, u_pts, v_pts)


def _integrand_factors(p4: Prefix4, sA: str, sB: str, branch: int,
                       a: float, b: float,
                       u_pts: Sequence[float], v_pts: Sequence[float]
                       ) -> List[Tuple[Step, Step]]:
    """The separable factors of the source mass on one region.

    Branch 1 gives a single product ``X(A') * Y(B')``; branch 2 gives the sum
    ``P(A') + Q(B')`` and therefore two primitives.  Each factor depends on at
    most one latent coordinate -- that is the separability the normal form
    needs, and it is what makes these representable as unary ``Step`` data.
    """
    one = Step.const(1.0)

    def as_u(fn) -> Step:
        xs = [u for u in u_pts if u > 0.0]
        return Step(tuple(xs), tuple([fn(0.0)] + [fn(u) for u in xs]))

    def as_v(fn) -> Step:
        xs = [v for v in v_pts if v > 0.0]
        return Step(tuple(xs), tuple([fn(0.0)] + [fn(v) for v in xs]))

    def A_bound_u(u):
        return cap(a, p4.rho3(u))

    def A_bound_v(v):
        return cap(a, p4.rho4(v))

    def B_bound_u(u):
        return cap(b, p4.sig3(u))

    def B_bound_v(v):
        return cap(b, p4.sig4(v))

    if branch == 1:
        # X(A') * Y(B'), each factor landing on whichever latent it depends on
        fu, fv = one, one
        if sA == "a":
            fu = _scale(fu, p4.X((a, True)))
        elif sA == "u":
            fu = _mul(fu, as_u(lambda u: p4.X(A_bound_u(u))))
        else:
            fv = _mul(fv, as_v(lambda v: p4.X(A_bound_v(v))))
        if sB == "b":
            fu = _scale(fu, p4.Y((b, True)))
        elif sB == "u":
            fu = _mul(fu, as_u(lambda u: p4.Y(B_bound_u(u))))
        else:
            fv = _mul(fv, as_v(lambda v: p4.Y(B_bound_v(v))))
        return [(fu, fv)]

    # branch 2: P(A') + Q(B'), one primitive each
    out: List[Tuple[Step, Step]] = []
    if sA == "a":
        out.append((Step.const(p4.P((a, True))), one))
    elif sA == "u":
        out.append((as_u(lambda u: p4.P(A_bound_u(u))), one))
    else:
        out.append((one, as_v(lambda v: p4.P(A_bound_v(v)))))
    if sB == "b":
        out.append((Step.const(p4.Q((b, True))), one))
    elif sB == "u":
        out.append((as_u(lambda u: p4.Q(B_bound_u(u))), one))
    else:
        out.append((one, as_v(lambda v: p4.Q(B_bound_v(v)))))
    return out


def _mul(s: Step, t: Step) -> Step:
    xs = tuple(sorted(set(s.xs) | set(t.xs)))
    vs = [s(0.0) * t(0.0)] + [s(x) * t(x) for x in xs]
    return Step(xs, tuple(vs))


def _scale(s: Step, k: float) -> Step:
    return Step(s.xs, tuple(v * k for v in s.vs))

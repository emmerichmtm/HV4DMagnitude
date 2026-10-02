"""The grounded 4-D six-staircase model and its claimed four-prefix normal form.

The object under test is

    F(x1,x2,x3,x4) = prod_{i<j} [ x_j < f_ij(x_i) ]

(strict, see ``identities``), and the quantity of interest is the four-prefix
integral ``K(a,b,c,d) = int_{[0,a]x[0,b]x[0,c]x[0,d]} F d mu``.

The report's claim, restated operationally
------------------------------------------
Pair the coordinates ``(x1,x2) | (x3,x4)``.  Integrating the target pair on the
outside leaves, for each target point ``q = (u,v)``, a two-dimensional source
mass ``S(u,v; a,b)`` below the staircase ``f = f12``, restricted to

    x1 in [0,a] cap [0, rho3(u)) cap [0, rho4(v))
    x2 in [0,b] cap [0, sig3(u)) cap [0, sig4(v))

The claim has three parts, and this module tests each separately:

1. **Separability.**  The target plane is cut by a *fixed* set of curves -- which
   of ``rho3(u)``/``rho4(v)`` is the active cap, likewise for ``sig``, and which
   branch of the two-branch staircase formula applies.  On each resulting region
   ``S`` is either ``X(A')*Y(B')`` (branch 1) or ``P(A') + Q(B')`` (branch 2),
   where ``X, Y, P, Q`` are *unary* functions precomputed once, and ``A'``,
   ``B'`` each depend on at most one of ``u, v``.

2. **The prefix arguments only cap.**  ``a`` and ``b`` enter solely through
   anchored rectangles in the target plane, so they change the limits of
   integration and never the region shapes or the unary functions.

3. **Constant region count.**  The number of regions depends on the dimension
   only, never on ``n`` or on the staircase complexity.

Claim 1 is the one that would sink the construction if false, so it is checked
pointwise against a direct sweep rather than inferred.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence, Tuple

from core import INF, Measure1D, Step, brute_integral
from identities import _cum, prefix_bound_gt

__all__ = ["SixStaircases", "Prefix4", "Bound"]


# A prefix of the half-line, carrying whether its endpoint is included.  The
# distinction is invisible under Lebesgue and decisive under magnitude at 0.
Bound = Tuple[float, bool]


def cap(inclusive: float, exclusive: float) -> Bound:
    """Intersect ``[0, inclusive]`` with ``[0, exclusive)``."""
    if inclusive < exclusive:
        return (inclusive, True)
    return (exclusive, False)


def mass(mu: Measure1D, bound: Bound) -> float:
    limit, closed = bound
    return mu.prefix_closed(limit) if closed else mu.prefix_open(limit)


@dataclass
class SixStaircases:
    """The easy residual of a grounded 4-D cell after slab removal."""

    f12: Step
    f13: Step
    f14: Step
    f23: Step
    f24: Step
    f34: Step

    def indicator(self, x1: float, x2: float, x3: float, x4: float) -> float:
        return 1.0 if (x2 < self.f12(x1) and x3 < self.f13(x1)
                       and x4 < self.f14(x1) and x3 < self.f23(x2)
                       and x4 < self.f24(x2) and x4 < self.f34(x3)) else 0.0

    def all_values(self) -> List[float]:
        vals = set()
        for f in (self.f12, self.f13, self.f14, self.f23, self.f24, self.f34):
            vals.update(v for v in f.vs if v != INF)
            vals.update(f.xs)
        return sorted(vals)


class Prefix4:
    """Four-prefix evaluator built from six staircases, with instrumentation.

    ``counters`` records what the structural claims are about: how many regions
    the target plane was cut into, and how many stored records the state uses.
    """

    def __init__(self, st: SixStaircases, measures: Sequence[Measure1D],
                 grids: Sequence[Sequence[float]]):
        self.st = st
        self.mu1, self.mu2, self.mu3, self.mu4 = measures
        self.g1, self.g2, self.g3, self.g4 = [list(g) for g in grids]
        self.counters: Dict[str, int] = {"regions": 0, "records": 0,
                                         "source_evals": 0}
        self._build_unary()

    # -- the unary records: X, Y, P, Q ------------------------------------- #

    def _build_unary(self) -> None:
        """Precompute the four unary functions of the two-branch formula.

        ``X`` and ``Y`` are prefix masses; ``P(A) = int_{[0,A]} Y_open(f(x)) dmu1``
        and ``Q`` is the branch-2 correction.  All four are functions of a single
        argument -- this is exactly the content of the separability claim.

        ``P`` is stored as an *interval record* exactly as the report specifies:
        on each grid interval the integrand is constant, so the cumulative is
        affine there and is recovered from its left-endpoint value plus the
        measure of a subinterval.  Evaluation is then one predecessor search
        rather than a fresh sweep, which is what makes a query cheap.
        """
        f = self.st.f12
        self._P_grid = sorted(set(self.g1) | set(f.xs))
        # cumulative of Y_open(f12(x)) at each grid point, plus the local rate
        self._P_x: List[float] = [0.0]
        self._P_y: List[float] = [self.mu1.atom0 * self.mu2.prefix_open(f(0.0))]
        self._P_rate: List[float] = []
        running = self._P_y[0]
        prev = 0.0
        for pt in [g for g in self._P_grid if g > 0.0]:
            rate = self.mu1.density * self.mu2.prefix_open(f((prev + pt) / 2.0))
            self._P_rate.append(rate)
            running += rate * (pt - prev)
            self._P_x.append(pt)
            self._P_y.append(running)
            prev = pt
        self._P_rate.append(
            self.mu1.density * self.mu2.prefix_open(f(prev + 1.0)))
        self.counters["records"] += len(self._P_grid) + len(self.g2)

    def X(self, bound: Bound) -> float:
        return mass(self.mu1, bound)

    def Y(self, bound: Bound) -> float:
        return mass(self.mu2, bound)

    def P(self, bound: Bound) -> float:
        """``int Y_open(f12(x)) dmu1`` over the prefix, from the stored record.

        One predecessor search plus an affine step inside the interval.  The
        open/closed flag only matters at 0, where the anchor atom lives.
        """
        limit, closed = bound
        if limit < 0.0:
            return 0.0
        if limit == 0.0:
            return self._P_y[0] if closed else 0.0
        if limit == INF:
            return INF if self._P_rate[-1] else self._P_y[-1]
        i = bisect_right(self._P_x, limit) - 1
        if i < 0:
            return 0.0
        return self._P_y[i] + self._P_rate[i] * (limit - self._P_x[i])

    def Q(self, bound_b: Bound) -> float:
        """Branch-2 correction ``X(tau)Y(B) - P(tau)``; unary in ``B``."""
        B, _ = bound_b
        tau = prefix_bound_gt(self.st.f12, B)
        return (self.mu1.prefix_open(tau) * self.Y(bound_b)
                - self.P((tau, False)))

    # -- target-plane geometry --------------------------------------------- #

    def rho3(self, u: float) -> float:
        return prefix_bound_gt(self.st.f13, u)

    def rho4(self, v: float) -> float:
        return prefix_bound_gt(self.st.f14, v)

    def sig3(self, u: float) -> float:
        return prefix_bound_gt(self.st.f23, u)

    def sig4(self, v: float) -> float:
        return prefix_bound_gt(self.st.f24, v)

    def region_of(self, u: float, v: float, a: float, b: float
                  ) -> Tuple[str, str, int]:
        """Label of the fixed region containing the target point ``(u, v)``.

        The first two entries say which cap is active for the source prefixes;
        the third is the branch of the two-branch formula.  ``'a'``/``'b'``
        label the rectangle regions where the prefix argument itself is active.
        """
        r3, r4 = self.rho3(u), self.rho4(v)
        s3, s4 = self.sig3(u), self.sig4(v)
        su = "a" if a < min(r3, r4) else ("u" if r3 <= r4 else "v")
        sv = "b" if b < min(s3, s4) else ("u" if s3 <= s4 else "v")
        A = cap(a, min(r3, r4))
        B = cap(b, min(s3, s4))
        branch = 1 if self._branch1(A, B) else 2
        return (su, sv, branch)

    def _branch1(self, A: Bound, B: Bound) -> bool:
        """Does the staircase cap never bite inside the source prefix box?

        Branch 1 needs every admissible ``y`` to lie strictly below every
        admissible value of ``f12``, i.e. ``sup y < inf f``.

        The prefix ``[0, 0]`` is *not* an empty source set: under the magnitude
        measure it is the anchor atom and carries mass 1, so the comparison
        against ``B`` must still be made.  Only ``[0, 0)`` is empty.  Treating
        the two alike is invisible under Lebesgue and wrong under magnitude; it
        is pinned by a regression test.
        """
        limit, closed = A
        if limit < 0.0 or (limit == 0.0 and not closed):
            return True                   # empty source prefix: vacuously true
        f = self.st.f12
        # infimum of f over the admissible x1-set
        finf = f(limit) if closed else self._left_value(f, limit)
        Bl, Bc = B
        return finf > Bl if Bc else finf >= Bl

    @staticmethod
    def _left_value(f: Step, t: float) -> float:
        """``lim_{x -> t^-} f(x)`` for a right-continuous step function."""
        below = [b for b in f.xs if b < t]
        return f(below[-1]) if below else f(0.0)

    # -- the source mass, two ways ----------------------------------------- #

    def source_separable(self, u: float, v: float, a: float, b: float) -> float:
        """``S(u,v;a,b)`` from the unary records and the region label.

        This is the formula the normal form depends on: within a region, the
        answer is built from ``X, Y, P, Q`` evaluated at arguments that each
        depend on at most one target coordinate.
        """
        self.counters["source_evals"] += 1
        A = cap(a, min(self.rho3(u), self.rho4(v)))
        B = cap(b, min(self.sig3(u), self.sig4(v)))
        if A[0] <= 0.0 and not A[1]:
            return 0.0
        if B[0] <= 0.0 and not B[1]:
            return 0.0
        if self._branch1(A, B):
            return self.X(A) * self.Y(B)
        return self.P(A) + self.Q(B)

    def source_sweep(self, u: float, v: float, a: float, b: float) -> float:
        """``S(u,v;a,b)`` by direct two-dimensional enumeration (ground truth)."""
        A = cap(a, min(self.rho3(u), self.rho4(v)))
        B = cap(b, min(self.sig3(u), self.sig4(v)))
        lim_a, closed_a = A
        lim_b, closed_b = B
        if lim_a < 0 or lim_b < 0:
            return 0.0

        def integrand(x1: float, x2: float) -> float:
            if not closed_a and x1 >= lim_a:
                return 0.0
            if not closed_b and x2 >= lim_b:
                return 0.0
            return 1.0 if x2 < self.st.f12(x1) else 0.0

        gx = list(self.g1) + [lim_a] + list(self.st.f12.xs)
        gy = list(self.g2) + [lim_b] + [w for w in self.st.f12.vs if w != INF]
        return brute_integral(integrand, [gx, gy], [self.mu1, self.mu2],
                              [lim_a, lim_b])

    # -- the four-prefix integral ------------------------------------------ #

    def K_decomposed(self, a: float, b: float, c: float, d: float) -> float:
        """``K(a,b,c,d)`` as an outer target integral of the separable source."""
        def integrand(u: float, v: float) -> float:
            if not (v < self.st.f34(u)):
                return 0.0
            return self.source_separable(u, v, a, b)

        gu = list(self.g3) + [c] + list(self.st.f34.xs)
        gv = list(self.g4) + [d] + [w for w in self.st.f34.vs if w != INF]
        gu += [w for w in self.st.f13.vs if w != INF]
        gu += [w for w in self.st.f23.vs if w != INF]
        gv += [w for w in self.st.f14.vs if w != INF]
        gv += [w for w in self.st.f24.vs if w != INF]
        return brute_integral(integrand, [gu, gv], [self.mu3, self.mu4], [c, d])

    def K_brute(self, a: float, b: float, c: float, d: float) -> float:
        """``K(a,b,c,d)`` by direct four-dimensional enumeration."""
        vals = self.st.all_values()
        grids = [list(self.g1) + vals, list(self.g2) + vals,
                 list(self.g3) + vals, list(self.g4) + vals]
        return brute_integral(self.st.indicator, grids,
                              [self.mu1, self.mu2, self.mu3, self.mu4],
                              [a, b, c, d])

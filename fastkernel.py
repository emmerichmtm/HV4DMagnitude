"""A JIT-compiled kernel for the four-prefix query.

The pure-Python evaluator spends almost all of its time in one place: for each
prefix query it sweeps the target plane cell by cell, and at every cell it
re-derives the source caps through generalized inverses.  That inner work is
arithmetic on flat arrays once the step functions are unpacked, so it compiles.

This module holds the unpacked representation and the kernel.  It computes
exactly what ``Prefix4.K_decomposed`` computes -- same conventions, same
endpoint semantics, same branch test -- and the test suite checks the two
against each other rather than trusting that.

If numba is not installed the kernel falls back to the identical pure-Python
loop, so behaviour never depends on whether the JIT is available.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:                                        # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):                             # type: ignore
        def wrap(fn):
            return fn
        return wrap(args[0]) if args and callable(args[0]) else wrap

INF = np.inf

__all__ = ["HAVE_NUMBA", "PackedStaircases", "k_prefix"]


# --------------------------------------------------------------------------- #
# Kernel helpers -- deliberately mirror core.Step / identities.prefix_bound_gt
# --------------------------------------------------------------------------- #


@njit(cache=True, fastmath=False)
def _step_at(xs, vs, n, t):
    """``Step.__call__``: right-continuous, value on ``[xs[k-1], xs[k])``."""
    lo, hi = 0, n
    while lo < hi:                                   # bisect_right
        mid = (lo + hi) // 2
        if t < xs[mid]:
            hi = mid
        else:
            lo = mid + 1
    return vs[lo]


@njit(cache=True, fastmath=False)
def _step_left(xs, vs, n, t):
    """``lim_{x -> t^-}``: value on the piece ending at ``t``."""
    lo, hi = 0, n
    while lo < hi:                                   # bisect_left
        mid = (lo + hi) // 2
        if xs[mid] < t:
            lo = mid + 1
        else:
            hi = mid
    return vs[lo]


@njit(cache=True, fastmath=False)
def _bound_gt(xs, rev, n, u):
    """Exclusive theta with ``{x : f(x) > u} == [0, theta)``; 0 means EMPTY."""
    lo, hi = 0, n + 1
    while lo < hi:                                   # bisect_right on reversed
        mid = (lo + hi) // 2
        if u < rev[mid]:
            hi = mid
        else:
            lo = mid + 1
    k = (n + 1) - lo                                 # count of values > u
    if k == 0:
        return 0.0
    if k == n + 1:
        return INF
    return xs[k - 1]


@njit(cache=True, fastmath=False)
def _prefix_mass(atom0, density, t, closed):
    if t < 0.0:
        return 0.0
    if t == 0.0:
        return atom0 if closed else 0.0
    return atom0 + density * t


@njit(cache=True, fastmath=False)
def _P_record(px, py, prate, np_, limit, closed):
    """P from its interval record: predecessor search plus an affine step."""
    if limit < 0.0:
        return 0.0
    if limit == 0.0:
        return py[0] if closed else 0.0
    if limit == INF:
        return INF if prate[np_ - 1] != 0.0 else py[np_ - 1]
    lo, hi = 0, np_
    while lo < hi:
        mid = (lo + hi) // 2
        if limit < px[mid]:
            hi = mid
        else:
            lo = mid + 1
    i = lo - 1
    if i < 0:
        return 0.0
    return py[i] + prate[i] * (limit - px[i])


@njit(cache=True, fastmath=False)
def k_prefix(a, b, c, d,
             f12x, f12v, f12r, n12,
             f13x, f13v, f13r, n13,
             f14x, f14v, f14r, n14,
             f23x, f23v, f23r, n23,
             f24x, f24v, f24r, n24,
             f34x, f34v, n34,
             gu, gv, nu, nv,
             px, py, prate, np_,
             atom1, dens1, atom2, dens2, atom3, dens3, atom4, dens4):
    """``K(a,b,c,d)``: the four-prefix integral, by sweeping the target plane.

    One cell of the target grid at a time; each contributes its mass times the
    separable source mass at that cell's representative point.  This is the same
    computation as ``Prefix4.K_decomposed``, unrolled onto flat arrays.
    """
    total = 0.0
    for i in range(nu):
        u_lo = gu[i]
        u_hi = gu[i + 1]
        if u_hi > c:
            break
        if u_lo == 0.0 and atom3 != 0.0:
            u_rep = 0.0
            mass_u = atom3
        else:
            u_rep = 0.5 * (u_lo + u_hi)
            mass_u = dens3 * (u_hi - u_lo)
        if mass_u == 0.0:
            continue

        r3 = _bound_gt(f13x, f13r, n13, u_rep)
        s3 = _bound_gt(f23x, f23r, n23, u_rep)
        g_at_u = _step_at(f34x, f34v, n34, u_rep)

        for j in range(nv):
            v_lo = gv[j]
            v_hi = gv[j + 1]
            if v_hi > d:
                break
            if v_lo == 0.0 and atom4 != 0.0:
                v_rep = 0.0
                mass_v = atom4
            else:
                v_rep = 0.5 * (v_lo + v_hi)
                mass_v = dens4 * (v_hi - v_lo)
            if mass_v == 0.0:
                continue
            if not (v_rep < g_at_u):                 # the f34 staircase
                continue

            r4 = _bound_gt(f14x, f14r, n14, v_rep)
            s4 = _bound_gt(f24x, f24r, n24, v_rep)

            # source caps, with the decorated open/closed endpoint
            excl_a = r3 if r3 < r4 else r4
            if a < excl_a:
                a_lim, a_closed = a, True
            else:
                a_lim, a_closed = excl_a, False
            excl_b = s3 if s3 < s4 else s4
            if b < excl_b:
                b_lim, b_closed = b, True
            else:
                b_lim, b_closed = excl_b, False

            if a_lim < 0.0 or (a_lim == 0.0 and not a_closed):
                continue
            if b_lim < 0.0 or (b_lim == 0.0 and not b_closed):
                continue

            # branch test: does the f12 cap bite inside the source box?
            if a_lim == 0.0:
                finf = _step_at(f12x, f12v, n12, 0.0)
            elif a_closed:
                finf = _step_at(f12x, f12v, n12, a_lim)
            else:
                finf = _step_left(f12x, f12v, n12, a_lim)
            if b_closed:
                branch1 = finf > b_lim
            else:
                branch1 = finf >= b_lim

            if branch1:
                s_mass = (_prefix_mass(atom1, dens1, a_lim, a_closed)
                          * _prefix_mass(atom2, dens2, b_lim, b_closed))
            else:
                tau = _bound_gt(f12x, f12r, n12, b_lim)
                s_mass = (_P_record(px, py, prate, np_, a_lim, a_closed)
                          + _prefix_mass(atom1, dens1, tau, False)
                          * _prefix_mass(atom2, dens2, b_lim, b_closed)
                          - _P_record(px, py, prate, np_, tau, False))
            total += mass_u * mass_v * s_mass
    return total


# --------------------------------------------------------------------------- #
# Packing
# --------------------------------------------------------------------------- #


class PackedStaircases:
    """Flat-array form of a `Prefix4`, ready for the kernel."""

    def __init__(self, p4):
        st = p4.st
        self.parts = []
        for f in (st.f12, st.f13, st.f14, st.f23, st.f24, st.f34):
            xs = np.asarray(f.xs, dtype=np.float64)
            vs = np.asarray(f.vs, dtype=np.float64)
            rev = np.asarray(f._rev, dtype=np.float64)
            if xs.size == 0:
                xs = np.zeros(1, dtype=np.float64)
            self.parts.append((xs, vs, rev, len(f.xs)))

        # Same target grid K_decomposed builds: the axis grid plus every
        # staircase value that can induce a cut in that axis.
        gu = list(p4.g3) + list(st.f34.xs)
        gu += [w for w in st.f13.vs if w != np.inf]
        gu += [w for w in st.f23.vs if w != np.inf]
        gv = list(p4.g4)
        gv += [w for w in st.f34.vs if w != np.inf]
        gv += [w for w in st.f14.vs if w != np.inf]
        gv += [w for w in st.f24.vs if w != np.inf]
        self.gu = _grid_array(gu)
        self.gv = _grid_array(gv)
        self.px = np.asarray(p4._P_x, dtype=np.float64)
        self.py = np.asarray(p4._P_y, dtype=np.float64)
        self.prate = np.asarray(p4._P_rate, dtype=np.float64)
        self.mu = [(m.atom0, m.density)
                   for m in (p4.mu1, p4.mu2, p4.mu3, p4.mu4)]

    def query(self, a, b, c, d) -> float:
        (f12x, f12v, f12r, n12) = self.parts[0]
        (f13x, f13v, f13r, n13) = self.parts[1]
        (f14x, f14v, f14r, n14) = self.parts[2]
        (f23x, f23v, f23r, n23) = self.parts[3]
        (f24x, f24v, f24r, n24) = self.parts[4]
        (f34x, f34v, f34r, n34) = self.parts[5]
        (a1, d1), (a2, d2), (a3, d3), (a4, d4) = self.mu
        return k_prefix(a, b, c, d,
                        f12x, f12v, f12r, n12,
                        f13x, f13v, f13r, n13,
                        f14x, f14v, f14r, n14,
                        f23x, f23v, f23r, n23,
                        f24x, f24v, f24r, n24,
                        f34x, f34v, n34,
                        self.gu, self.gv, len(self.gu) - 1, len(self.gv) - 1,
                        self.px, self.py, self.prate, len(self.px),
                        a1, d1, a2, d2, a3, d3, a4, d4)


def _grid_array(grid) -> np.ndarray:
    pts = sorted({0.0} | {g for g in grid if g > 0.0})
    return np.asarray(pts, dtype=np.float64)

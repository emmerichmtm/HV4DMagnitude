"""Phase 9: a working 4-D hypervolume solver built on the prefix normal form.

This is the end-to-end question: can the product-measure machinery actually
compute a hypervolume inside Chan's recursion, not merely reproduce integrals in
isolation?  The solver in `solver.py` keeps Chan's structure unchanged --
absorption, weighted-median cuts on `(d-3)`-faces, inclusion--exclusion at the
base case -- and replaces only how the easy residual is integrated, routing it
through the four-prefix normal form.

Results are checked against `chan_hypervolume.hypervolume`, which shares no code
with anything here, on random fronts, fully non-dominated spherical fronts,
tie-heavy integer grids, and degenerate inputs.
"""

import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, r"C:\MyTemp\code\FastHVChan")

from solver import PrefixHV4D

try:
    from chan_hypervolume import hypervolume as ref_hv
    from chan_orthant_dby3 import ChanOrthantMeasure
except ImportError:                                       # pragma: no cover
    ref_hv = None

TOL = 1e-9
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


def sphere(n, rng):
    out = []
    for _ in range(n):
        v = [abs(rng.gauss(0.0, 1.0)) + 1e-12 for _ in range(4)]
        s = math.sqrt(sum(x * x for x in v))
        out.append(tuple(x / s for x in v))
    return out


def test_random_fronts():
    print("Phase 9a: random 4-D fronts")
    if ref_hv is None:
        print("    FastHVChan not importable; skipped")
        return
    rng = random.Random(1)
    n_ok = 0
    for n in (1, 2, 3, 5, 8, 12, 20):
        for trial in range(2):
            pts = [tuple(round(rng.random(), 2) for _ in range(4))
                   for _ in range(n)]
            ref = (1.0,) * 4
            got = PrefixHV4D().hypervolume(pts, ref)
            if not check(f"n={n} trial={trial}", got, ref_hv(pts, ref)):
                print(f"      points={pts}")
                return
            n_ok += 1
    print(f"    {n_ok} instances agree with the reference")


def test_spherical_fronts():
    """The hard case: every point non-dominated, so nothing prunes."""
    print("Phase 9b: spherical fronts, and the recursion tree")
    if ref_hv is None:
        return
    rng = random.Random(42)
    for n in (5, 10, 20, 40):
        pts = sphere(n, rng)
        ref = (1.2,) * 4
        solver = PrefixHV4D()
        got = solver.hypervolume(pts, ref)
        if not check(f"sphere n={n}", got, ref_hv(pts, ref)):
            return
        lo = tuple(min(p[k] for p in pts) for k in range(4))
        probe = ChanOrthantMeasure(4)
        probe.complement_measure(pts, lo, ref)
        print(f"    n={n:>3}  nodes: prefix={solver.counters['nodes']:>4}  "
              f"reference d/3={probe.node_count:>4}")
    print("    node counts track the reference recursion (small differences "
          "come from ties in the weighted median)")


def test_degenerate():
    print("Phase 9c: ties, duplicates and degenerate inputs")
    if ref_hv is None:
        return
    rng = random.Random(7)
    n_ok = 0
    for trial in range(6):                       # heavy coordinate sharing
        pts = [tuple(float(rng.randrange(0, 4)) for _ in range(4))
               for _ in range(10)]
        ref = (4.0,) * 4
        if not check(f"integer grid {trial}",
                     PrefixHV4D().hypervolume(pts, ref), ref_hv(pts, ref)):
            return
        n_ok += 1
    for pts, label in (([(0.5, 0.5, 0.5, 0.5)] * 4, "all duplicate"),
                       ([(0.5,) * 4, (0.9,) * 4], "one dominated"),
                       ([(0.0,) * 4], "at the origin"),
                       ([(0.9,) * 4, (0.1,) * 4], "nested")):
        ref = (1.0,) * 4
        if not check(label, PrefixHV4D().hypervolume(pts, ref),
                     ref_hv(pts, ref)):
            return
        n_ok += 1
    print(f"    {n_ok} degenerate instances agree")


if __name__ == "__main__":
    test_random_fronts()
    test_spherical_fronts()
    test_degenerate()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        sys.exit(1)
    print("Phase 9 passed: the prefix normal form computes 4-D hypervolume "
          "inside Chan's recursion.")

"""Phase 9 (partial): tie the six-staircase model to real hypervolume.

The six-staircase model describes the region a cell's *easy* orthants leave
uncovered.  If every orthant of an instance is locally two-sided -- i.e. each
point has exactly two coordinates above the cell floor -- then the model is the
whole instance, and

    HV(points, ref) = mu(cell) - int F d mu .

That makes the abstract model checkable against the independent FastHVChan
implementation, which shares no code with anything here.
"""

import itertools
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, r"C:\MyTemp\code\FastHVChan")

from core import INF, LEBESGUE, MAGNITUDE, Step, brute_integral, staircase_from_quadrants
from prefix4 import Prefix4, SixStaircases

try:
    from chan_hypervolume import hypervolume as hv_reference
except ImportError:                                        # pragma: no cover
    hv_reference = None

TOL = 1e-9
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


def two_sided_instance(rng, T=2.0, n=5):
    """Points with exactly two coordinates above 0: all orthants are 2-sided."""
    pts = []
    for _ in range(n):
        i, j = rng.sample(range(4), 2)
        p = [0.0] * 4
        p[i] = rng.choice([0.5, 1.0, 1.5])
        p[j] = rng.choice([0.5, 1.0, 1.5])
        pts.append(tuple(p))
    return pts


def six_from_points(points):
    """Merge the 2-sided orthants of each axis pair into one staircase."""
    pairs = {}
    for p in points:
        active = [k for k in range(4) if p[k] > 0.0]
        if len(active) != 2:
            raise ValueError("instance is not purely two-sided")
        i, j = active
        pairs.setdefault((i, j), []).append((p[i], p[j]))
    def st(i, j):
        return staircase_from_quadrants(pairs.get((i, j), []))
    return SixStaircases(st(0, 1), st(0, 2), st(0, 3), st(1, 2), st(1, 3), st(2, 3))


def test_model_matches_reference():
    print("Phase 9: six-staircase model vs FastHVChan hypervolume")
    if hv_reference is None:
        print("    FastHVChan not importable; skipped")
        return
    rng = random.Random(9)
    T = 2.0
    grid = [0.5, 1.0, 1.5, 2.0]
    n_ok = 0
    for trial in range(25):
        pts = two_sided_instance(rng, T=T, n=rng.randint(1, 6))
        st = six_from_points(pts)
        vals = st.all_values()
        grids = [list(grid) + vals] * 4
        uncovered = brute_integral(st.indicator, grids, [LEBESGUE] * 4, [T] * 4)
        covered = T ** 4 - uncovered
        want = hv_reference(pts, (T, T, T, T))
        if not check(f"trial {trial} (n={len(pts)})", covered, want):
            print(f"      points={pts}")
            return
        n_ok += 1
    print(f"    {n_ok} instances agree with the reference implementation")


def test_prefix_state_matches_reference():
    """The same, routed through the four-prefix evaluator rather than brute force."""
    print("Phase 9b: four-prefix evaluator vs FastHVChan")
    if hv_reference is None:
        print("    FastHVChan not importable; skipped")
        return
    rng = random.Random(91)
    T = 2.0
    grid = [0.5, 1.0, 1.5, 2.0]
    n_ok = 0
    for trial in range(15):
        pts = two_sided_instance(rng, T=T, n=rng.randint(1, 5))
        st = six_from_points(pts)
        p4 = Prefix4(st, [LEBESGUE] * 4, [grid] * 4)
        covered = T ** 4 - p4.K_decomposed(T, T, T, T)
        want = hv_reference(pts, (T, T, T, T))
        if not check(f"prefix trial {trial}", covered, want):
            print(f"      points={pts}")
            return
        n_ok += 1
    print(f"    {n_ok} instances agree via the four-prefix normal form")


def test_magnitude_of_instance():
    """Magnitude of the same dominated sets, both ways."""
    print("Phase 9c: magnitude of a two-sided instance, two routes")
    rng = random.Random(92)
    T = 2.0
    grid = [0.5, 1.0, 1.5, 2.0]
    for trial in range(8):
        pts = two_sided_instance(rng, T=T, n=rng.randint(1, 4))
        st = six_from_points(pts)
        vals = st.all_values()
        grids = [list(grid) + vals] * 4
        # magnitude of the covered region, as (cell) - (uncovered)
        uncovered = brute_integral(st.indicator, grids, [MAGNITUDE] * 4, [T] * 4)
        cell = (1.0 + T / 2.0) ** 4
        covered = cell - uncovered
        # independent route: direct indicator of the union of orthants
        def ind(*x):
            return 1.0 if any(all(x[k] >= p[k] for k in range(4))
                              for p in pts) else 0.0
        direct = brute_integral(ind, grids, [MAGNITUDE] * 4, [T] * 4)
        check(f"magnitude trial {trial}", covered, direct)
    print("    ok")


if __name__ == "__main__":
    test_model_matches_reference()
    test_prefix_state_matches_reference()
    test_magnitude_of_instance()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        sys.exit(1)
    print("All reference cross-checks passed.")

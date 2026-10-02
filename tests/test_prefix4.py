"""Phases 5-6: the four-prefix normal form against exact 4-D enumeration.

Three things are checked, in increasing order of what they would cost the
construction if false:

* the boundary class ``min(h_up, h_down)`` is closed under new monotone masks;
* the source mass ``S(u,v;a,b)`` really is separable on each fixed region, i.e.
  the unary records ``X, Y, P, Q`` suffice (this is the normal form's core);
* the assembled ``K(a,b,c,d)`` matches brute force, under both measures.
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import INF, LEBESGUE, MAGNITUDE, Step, staircase_from_quadrants
from prefix4 import Prefix4, SixStaircases

TOL = 1e-9
failures = []


def check(name, got, want, tol=TOL):
    if abs(got - want) > tol * max(1.0, abs(want)):
        failures.append((name, got, want))
        print(f"  [FAIL] {name}: got {got!r} want {want!r}")
        return False
    return True


GRID = [0.5, 1.0, 1.5, 2.0]
VALUES = [0.0, 0.5, 1.0, 1.5, 2.0, INF]


def random_staircase(rng, n_break=2, nonincreasing=True):
    """Random monotone step function.

    Built in one go: ``Step`` drops breakpoints that do not change the value,
    so ``xs`` and ``vs`` of two separately generated objects need not match.
    """
    xs = sorted(rng.sample(GRID, n_break))
    vs = sorted((rng.choice(VALUES) for _ in range(n_break + 1)),
                reverse=nonincreasing)
    return Step(tuple(xs), tuple(vs))


def random_six(rng):
    return SixStaircases(*[random_staircase(rng) for _ in range(6)])


# --------------------------------------------------------------------------- #
# Phase 5: boundary class closure
# --------------------------------------------------------------------------- #

def test_boundary_closure():
    print("Phase 5: boundary class min(h_up, h_down) closure")
    rng = random.Random(5)
    fine = [i * 0.25 for i in range(1, 17)]
    n = 0
    for trial in range(200):
        h_up = random_staircase(rng, 2, nonincreasing=False)     # nondecreasing
        h_down = random_staircase(rng, 2)                        # nonincreasing
        h = h_up.minimum(h_down)
        # absorb a sequence of new monotone masks, as repeated compression would
        for _ in range(3):
            if rng.random() < 0.5:
                k = random_staircase(rng, 2, nonincreasing=False)
                h_up = h_up.minimum(k)
            else:
                k = random_staircase(rng, 2)
                h_down = h_down.minimum(k)
            h = h.minimum(k)
            merged = h_up.minimum(h_down)
            for t in fine:
                if not check(f"closure trial {trial} t={t}", merged(t), h(t)):
                    return
                n += 1
            if not (h_up.is_nondecreasing() and h_down.is_nonincreasing()):
                failures.append(("orientation lost", 0, 0))
                print("  [FAIL] orientation lost after absorption")
                return
    print(f"    {n} evaluations ok; orientations preserved")


# --------------------------------------------------------------------------- #
# Phase 6a: separability of the source mass (the core claim)
# --------------------------------------------------------------------------- #

def test_source_separability():
    print("Phase 6a: source mass separable on each fixed region")
    rng = random.Random(6)
    n = 0
    for measures in ([LEBESGUE] * 4, [MAGNITUDE] * 4):
        label = measures[0].name
        for trial in range(25):
            st = random_six(rng)
            p4 = Prefix4(st, measures, [GRID] * 4)
            for a in (0.5, 1.0, 2.0):
                for b in (0.5, 1.5, 2.0):
                    for u in (0.0, 0.25, 0.75, 1.25, 1.75):
                        for v in (0.0, 0.25, 0.75, 1.25, 1.75):
                            got = p4.source_separable(u, v, a, b)
                            want = p4.source_sweep(u, v, a, b)
                            if not check(
                                    f"S {label} trial {trial} a={a} b={b} "
                                    f"u={u} v={v}", got, want):
                                print(f"      region={p4.region_of(u,v,a,b)}")
                                print(f"      f12={st.f12}")
                                print(f"      f13={st.f13} f14={st.f14}")
                                print(f"      f23={st.f23} f24={st.f24}")
                                return
                            n += 1
    print(f"    {n} source evaluations ok")


# --------------------------------------------------------------------------- #
# Phase 6b: region count is constant
# --------------------------------------------------------------------------- #

def test_region_count():
    print("Phase 6b: number of distinct regions is dimension-bounded")
    rng = random.Random(61)
    worst = 0
    for trial in range(40):
        st = random_six(rng)
        p4 = Prefix4(st, [LEBESGUE] * 4, [GRID] * 4)
        labels = set()
        for a in (0.5, 1.0, 2.0):
            for b in (0.5, 1.5, 2.0):
                for u in [i * 0.2 for i in range(11)]:
                    for v in [i * 0.2 for i in range(11)]:
                        labels.add(p4.region_of(u, v, a, b))
        worst = max(worst, len(labels))
    # 3 cap-selectors x 3 x 2 branches = 18 is the structural maximum
    print(f"    max distinct region labels observed: {worst} (structural max 18)")
    if worst > 18:
        failures.append(("region count", worst, 18))
        print("  [FAIL] region count exceeds the structural bound")


# --------------------------------------------------------------------------- #
# Phase 6c: the assembled four-prefix integral
# --------------------------------------------------------------------------- #

def test_four_prefix():
    print("Phase 6c: K(a,b,c,d) against exact 4-D enumeration")
    rng = random.Random(600)
    n = 0
    for measures in ([LEBESGUE] * 4, [MAGNITUDE] * 4):
        label = measures[0].name
        for trial in range(12):
            st = random_six(rng)
            p4 = Prefix4(st, measures, [GRID] * 4)
            for (a, b, c, d) in [(2.0, 2.0, 2.0, 2.0), (1.0, 2.0, 1.5, 0.5),
                                 (0.5, 0.5, 2.0, 2.0), (2.0, 1.0, 0.5, 1.5),
                                 (0.0, 2.0, 2.0, 2.0), (1.5, 1.5, 1.5, 1.5)]:
                got = p4.K_decomposed(a, b, c, d)
                want = p4.K_brute(a, b, c, d)
                if not check(f"K {label} trial {trial} ({a},{b},{c},{d})",
                             got, want):
                    print(f"      f12={st.f12}\n      f13={st.f13}\n"
                          f"      f14={st.f14}\n      f23={st.f23}\n"
                          f"      f24={st.f24}\n      f34={st.f34}")
                    return
                n += 1
    print(f"    {n} four-prefix queries ok")


def test_anchor_atom_regression():
    """Regression: source prefix [0,0] is the anchor atom, not an empty set.

    Minimized from a random magnitude failure.  With ``a = 0`` the admissible
    x1-set is the single point 0, which has mass 1 under magnitude.  Branch 1
    of the two-branch formula must still compare ``B`` against ``f12(0)``;
    short-circuiting on ``limit <= 0`` returns the whole rectangle instead.
    Under Lebesgue the bug is invisible, because {0} is null.
    """
    print("Phase 6d: anchor-atom regression (a = 0 under magnitude)")
    st = SixStaircases(
        f12=Step((0.5, 1.5), (1.0, 0.5, 0.0)),
        f13=Step((1.0,), (2.0, 1.0)),
        f14=Step((0.5,), (1.0, 0.5)),
        f23=Step((0.5, 1.5), (INF, 2.0, 1.5)),
        f24=Step((0.5, 1.0), (INF, 1.5, 1.0)),
        f34=Step((1.5, 2.0), (1.5, 0.5, 0.0)))
    p4 = Prefix4(st, [MAGNITUDE] * 4, [GRID] * 4)
    check("anchor atom S(0,0;a=0,b=2)", p4.source_separable(0.0, 0.0, 0.0, 2.0),
          p4.source_sweep(0.0, 0.0, 0.0, 2.0))
    check("anchor atom K(0,2,2,2)", p4.K_decomposed(0.0, 2.0, 2.0, 2.0),
          p4.K_brute(0.0, 2.0, 2.0, 2.0))
    # and the Lebesgue twin, where the same instance must be identically 0
    p4L = Prefix4(st, [LEBESGUE] * 4, [GRID] * 4)
    check("anchor atom Lebesgue K(0,...)", p4L.K_decomposed(0.0, 2.0, 2.0, 2.0),
          0.0)
    print("    ok")


if __name__ == "__main__":
    test_anchor_atom_regression()
    test_boundary_closure()
    test_source_separability()
    test_region_count()
    test_four_prefix()
    print()
    if failures:
        print(f"{len(failures)} FAILURES")
        sys.exit(1)
    print("All prefix4 tests passed.")

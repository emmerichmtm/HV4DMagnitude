# HV468DMagnitude

Verification of the **general pair-elimination** simplification of Chan's
hypervolume algorithm, for fixed dimensions 4, 6 and 8.

The 4-D case is verified end to end in the parent project. This one tests the
generalization: that eliminating variables from a *weighted monotone-pair term*
stays inside the class, with branching bounded by dimension alone.

**Status: no counterexample found. The elimination calculus is implemented and
exact at `p = 2,3,4,5,6,8` under both Lebesgue and magnitude product measures,
with zero abstentions.** Four blocking defects in the report's written
statements were found and repaired to get there; the compression that *uses*
the calculus is still validated only at `d = 4`.

* [`notes/REPORT_PROBLEMS.md`](notes/REPORT_PROBLEMS.md) — itemized problems
  with the report, by severity
* [`notes/PAIR_ELIMINATION_AUDIT.md`](notes/PAIR_ELIMINATION_AUDIT.md) — the
  audit, finding by finding
* [`PAIR_ELIMINATION_COMPLEXITY.md`](PAIR_ELIMINATION_COMPLEXITY.md) — measured
  growth tables

## Results

Against exact product-grid enumeration, anchor atom as its own cell so
magnitude is exact rather than approximated. Equal instance counts under both
measures.

| `p` | instances | exactness | abstentions | max branching / elimination |
|---|---|---|---|---|
| 2 | 30 | exact | 0 | 2 |
| 3 | 30 | exact | 0 | 3 |
| 4 | 30 | exact | 0 | 6 |
| 5 | 12 | exact | 0 | 4 |
| 6 | 12 | exact | 0 | 14 |
| 8 | 8 | exact | 0 | 3 |

### Which dimension that corresponds to

The report uses these terms in two settings needing different variable counts,
and the distinction is the difference between "validated at `d = 8`" and
"validated at `d = 4`":

| setting | variables | `d=4` | `d=6` | `d=8` |
|---|---|---|---|---|
| static prefix chain | `p = d` | 4 | 6 | 8 |
| compression numerator | `p = 2d` | 8 | 12 | 16 |

So the static chains behind `T_6 = Õ(n²)` and `T_8 = Õ(n^{8/3})` are tested at
their own `p`. Compression is reached only for `d = 4`.

### The claim that matters

The branching constant must not depend on instance size. Holding `p` and the
predicate graph fixed and varying only staircase resolution:

| `p = 4` breakpoints | 6 | 10 | 13 | 22 | 29 | 39 | 53 |
|---|---|---|---|---|---|---|---|
| max branching | 2 | 8 | 10 | 8 | 10 | 9 | 9 |

A 9x increase in breakpoints leaves branching at 9–10. At `p = 6` the maximum
also saturates (18, 19, 16) but the mean still creeps, so saturation is clean at
`p = 4` and only probable at `p = 6`.

Branching stays modest in `p` (2–15 through `p = 9`); **peak state size** is the
real cost driver, reaching 499 terms at `p = 8`.

### Coverage ceiling

The constructor abstains rather than guessing when a feasible set is not an
interval. That is the superlevel set of a *valley* — a pair carrying two
opposite-orientation predicates on the lower side. A *hill* is fine, and is what
the one-turn extension admits.

Coverage is complete through `p = 8`. From `p = 9` it degrades as graphs
densify (4/6, 2/6, 2/6, 1/6 at `p = 9..12`). Closing it means splitting a
two-run feasible set into two signed terms — routine, not implemented.

## Running

```bash
python tests/test_elimination.py
python complexity_probe.py
```

## Layout

```
pair_state.py             weighted monotone-pair terms, decorated bounds
pair_elimination.py       the elimination constructor plus instrumentation
pair_oracles.py           exact product-grid enumeration (atom as its own cell)
complexity_probe.py       branching vs instance size, and vs dimension
tests/                    elimination against the oracles
notes/                    the audit, the problems list, and the task prompt
```

## Two methodology notes worth carrying forward

**The magnitude measure earns its place in the test suite.** Of the six genuine
bugs across the two sprints, four were invisible under Lebesgue and exposed only
by magnitude — all four endpoint or atom errors at `0`. That is an argument for
the product-measure framing on its own terms, independent of the complexity
result.

**A NaN passes a tolerance check.** `abs(got - want) > tol` is `False` when
`got` is NaN, so a full test cycle reported "no mismatches" while returning NaN
throughout. The suite now rejects NaN explicitly.

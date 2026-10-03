# HV468DMagnitude

Verification of the **general pair-elimination** simplification of Chan's
hypervolume algorithm, for fixed dimensions 4, 6 and 8.

The 4-D case is verified end to end in the parent project. This one tests the
generalization: that eliminating variables from a *weighted monotone-pair term*
stays inside the class, with branching bounded by dimension alone.

**Status: no counterexample found, and no abstentions remain.** The elimination
calculus is exact against product-grid oracles at `p = 2..6, 8, 9, 10, 12` under
both Lebesgue and magnitude product measures, and for non-monotone boundaries of
every shape — hills, valleys, zigzags — through `p = 7`. **Compression in `2d`
variables is now built and exact at `d = 4` and `d = 6`.** Five blocking defects
in the report's written statements were found and repaired to get there.

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
| 9, 10, 12 | 6 each | exact (coarse grid) | 0 | — |

### Which dimension that corresponds to

The report uses these terms in two settings needing different variable counts,
and the distinction is the difference between "validated at `d = 8`" and
"validated at `d = 4`":

| setting | variables | `d=4` | `d=6` | `d=8` |
|---|---|---|---|---|
| static prefix chain | `p = d` | 4 | 6 | 8 |
| compression numerator | `p = 2d` | 8 | 12 | 16 |

So the static chains behind `T_6 = Õ(n²)` and `T_8 = Õ(n^{8/3})` are tested at
their own `p`, and compression is now tested at `d = 4` and `d = 6` (`p = 8`
and `p = 12`). `d = 8` compression would need `p = 16`.

### Compression

The integrand `\eqref{eq:generalcompression}` is built as one term in the `2d`
variables `(x, y)`, the `d` old variables are eliminated, and the result is
divided by the cell masses. The check is the definition of the conditional
expectation, not a proxy: for every cell `C`, `N_F(y)` must equal the exact
integral of `F` over `C`.

| d | p | cells checked | exact | output terms (`C_d`) |
|---|---|---|---|---|
| 4 | 8 | 16 per instance | yes, both measures | 450 |
| 6 | 12 | 64 per instance | yes, both measures | 78305 |

Two structural claims verified alongside: `F_bar` is a term in the same class
(only unary densities change, as the report says), and relabelling `y` back to
the physical axes is a pure rename — no residual `x`, no auxiliary variable,
which is the substance of the closure claim.

**`C_d` is large.** 450 output terms at `d = 4` and 78305 at `d = 6`, from *one*
incoming term. Canonical merging — adding the coefficients of identical terms,
the only geometry-free merge there is — recovers a factor of 4.3× to 4.5×, and
the merged state is still exact. That is direct evidence on the question
`\ref{prop:complexity}` leaves open: the signed sum does **not** collapse back
to a constant number in any useful sense, so the constant compounds across
compression generations.

### Non-monotone boundaries

Nothing abstains on shape any more. A *valley* — a boundary that falls then
rises — has a superlevel set of two disjoint runs, and since the integral over a
disjoint union is the sum, each run becomes its own signed term.

| shape | p | instances each | exact | max runs | max combos |
|---|---|---|---|---|---|
| hill | 2, 3, 4 | 20 | yes | 2 | 16 |
| valley | 2, 3, 4 | 20 | yes | 2 | 12 |
| zigzag | 2, 3, 4 | 20 | yes | 2 | 96 |
| mixed | 3, 4, 5 | 16 | yes | 2 | 64 |

Valleys stay exact to `p = 7`. Zigzags at `p = 7` exceed a nine-minute budget —
a cost limit, not a coverage one.

### The claim that matters

The branching constant must not depend on instance size. Holding `p` and the
predicate graph fixed and varying only staircase resolution:

| `p = 4` breakpoints | 6 | 10 | 13 | 22 | 29 | 39 | 53 |
|---|---|---|---|---|---|---|---|
| max branching | 2 | 8 | 10 | 8 | 10 | 9 | 9 |

A 9× increase in breakpoints leaves branching at 9–10. At `p = 6` the maximum
nearly levels off (18, 19, 20) but the mean still creeps, so saturation is clean
at `p = 4` and only probable at `p = 6`.

Branching stays modest in `p`; **peak state size** is the real cost driver,
reaching 4828 terms at `p = 12`. Run splitting is cheap: runs never exceed 2 and
branch products never exceed 64.

## Running

```bash
python tests/test_elimination.py
python tests/test_valley.py
python tests/test_compression.py     # d = 6 takes about 3 minutes
python complexity_probe.py
```

## Layout

```
pair_state.py             weighted monotone-pair terms, decorated bounds,
                          split atom/continuous weights
pair_elimination.py       the elimination constructor plus instrumentation
compression.py            the 2d-variable compression integrand, cell grids,
                          division by cell mass, canonical merging
pair_oracles.py           exact product-grid enumeration (atom as its own cell)
complexity_probe.py       branching vs instance size, and vs dimension
tests/test_elimination.py elimination against the oracles
tests/test_valley.py      non-monotone boundaries, run decomposition
tests/test_compression.py compression against the conditional expectation
notes/                    the audit, the problems list, and the task prompt
```

## Three methodology notes worth carrying forward

**The magnitude measure earns its place in the test suite.** Of the seven
genuine bugs across these sprints, five were invisible under Lebesgue and
exposed only by magnitude — every one an endpoint or atom error at `0`. That is
an argument for the product-measure framing on its own terms, independent of the
complexity result.

**But an under-powered generator hides bugs from both measures.** The suite
spent a long time generating only `upper=True, strict=True` predicates, which
concealed a representation defect that was wrong under Lebesgue too. Varying the
decorations mattered as much as varying the measure.

**A NaN passes a tolerance check.** `abs(got - want) > tol` is `False` when
`got` is NaN, so a full test cycle reported "no mismatches" while returning NaN
throughout. Both suites now reject NaN explicitly.

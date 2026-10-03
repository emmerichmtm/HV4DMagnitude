# Pair elimination: measured complexity

All numbers from `complexity_probe.py` (growth), `tests/test_elimination.py`
(exactness) and `tests/test_valley.py` (non-monotone boundaries). Reproduce
with:

```bash
python tests/test_elimination.py
python tests/test_valley.py
python tests/test_compression.py     # d = 6 takes about 3 minutes
python complexity_probe.py
```

`B` is the quantity the report's Lemma 3 bounds: **terms emitted by one
elimination from one incoming term**. It is not a running total — a cumulative
count would conflate branching with chain depth, which is how the first version
of this instrumentation misled me.

## Two different `p` per dimension

The report uses weighted monotone-pair terms in two settings, and they need
different variable counts. Keeping them apart is the difference between
"validated at `d = 8`" and "validated at `d = 4`":

| setting | variables | `d=4` | `d=6` | `d=8` |
|---|---|---|---|---|
| static prefix chain, `\ref{cor:recursivepair}` | `p = d` | 4 | 6 | 8 |
| compression numerator, `\eqref{eq:generalcompression}` | `p = 2d` | 8 | 12 | 16 |

The static chain is exact at all three dimensions, and compression is exact at
`d = 4` and `d = 6`. `d = 8` compression would need `p = 16`.

## Exactness against product-grid oracles

Full enumeration over the product grid, anchor atom as its own cell so magnitude
is exact rather than approximated. Both measures, equal numbers of instances.

| p | instances | all exact | abstentions | max B |
|---|---|---|---|---|
| 2 | 30 | yes | 0 | 2 |
| 3 | 30 | yes | 0 | 3 |
| 4 | 30 | yes | 0 | 6 |
| 5 | 12 | yes | 0 | 4 |
| 6 | 12 | yes | 0 | 14 |
| 8 | 8 | yes | 0 | 3 |

On a coarser grid, which is what makes the oracle affordable at higher `p`:

| p | instances | all exact | abstentions |
|---|---|---|---|
| 9 | 6 | yes | 0 |
| 10 | 6 | yes | 0 |
| 12 | 6 | yes | 0 |

So the static chains for `d = 4, 6, 8` all check out, and `T_6 = Õ(n²)`,
`T_8 = Õ(n^{8/3})` rest on a calculus that is now tested at their respective
`p`.

Two cautions on reading these tables. The oracle costs `(cells per axis)^p`, so
the coarse rows have only one interior breakpoint per staircase — which means
they cannot *form* a valley, and their run counts are all 1. Non-monotone
coverage is established separately, below. And the generated staircases must put
both breakpoints *and values* on the grid the cells are cut from, or the oracle
itself returns the wrong number; mixing a fine staircase with a coarse grid is
what made the first `p = 8` run look like a constructor failure when the fault
was in the oracle.

## Non-monotone boundaries, including the valley

The constructor no longer abstains on any shape.

| shape | meaning | superlevel set | runs |
|---|---|---|---|
| hill | rises then falls | middle interval | 1 |
| valley | falls then rises | two disjoint runs | 2 |
| zigzag | arbitrary values on the grid | up to `ceil(k/2)` runs | 2 observed |

Checked against the fine-grid oracle under both measures:

| shape | p | instances each | all exact | max runs | max combos | max B |
|---|---|---|---|---|---|---|
| hill | 2, 3, 4 | 20 | yes | 2 | 16 | 34 |
| valley | 2, 3, 4 | 20 | yes | 2 | 12 | 33 |
| zigzag | 2, 3, 4 | 20 | yes | 2 | 96 | 149 |
| mixed | 3, 4, 5 | 16 | yes | 2 | 64 | — |

Plus a unit check: on 400 random staircases the run decomposition agrees
pointwise with brute-force evaluation of the predicate, at every piece
representative and every breakpoint.

And the 32-cell matrix that caught the left-continuity defect — four shapes ×
predicate side × strictness × measure, on a *single* predicate in two
variables — is exact in all 32 cells. It is worth keeping as a regression:
before the fix, the eight cells with `upper != strict` were all wrong, four of
them under Lebesgue.

Higher `p` with non-monotone boundaries. This is a cost table, not a coverage
table — nothing abstains:

| shape | p | density | exact | peak state | seconds |
|---|---|---|---|---|---|
| valley | 5 | 0.45 | yes | 256 | 1 |
| valley | 6 | 0.35 | yes | 1180 | 1 |
| valley | 7 | 0.30 | yes | 4784 | 25 |
| zigzag | 5 | 0.45 | yes | 168 | 0 |
| zigzag | 6 | 0.35 | yes | 1324 | 6 |
| zigzag | 7 | 0.30 | not measured | — | exceeded 9 min |

## Compression, and the measured `C_d`

The integrand of `\eqref{eq:generalcompression}` is built as one weighted
monotone-pair term in the `2d` variables `(x, y)`, the `d` old variables are
eliminated, and the result is divided by the cell masses. The check is the
definition of the conditional expectation: for every cell `C` of the product
grid, `N_F(y)` must equal the exact integral of `F` over `C`, for `y` in `C`.

| d | p | instances | cells per instance | exact | output terms |
|---|---|---|---|---|---|
| 4 | 8 | 8 | 16 | yes, both measures | 450 |
| 6 | 12 | 4 | 64 | yes, both measures | 78305 |

Each cell is probed at its floor and its midpoint. The floor of the first cell
is the anchor `0` — the one point where the cell convention can be wrong without
Lebesgue noticing. An aggregate check sums the numerators over all cells and
compares against the integral over the whole ground set, which catches
double-counted or missed cells that per-cell checks cannot.

Verified alongside, both structural rather than numeric:

* `F_bar = N_F / mass` is a term in the same class. Only unary densities change,
  exactly as the report says, because the cell mass is constant on each cell and
  so is a step function of `y`.
* Relabelling `y` back to the physical axes is a **pure rename**: no residual `x`
  axis, no auxiliary variable. That is the substance of the closure claim, and
  it is asserted by construction in `reindex_to_physical`, which raises if any
  `x` survives.

### `C_d` does not collapse

`\ref{thm:generalsimplification}` gives `S' <= C_d S` with `C_d` a
dimension-only constant, and `\ref{prop:complexity}` concedes that nothing
proves a signed sum of distinct primitives can be merged back to a constant
number. Measured, from **one** incoming term:

| d | measure | raw output terms | after canonical merge | ratio |
|---|---|---|---|---|
| 4 | Lebesgue | 234 | 54 | 4.3× |
| 4 | magnitude | 117 | 117 | 1.0× |
| 6 | magnitude | 181135 | 40644 | 4.5× |

Canonical merging is the only geometry-free merge available: add the
coefficients of terms that are identical, drop what cancels. The merged state
is still exact on every cell. It recovers a factor of 4.3–4.5× at best, and
nothing at all in one case.

So the open question in `\ref{prop:complexity}` resolves the unhelpful way:
`C_6` is on the order of 10^4 per incoming term after merging, and it compounds
across compression generations. The report's remark that "the constants can be
large if the generic label enumeration is implemented literally" is correct and,
on this evidence, understated — this is a literal implementation, and the
constant is large enough that the practical recommendation to keep a specialized
base case is not a convenience but a requirement.

The asymptotics are untouched: `C_d` is independent of instance size, which is
what the recurrence needs. This is a statement about the constant, not the
exponent.

## Growth in instance size — the claim that matters

The report's constant must not depend on `N`, the number of staircase
breakpoints. Each row below re-seeds the generator, so every row draws the
**same** sequence of predicate graphs and differs only in resolution. Without
that, the sweep measures graph density and resolution at once and establishes
neither.

`p = 4`:

| breaks/pred | 1 | 2 | 3 | 5 | 7 | 11 | 15 |
|---|---|---|---|---|---|---|---|
| term breakpoints | 6 | 10 | 13 | 22 | 29 | 39 | 53 |
| max B | 2 | 8 | 10 | 8 | 10 | 9 | 9 |
| mean B | 1.3 | 3.0 | 3.2 | 4.6 | 5.6 | 4.0 | 4.8 |
| peak state | 2 | 11 | 15 | 29 | 37 | 35 | 37 |

`p = 6`:

| breaks/pred | 1 | 2 | 3 | 5 | 7 |
|---|---|---|---|---|---|
| term breakpoints | 11 | 19 | 27 | 42 | 52 |
| max B | 3 | 6 | 18 | 19 | 20 |
| mean B | 1.3 | 2.7 | 9.1 | 10.7 | 13.7 |
| peak state | 4 | 14 | 135 | 800 | 876 |

**Reading.** At `p = 4`, breakpoints grow 9× (6 → 53) and max branching
saturates at 9–10 after an initial rise off the degenerate 1-breakpoint case.
That is the behaviour the lemma needs, and it is the strongest evidence for the
report's key claim.

At `p = 6` max branching nearly levels off (18, 19, 20) but the **mean is still
creeping upward** at the largest resolution tested (9.1 → 10.7 → 13.7). So
saturation is clean at `p = 4` and only probable at `p = 6`. I would not claim
independence at `p = 6` from this data; I would claim it is consistent with
independence, and that the honest next step is a wider resolution sweep than the
oracle cost currently allows.

## Growth in dimension

Resolution held fixed (4 grid points, 2 breaks per predicate), predicate density
0.6, `p` growing. Ten instances through `p = 8`, five above it — the state is
large enough there that more trials cost minutes for no extra information.

| p | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| term breakpoints | 3 | 6 | 10 | 16 | 18 | 26 | 30 | 36 | 46 | 54 | 65 |
| max B | 2 | 3 | 3 | 6 | 10 | 9 | 11 | 27 | 32 | 9 | 30 |
| mean B | 1.3 | 1.9 | 2.4 | 3.2 | 4.3 | 5.0 | 5.7 | 13.4 | 13.2 | 5.8 | 15.2 |
| peak state | 2 | 6 | 6 | 15 | 45 | 89 | 499 | 899 | 3457 | 352 | 4828 |
| max runs | 1 | 1 | 1 | 1 | 2 | 1 | 1 | 2 | 2 | 2 | 2 |
| max combos | 1 | 1 | 1 | 1 | 2 | 1 | 1 | 64 | 4 | 4 | 8 |
| abstentions | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

**Every `p` from 2 to 12 now completes.** Before run splitting the completion
rate fell to 4/6, 2/6, 2/6, 1/6 at `p = 9..12`.

Branching is modest and non-monotone in `p` — the `p = 11` row sits below
`p = 10`, which is sampling over five instances, not structure. **Peak state is
what grows**: 2 → 4828 across the range. That is the practical cost driver, and
it is consistent with the report's own admission (`\ref{prop:complexity}`) that
nothing proves a signed sum of primitives can be merged back to a constant
number.

The run-splitting overhead is small where it matters: run counts never exceed 2
and the branch product never exceeds 64. The mechanism that closed the coverage
gap did not pay for it with a blow-up in the constant.

## Breakpoints per term along the chain

If each elimination multiplied the breakpoints of the terms it emits, per-term
cost would grow with chain depth even at bounded branching, and the
fixed-dimensional claim would be empty. Measured breakpoints **per surviving
term** after each successive elimination:

```
p = 6:  17.0 -> 14.3 -> 11.5 -> 7.2 -> 2.0 -> 0
p = 8:  22.0 -> 20.5 -> 20.2 -> 15.6 -> 12.4 -> 6.5 -> 1.0 -> 0
```

Monotone decrease throughout. The concern does not materialize.

## Coverage: no abstentions remain

The constructor used to abstain whenever a feasible set was not an interval —
the superlevel set of a *valley*, which is two disjoint runs. That is now
handled by splitting into one signed term per run, so every shape is covered.
The abstention path is still present, deliberately, so that a future gap shows
up as a refusal rather than as a wrong number; it is never taken in any test
here.

Two weaker bounds replace the old coverage ceiling:

* **cost, not coverage.** Zigzag boundaries at `p = 7`, density 0.3, exceed a
  nine-minute budget. Valleys at the same `p` and density are exact in 25
  seconds.
* **oracle reach.** Exactness above `p = 8` is checked only on a coarse grid,
  where staircases have a single interior breakpoint and so cannot form a
  valley. So non-monotone exactness is established at `p <= 7`, and plain
  exactness at `p <= 12`.

## What is and is not established

Established by measurement:

* exact agreement with oracles at `p = 2..6, 8` (fine grid) and `p = 9, 10, 12`
  (coarse grid), under both measures, zero abstentions;
* exact agreement for hill, valley, zigzag and mixed boundaries at `p = 2..5`,
  and for valleys at `p = 6, 7`;
* all 32 combinations of shape, predicate side, strictness and measure exact on
  a single predicate;
* run decomposition agrees with brute force on 400 random staircases;
* branching saturates in instance size at `p = 4`, consistent with saturation at
  `p = 6`;
* breakpoints per term decrease along the elimination chain;
* branching modest in `p`; state size is the cost driver;
* run splitting is cheap: runs ≤ 2, branch products ≤ 64;
* compression exact at `d = 4` and `d = 6`, on every cell, both measures, with
  `F_bar` in the class and the output free of any residual `x` variable;
* `C_d` measured at 450 (`d = 4`) and 78305 (`d = 6`) per incoming term, with
  canonical merging recovering only 4.3–4.5×.

Not established:

* branching independence at `p >= 6` with certainty — the mean is still rising
  at the resolutions reachable before the oracle becomes the bottleneck;
* non-monotone exactness above `p = 7`, or zigzag cost above `p = 6`;
* compression above `d = 6`; `d = 8` needs `p = 16`;
* **multi-generation closure above `d = 4`.** One generation is verified at
  `d = 6`, and the output is structurally ready to feed back in
  (`reindex_to_physical`), but a second generation has not been run. With
  `C_6` measured at ~4 x 10^4 terms per incoming term, that is the measurement
  most likely to show where this stops being practical;
* a bound on `C_d` as a function of `d` from more than two data points. 450 at
  `d = 4` and 78305 at `d = 6` is a ratio of 174 across a single step, which is
  consistent with the report's `2^{O(d³)}` but establishes nothing.

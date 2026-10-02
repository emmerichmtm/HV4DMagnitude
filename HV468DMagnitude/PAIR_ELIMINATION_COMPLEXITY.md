# Pair elimination: measured complexity

All numbers from `complexity_probe.py` (growth), `tests/test_elimination.py`
(exactness) and `tests/test_valley.py` (non-monotone boundaries). Reproduce
with:

```bash
python tests/test_elimination.py
python tests/test_valley.py
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

The static chain is exact at all three dimensions. The *variable count* for
`d = 6` compression (`p = 12`) is now reached too, but the compression
integrand itself has not been built.

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
* run splitting is cheap: runs ≤ 2, branch products ≤ 64.

Not established:

* branching independence at `p >= 6` with certainty — the mean is still rising
  at the resolutions reachable before the oracle becomes the bottleneck;
* non-monotone exactness above `p = 7`, or zigzag cost above `p = 6`;
* compression in `2d` variables above `d = 4` (`\eqref{eq:generalcompression}`),
  multi-generation closure above `d = 4`, or the `S' <= C_d S` bound of
  `\ref{thm:generalsimplification}`. `p = 12` reaches the *variable count* for
  `d = 6`, but the compression integrand itself is not yet built;
* `C_d` itself. With the report's own `B_p = 2^{O(p²)}` and `p <= 2d`, the chain
  gives `C_d = 2^{O(d³)}`, which the report should state because it compounds
  across compression generations.

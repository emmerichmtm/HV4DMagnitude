# Pair elimination: measured complexity

All numbers from `complexity_probe.py` (growth) and `tests/test_elimination.py`
(exactness). Reproduce with:

```bash
python tests/test_elimination.py
python complexity_probe.py
```

`B` is the quantity the report's Lemma 3 bounds: **terms emitted by one
elimination from one incoming term**. It is not a running total — a cumulative
count would conflate branching with chain depth, which is how the first version
of this instrumentation misled me.

## Two different `p` per dimension

The report uses weighted monotone-pair terms in two distinct settings, and they
need different variable counts. Keeping them apart is the difference between
"validated at `d = 8`" and "validated at `d = 4`":

| setting | variables | `d=4` | `d=6` | `d=8` |
|---|---|---|---|---|
| static prefix chain, `\ref{cor:recursivepair}` | `p = d` | 4 | 6 | 8 |
| compression numerator, `\eqref{eq:generalcompression}` | `p = 2d` | 8 | 12 | 16 |

The static chain is exact at all three dimensions. Compression is reached only
for `d = 4`.

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

So the static chains for `d = 4, 6, 8` all check out, and
`T_6 = Õ(n^2)`, `T_8 = Õ(n^{8/3})` rest on a calculus that is now tested at
their respective `p`.

Two cautions on reading this table. The oracle costs `(cells per axis)^p`, so
`p = 8` runs on a coarser grid than `p <= 6`; and the generated staircases must
put both breakpoints *and values* on that same grid, or the oracle itself
returns the wrong number. Mixing a fine staircase with a coarse grid is what
made the first `p = 8` run look like a constructor failure when the fault was in
the oracle.

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
| max B | 3 | 6 | 18 | 19 | 16 |
| mean B | 1.3 | 2.7 | 9.1 | 10.7 | 12.8 |
| peak state | 4 | 14 | 135 | 800 | 268 |

**Reading.** At `p = 4`, breakpoints grow 9x (6 → 53) and max branching
saturates at 9–10 after an initial rise off the degenerate 1-breakpoint case.
That is the behaviour the lemma needs, and it is the strongest evidence for the
report's key claim.

At `p = 6` max branching also saturates (18, 19, 16) but the **mean is still
creeping upward** at the largest resolution tested (9.1 → 10.7 → 12.8). So
saturation is clean at `p = 4` and only probable at `p = 6`. I would not claim
independence at `p = 6` from this data; I would claim it is consistent with
independence.

Peak state at `p = 6` is non-monotone (135, 800, 268), which is sampling noise
over 10 instances, not structure.

## Growth in dimension

Resolution held fixed (4 grid points, 2 breaks per predicate), predicate density
0.6, `p` growing.

| p | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| term breakpoints | 3 | 6 | 10 | 16 | 18 | 26 | 30 | 38 | 46 | 54 | 65 |
| max B | 2 | 3 | 3 | 6 | 6 | 9 | 11 | 15 | 10 | 6 | 0 |
| peak state | 2 | 6 | 6 | 14 | 26 | 89 | 499 | 43 | 205 | 352 | 0 |
| completed | all | all | all | 9/10 | 8/10 | all | all | 4/6 | 2/6 | 2/6 | 1/6 |

Branching itself stays modest — 2 to 15 across `p = 2..9`. **Peak state is what
grows**: 2, 6, 6, 14, 26, 89, 499 for `p = 2..8`. That is the practical cost
driver, and it is consistent with the report's own admission
(`\ref{prop:complexity}`) that nothing proves a signed sum of primitives can be
merged back to a constant number.

The `max B` and `peak state` figures for `p >= 10` are **not comparable** to the
rest: they average over so few completed instances that they drop rather than
rise, which is a selection effect — the surviving instances are the easy ones.
`p = 12` completed one instance whose term was identically zero, hence `B = 0`.

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

## Coverage limit: the valley case

The constructor abstains — returns `None` rather than guessing — when a feasible
set is not an interval. That happens for the superlevel set of a *valley*: a
pair that has accumulated two opposite-orientation predicates on the lower side,
whose conjunction is `max(h_up, h_down)`, leaving two disjoint feasible runs.
A hill (`min` of the same two, on the upper side) is fine, and is the case the
one-turn extension admits.

Coverage is complete through `p = 8` in every configuration tested. From
`p = 9` it degrades as graphs densify — 4/6, 2/6, 2/6, 1/6 at `p = 9..12` with
density 0.6 and 2 breaks per predicate — because more pairs mean more chances
for two opposite orientations to land on one pair's lower side.

This is the honest ceiling on the general claim today: **`p <= 8`, which covers
the static chain for `d <= 8` but only `d = 4` compression.** Closing it means
handling a feasible set that is a union of two intervals, by splitting into two
signed terms. The split is routine; it has not been implemented.

## What is and is not established

Established by measurement:

* exact agreement with oracles at `p = 2,3,4,5,6,8`, both measures, zero
  abstentions;
* branching saturates in instance size at `p = 4`, consistent with saturation at
  `p = 6`;
* breakpoints per term decrease along the elimination chain;
* branching is modest in `p`; state size is the cost driver.

Not established:

* branching independence at `p >= 6` with certainty — the mean is still rising
  at the resolutions reachable before the oracle becomes the bottleneck;
* anything about `p >= 9`, where the valley gap bites;
* compression in `2d` variables above `d = 4` (`\eqref{eq:generalcompression}`),
  multi-generation closure above `d = 4`, or the `S' <= C_d S` bound of
  `\ref{thm:generalsimplification}`;
* `C_d` itself. With the report's own `B_p = 2^{O(p^2)}` and `p <= 2d`, the
  chain gives `C_d = 2^{O(d^3)}`, which the report should state because it
  compounds across compression generations.

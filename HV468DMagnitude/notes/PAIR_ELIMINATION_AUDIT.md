# Audit: the general pair-elimination simplification

Report audited: *A Product-Measure and Magnitude Simplification of Chan's
Hypervolume Algorithm* — generalized version with Section "General
fixed-dimensional pair elimination"
(`paper/magnitude_chan_pair_elimination_report.tex`).

**Status: no counterexample found. The elimination calculus is implemented and
exact at `p = 2..6, 8, 9, 10, 12` under both measures, non-monotone boundaries
included; compression in `2d` variables is exact at `d = 4` and `d = 6`. Six
blocking defects in the written statements were repaired along the way.**

The full itemized list of problems with the report, with severities, is in
[`REPORT_PROBLEMS.md`](REPORT_PROBLEMS.md). The measured growth tables are in
[`../PAIR_ELIMINATION_COMPLEXITY.md`](../PAIR_ELIMINATION_COMPLEXITY.md). This
file records what was tested and the findings in the order they were found.

---

## What was tested

Elimination checked against exact product-grid enumeration, anchor atom as its
own cell so magnitude is exact rather than approximated. Equal instance counts
under Lebesgue and under magnitude.

| `p` | instances | exact | abstentions | max branching |
|---|---|---|---|---|
| 2 | 30 | yes | 0 | 2 |
| 3 | 30 | yes | 0 | 3 |
| 4 | 30 | yes | 0 | 6 |
| 5 | 12 | yes | 0 | 4 |
| 6 | 12 | yes | 0 | 14 |
| 8 | 8 | yes | 0 | 3 |

This covers the static prefix chain for `d = 4, 6, 8`. Compression, which needs
`p = 2d`, is covered separately at `d = 4` and `d = 6` -- see Finding 9.

An earlier version of this table read 30/30, 18/30, 11/30, 1/6, **0/6** with 42
abstentions. The difference is Findings 1–3 below.

---

## Finding 1 — the term class is not closed, and the abstentions proved it

Definition 1 requires every `f_e` to be **monotone** and justifies "at most one
lower and one upper predicate per pair" by saying same-orientation predicates
merge via min/max. That merge exists only when the orientations agree.
Elimination routinely produces, for one pair, an upper predicate from a
nondecreasing map *and* one from a nonincreasing map, and their conjunction

```
xi_j  <|  min( h_up(xi_i), h_down(xi_i) )
```

has a turning point. The class as defined does not contain it.

The abstention rate was the proof: 0% at `p = 2`, 40% at `p = 3`, 63% at
`p = 4`, 83% at `p = 5`, 100% at `p = 6` — growing with dimension precisely
because more eliminations mean more chances for two opposite orientations to
land on the same pair.

**Repaired, two ways at once.** A pair now carries a *list* of predicates
(`pair_state.merge_same_orientation` merges only within an orientation), and
one-turn boundaries are admitted, mirroring the 4-D "latent one-turn normal
form". Abstentions went to zero through `p = 8`.

The invariant was load-bearing, which is the part worth emphasizing in the
paper: the `O(p)` candidate count in Lemma 2's proof *is* "one lower and one
upper per incident pair". With a list it becomes `O_p(1)` with a worse constant
— still fine for fixed `p`, but not what the proof says.

## Finding 2 — inversion is interval-valued, not one-sided

Lemma 2's proof says each incident constraint becomes "a lower **or** upper
interval bound". That disjunction is not exhaustive. Even for monotone `f` the
degenerate pieces — all `x` feasible, none feasible — belong to the opposite
family, so forcing one side per predicate misrepresents the rest; and once
one-turn maps are admitted, a hill's superlevel set is a genuine middle interval.

**Repaired.** `_invert_to_bounds` returns a decorated interval: a lower *and* an
upper candidate plus a `0/1` feasibility weight for the empty pieces. The side
condition the construction actually needs is that the feasible set be an
**interval**, tested by contiguity of the feasible run — not by orientation.
This single change took `p = 6` from 0% coverage to 100%.

The *valley* -- `max(h_up, h_down)` on the lower side, whose feasible set is two
disjoint runs -- was the remaining gap, and is now closed. `_feasible_runs`
returns every maximal run, and since the integral over a disjoint union is the
sum, a set with `r` runs becomes `r` signed terms. The run count varies per
piece of the surviving variable, so the split is by run *index* with short
pieces padded by the degenerate `[0, 0)`; and because intersecting two unions of
runs gives the union over choices of one run from each, the branches multiply
across the predicates incident to the eliminated variable. Measured run counts
stay at 2 and branch products at 96 or below. This removed the `p = 9..12`
abstentions entirely.

## Finding 3 — the decorated lower bound inverts in the cumulative difference

The single genuine bug of the first pass. The mass of a feasible interval is

```
Q(U) - Q(L^-)
```

so on the **lower** side the endpoint flag must be *inverted*: if `L` is
attained, what gets subtracted is the mass strictly below it; if `L` is
excluded, the mass up to and including it. Using the flag directly deletes the
anchor atom whenever `L = 0`.

Invisible under Lebesgue — every test passed. Under magnitude it is a mass-1
error that broke every instance from `p = 2` upward.

## Finding 4 — the label partition double-counts ties

"Partition by which lower endpoint is maximal and which upper is minimal" is not
a partition. When two candidates coincide on a region of positive measure both
labels hold and the region is integrated twice — and ties are the common case,
not an edge case, because distinct predicates routinely invert to the same
domain limit.

**Repaired** by first-winner tie-breaking: the candidate at index `k` compares
strictly against earlier candidates and non-strictly against later ones. Before
this, integrals came out exactly 2x and 4x too large.

## Finding 5 — abstention has to be total, or it becomes an overcount

`_impose` previously returned "imposed" when the comparison could not be
expressed, which **deletes a constraint** and inflates the integral. Dropping a
constraint is not a conservative approximation; it is a wrong answer that looks
like a passing test until the class is wide enough to reach it. Latent while
only monotone predicates existed, dominant the moment one-turn predicates were
admitted — every resulting failure had `got > want`.

`_impose` is now tri-state (imposed / identically false / abstain) and the
abstention propagates out of `eliminate_variable`.

## Finding 6 — NaN is not annihilated by multiplying with zero

An inverted bound carries `+inf` where its predicate constrains nothing, and the
cumulative of an unbounded endpoint integrates `density * (inf - g) *
weight(inf)`: `inf`, or `0 * inf = NaN`. The label indicators already exclude
those pieces, but zero times NaN is NaN, so the whole signed sum is poisoned.

**Repaired** by clamping endpoints to the ground set's ceiling at the measure
step — exact, because `x` is confined to `[lo, hi]`, so "unconstrained above"
*is* `hi`. The clamp must not be applied to the candidate before the label
comparisons, where `+inf` correctly means "weaker than every other bound"; the
stronger clamp was tried and broke `p = 4` and `p = 5`.

This one hid behind a test defect: `abs(got - want) > tol` is `False` when `got`
is NaN, so **every NaN silently passed** and one full cycle reported "no
mismatches" while returning NaN throughout. The suite now rejects NaN
explicitly. Any audit of this report should assert non-NaN separately.

Removing the NaNs also cut branching sharply — `p = 4`: 21 → 6, `p = 6`: 50 → 14
— because NaN weights were never recognized as dead terms.

## Finding 7 — `upper != strict` inverts left-continuously, and the anchor needs its own weight

Found by the valley tests, but not a valley bug at all -- a pre-existing one the
earlier suite could not reach, because its generators only ever produced
`upper = strict = True` predicates. With all four decorated comparisons
generated, exactly the two *mixed* combinations failed.

Inverting `v <| f(x)` gives a feasible `x`-set that moves as `v` crosses the
values `f` attains, and the direction of continuity follows the decoration:
`v < f(x)` and `v >= f(x)` are constant on `[a, b)`, matching a right-continuous
`Step`, while `v <= f(x)` and `v > f(x)` are constant on `(a, b]` and cannot be
written that way at all. Sampling the first piece at `0` then spread the `v = 0`
value across the whole slab -- an error of full measure, which is why this one
was wrong under **Lebesgue too**, not only under magnitude.

Interior sampling repairs every piece except `v = 0` itself. That point is
measure-zero under Lebesgue but carries mass 1 under magnitude, and the obvious
fix is unavailable: a step function cannot distinguish `{0}` from `[0, x_0)`, so
neither `[v = 0]` nor `[v > 0]` is expressible as a unary step density.

**Repaired** by giving each axis's weight two components, as the measure itself
has: a density against the continuous part and a separate weight at the anchor
(`PairTerm.atoms`, `mul_split`, and a local `_cum_atom`). A left-continuous
inversion then emits an extra branch owning the anchor alone, where the
predicate's truth is independent of `x` -- `0 <= f(x)` holds for every `x`,
`0 > f(x)` for none -- so that branch either drops the constraint or does not
exist, while the run branches disown the anchor. The two families partition the
surviving variable, so the atom is neither lost nor double-counted.

All 32 combinations of shape, side, strictness and measure are now exact on a
single predicate.

## Finding 8 — a constant bound is the degenerate case of Finding 7, and was missed

Same family as Finding 7, found later and only by building compression. When a
staircase has **no breakpoints at all**, the sampling helper ignored its
`interior` flag and returned the anchor `0` as the representative of
`[0, inf)`. With no breakpoints the single piece *is* `[0, inf)`, whose interior
is `(0, inf)`, and any positive representative does -- returning `0` describes
the anchor and then claims it holds everywhere.

The consequence was concrete: a constant bound `f == 0`, i.e. the constraint
`x_j <= 0`, was read as "every `x` is feasible" for every value of `x_j`. That
invented mass on grid cells where the term is identically zero, so the
compressed numerator was positive where the exact integral was `0`.

Worth recording as its own finding because the general case was already fixed
and tested: the shaped-boundary suite covers hills, valleys and zigzags, but
every generated staircase there has at least one breakpoint. A degenerate
constant slipped through a test suite that was otherwise looking straight at it.

**Localized** by comparing partial states against a partially integrated oracle
after each elimination, which named the exact variable whose elimination first
disagreed. That technique is worth keeping: with 15966 terms in the state, no
amount of staring at the output would have found it.

## Finding 9 — compression in `2d` variables, and the size of `C_d`

Not a defect; the construction the general theorem needs, now built and checked.
The integrand of `\eqref{eq:generalcompression}` is one weighted monotone-pair
term in the `2d` variables `(x, y)`; eliminating the `d` old variables and
dividing by the cell masses gives the conditional expectation. Verified against
the *definition* -- for every cell `C`, `N_F(y)` equals the exact integral of `F`
over `C` -- at `d = 4` (16 cells) and `d = 6` (64 cells), both measures, zero
abstentions. `F_bar` is a term in the same class, and relabelling `y` back to the
physical axes is a pure rename with no residual `x`.

One deliberate departure from the report. It writes cells as `(pi^-, pi^+]`,
left-open; a family of `(a, b]` cells is left-continuous in `y`, so a
right-continuous step cannot say which cell `y` falls in. That is harmless for
`y > 0` and wrong at `y = 0`, where the magnitude atom has mass 1. Cells here
are `[g_m, g_{m+1})` with the top boundary placed beyond the ceiling, so every
cell carries the same decoration, the last closes by the ground set, and the
anchor sits inside the first cell. A representation choice, not a change of
content.

The cost is the finding. From **one** incoming term, the output is 450 terms at
`d = 4` and 78305 at `d = 6`. Canonical merging -- adding the coefficients of
identical terms and dropping what cancels, the only geometry-free merge
available -- recovers 4.3x to 4.5x, and nothing at all in one case. So the
question `\ref{prop:complexity}` leaves open resolves the unhelpful way: the
signed sum does not collapse to a constant number in any useful sense, and
`C_d` compounds across generations. The report's remark that the constants "can
be large if the generic label enumeration is implemented literally" is right and
understated; its advice to keep a specialized base case is a requirement, not a
convenience.

## Finding 10 — the oracle has to be grid-aligned in values, not just breakpoints

Not a report defect; a methodology note that cost real time twice across the two
sprints. The oracle evaluates one representative per cell, which is exact only
if every function in the term is constant on each cell. So a generated staircase
must put **both its breakpoints and its values** on the grid the cells are cut
from: a value off the grid becomes a region boundary *inside* a cell, and the
oracle silently returns the wrong number. Mixing fine staircases with the coarse
`p = 8` grid produced four failures that looked like constructor bugs.

## Observation — Lemma 3 follows from Lemma 2

For closure purposes Lemma 3 is not needed: from `x`'s point of view the mutual
`x`–`y` staircase is just one more bound candidate that happens to depend on
`y`, so two applications of Lemma 2 already give the pair result. The
implementation does exactly this and is exact wherever it applies. Lemma 3 buys
a better constant and a closed form, not extra generality — worth saying,
because it makes the general argument rest on the simpler lemma.

---

## The prompt's seven obligations

| # | obligation | status |
|---|---|---|
| 1 | finite active-bound candidate count | holds: domain bound plus up to two per incident predicate, `O_p(1)` per side |
| 2 | finite label count | holds: product of the two candidate sets; measured max branching 2–14 for `p <= 8` |
| 3 | separable branch formulas | verified at `p = 2..8`, both measures, zero abstentions |
| 4 | closure of branch conditions | **fails as written** (Finding 1), repaired by admitting one-turn boundaries |
| 5 | complement handling | holds: hills are one interval, valleys split into runs, and arbitrary-turn boundaries are handled on the same footing |
| 6 | breakpoint complexity vs primitive count | holds: breakpoints per term *decrease* along the chain (`p=8`: 22 → 20.5 → 20.2 → 15.6 → 12.4 → 6.5 → 1 → 0) |
| 7 | no hidden instance-size dependence in `B_p` | holds at `p = 4` (9x breakpoints, branching saturates at 9–10); consistent but not established at `p = 6`, where the mean still creeps |

---

## What would close the rest

1. **Multi-generation closure above `d = 4`.** One generation is verified at
   `d = 6` and the output is structurally ready to feed back in, but a second
   generation has not been run. With `C_6` around `4 x 10^4` per incoming term,
   this is the measurement most likely to show where the construction stops
   being practical.
2. Compression at `d = 8`, which needs `p = 16`.
3. A real merge. Canonical merging recovers only ~4.5x; whether a
   geometry-aware merge does better is the open question that decides whether
   repeated compression is viable at all.
4. FastHVChan backend behind a `compression_backend` flag.

The honest statement is now: **the elimination calculus is validated to
`p = 12`, the static chain for `d <= 8`, and compression for `d <= 6` -- for a
single generation.** What is not established is that the state stays bounded
across generations, and the measured `C_d` is the reason to doubt it.

Of the eight genuine bugs across these sprints, **five were invisible under
Lebesgue and exposed only by magnitude** — every one an endpoint or atom error
at `0`. That is the strongest practical argument for keeping the magnitude
oracle in the suite, and an argument for the product-measure framing independent
of the complexity result.

Finding 7 is the instructive exception: it was wrong under Lebesgue as well, and
survived because the generators only ever produced one of the four decorated
comparison types. A representation defect hid behind an under-powered generator
rather than behind the choice of measure. Both lessons are needed — vary the
decorations, and test against the measure that gives the anchor mass.

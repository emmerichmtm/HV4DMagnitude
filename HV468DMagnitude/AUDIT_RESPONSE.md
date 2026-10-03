# Audit response

**Paper.** *A Product-Measure and Magnitude Simplification of Chan's Hypervolume
Algorithm*, the generalized version with the section "General fixed-dimensional
pair elimination".

**What this is.** A response to the paper after implementing it. Every claim
below is backed by a run against an exact oracle, and the numbers are
reproducible from the repository. Where I disagree with the text I say what to
write instead.

---

## Verdict

**The construction works. I tried to break it and could not.** No
counterexample was found at any dimension I could reach.

Six statements in the text are, however, wrong or unusable as written: a reader
implementing from the paper will get wrong numbers, not error messages. All six
are repairable and none threatens the result. They are listed as R1–R6.

One further change is mathematical rather than editorial: **Theorem 3 bounds the
wrong quantity.** `C_d` bounds what a single compression *emits*; what the
recursion needs bounded is what the state *holds*. Measured, these differ by a
factor of up to 464. The gap is not a hole in the argument so much as a
substitution of one constant for another, and the replacement looks easier to
prove than the original. This is item **M1** and it is the one I would most
want addressed before submission.

A pattern worth recording in its own right: of eight genuine implementation
bugs, **five were invisible under Lebesgue measure and exposed only by the
magnitude measure**, every one an endpoint or atom error at `0`. That is a
practical argument for the product-measure framing that is independent of the
complexity result, and I think the paper should make it.

---

## What was built and checked

Everything is checked against exact enumeration over a product grid, with the
anchor `{0}` given its own cell so that the magnitude measure is computed
exactly rather than approximated. Every test runs twice, under
`mu = lambda` and under `mu = delta_0 + (1/2) lambda`.

| construction | reached | exact |
|---|---|---|
| one-variable elimination (Lemma 1) | `p = 2..6, 8` fine grid; `p = 9, 10, 12` coarse | yes, both measures, no abstentions |
| non-monotone boundaries | hill, valley, zigzag, mixed; `p = 2..5`, valleys to `p = 7` | yes, both measures |
| compression (eq. 42) | `d = 4` (16 cells), `d = 6` (64 cells) | yes, both measures |
| two generations (Theorem 3, recursive use) | `d = 2, 3, 4, 6` | yes, both measures |

Supporting unit checks: the run decomposition agrees pointwise with brute-force
evaluation on 400 random staircases; and all 32 combinations of boundary shape,
predicate side, strictness and measure are exact on a single predicate in two
variables. That 32-cell matrix is what caught R5, and I would suggest it as the
minimum acceptance test for any other implementation of this calculus.

### Two different `p` per dimension

The paper uses weighted monotone-pair terms in two settings that need different
variable counts, and does not distinguish them. The distinction decides what
"validated at `d = 8`" means:

| setting | variables | `d=4` | `d=6` | `d=8` |
|---|---|---|---|---|
| static prefix chain (Cor. 2) | `p = d` | 4 | 6 | 8 |
| compression numerator (eq. 42) | `p = 2d` | 8 | 12 | 16 |

The static chains behind `T_6 = O~(n^2)` and `T_8 = O~(n^{8/3})` are therefore
tested at their own `p`. Compression is tested at `d = 4` and `d = 6`.
`d = 8` compression needs `p = 16` and was not reached.

---

## R1–R6: the text is wrong as written

### R1. Definition 1's normalization is not closed under the lemmas applied to it

Definition 1 requires every `f_e` to be **monotone**, and justifies "at most one
lower and one upper monotone predicate per pair" by saying same-orientation
predicates merge by min/max.

That merge exists only when the orientations agree. Elimination routinely
produces, for one pair, an upper predicate from a nondecreasing map *and* one
from a nonincreasing map. Their conjunction is

```
xi_j  <|  min( h_up(xi_i), h_down(xi_i) )
```

which has a turning point and is not monotone. So the class is not closed under
the operation both lemmas perform on it, and the proof step "each labelled
region is itself a weighted monotone-pair term" does not follow from Definition
1 as stated.

This matters beyond tidiness because the invariant is load-bearing: the `O(p)`
candidate count in the proof of Lemma 1 *is* "one lower and one upper per
incident pair". With a list per pair the count becomes `O_p(1)` with a worse
constant — still fine for fixed `p`, but not what the proof says.

**Suggested change.** Import the four-dimensional section's own repair.
Definition 4 ("latent one-turn normal form") is already exactly
`min(h_up, h_down)`; let the general definition admit one-turn boundaries, or
carry a list of predicates per pair, and restate the candidate count as
`O_p(1)`.

**Evidence.** With monotone-only predicates the constructor had to decline 42 of
102 instances — 0% at `p = 2`, 40% at `p = 3`, 63% at `p = 4`, 83% at `p = 5`,
**100% at `p = 6`**. Admitting one-turn boundaries took that to zero at every
`p` up to 12.

### R2. Inversion yields a union of intervals, not a one-sided bound

The proof of Lemma 1 says each constraint incident to the eliminated variable
becomes "a lower **or** upper interval bound whose endpoint is a monotone unary
function of one surviving variable". That disjunction is not exhaustive, for
three separate reasons:

1. Even for monotone `f`, the degenerate pieces — every `x` feasible, no `x`
   feasible — belong to the opposite family, so forcing one side per predicate
   misrepresents the rest.
2. Once R1 is repaired, the superlevel set of a **hill** is a genuine middle
   interval with two nontrivial endpoints.
3. The superlevel set of a **valley** is the union of **two disjoint runs**, and
   is not an interval at all.

**Suggested change.** State the inversion as returning a *finite union of
decorated intervals*, each endpoint a unary function of one surviving variable,
plus a `0/1` unary filter for the pieces where a run is absent. Since the
integral over a disjoint union is the sum, a feasible set with `r` runs becomes
`r` signed terms. Two details belong in the text because they are easy to get
wrong: the run count varies **per piece** of the surviving variable, so the
split must be by run *index* with short pieces padded by a degenerate empty
interval; and intersecting two unions of runs gives the union over choices of
one run from each, so the branches **multiply** across the predicates incident
to the eliminated variable. All of this stays `O_p(1)`.

**Evidence.** Interval-valued inversion moved `p = 6` from 0% coverage to 100%.
Run splitting then removed the remaining declines at `p = 9..12`. Measured run
counts never exceed 2 and branch products never exceed 96, so the repair costs
very little.

### R3. EMPTY has no safe numeric encoding, and the paper does not supply one

Definition 1 says EMPTY/FULL are "carried explicitly", but both lemma proofs
then write `L_x(z) < x <= U_x(z)` as though the endpoints were numbers. With an
infinite sentinel the downstream cumulative produces `0 * inf = NaN`, and **NaN
is not annihilated by multiplying with zero** — so a region correctly excluded
by a label indicator still poisons the whole signed sum.

**Suggested change.** Make it part of the normalization that an empty feasible
interval is encoded as a degenerate `[c, c)` — lower closed, upper open at the
same point — which carries mass exactly `0` under both measures. And say
explicitly that `sup (empty set)` must **not** be encoded as `0`: under
magnitude the singleton `{0}` has mass 1 while the empty set has mass 0, and
conflating them is a mass-1 error.

**Evidence.** Before this, every instance at every `p` returned NaN. It went
unnoticed for a full test cycle because the natural tolerance check
`|got - want| > tol` is **false** when `got` is NaN — a NaN silently passes a
float comparison. Any audit of this paper should assert non-NaN separately.

### R4. The label partition is not a partition — ties are unhandled

Both proofs say: "Partition the surviving-variable space according to which
lower endpoint is maximal and which upper endpoint is minimal."

When two candidates coincide on a region of positive measure, both labels hold
and that region is integrated twice. This is not a measure-zero edge case: in
the compression setting distinct predicates routinely invert to the *same*
domain limit, so ties are the common case.

**Suggested change.** Specify a tie-break and carry it into the label
conditions. First index wins is enough: the candidate at index `k` compares
**strictly** against every earlier candidate and non-strictly against every
later one. Worth stating, because the natural implementation — comparing all
others non-strictly — is the wrong one.

**Evidence.** Before the tie-break, integrals came out exactly 2x and 4x too
large.

### R5. `upper != strict` inverts left-continuously, and the anchor needs its own weight

Definition 1 says the decorated comparison carries "endpoint openness
explicitly", but neither proof says what that costs, and it costs more than
bookkeeping. Inverting `v <| f(x)`, the direction of continuity follows the
decoration:

| predicate | feasible set | constant on | right-continuous? |
|---|---|---|---|
| `v < f(x)` | `{f > v}` | `[a, b)` | yes |
| `v >= f(x)` | `{f <= v}` | `[a, b)` | yes |
| `v <= f(x)` | `{f >= v}` | `(a, b]` | **no** |
| `v > f(x)` | `{f < v}` | `(a, b]` | **no** |

The two left-continuous cases — exactly those with `upper != strict` — cannot be
written as right-continuous step functions of `v`. A right-continuous
representation is then wrong only at the jump points, which is harmless at every
jump except `v = 0`: under magnitude the anchor has mass 1, so that single point
is a mass-1 error.

Worse, the obvious repair is not available in the stated class. A step function
cannot distinguish `{0}` from `[0, x_0)`, so **neither `[v = 0]` nor `[v > 0]`
is expressible** as a unary step density, and Definition 1 offers nothing else.

**Suggested change.** Give the weight on each axis two components, as the
measure itself has: a density against the continuous part and a separate weight
at the anchor. Then state that a left-continuous inversion emits an extra branch
owning the anchor alone, where the predicate's truth is independent of `x` —
`0 <= f(x)` holds for every `x` and `0 > f(x)` for none, since every value in
this class is non-negative — so that branch either drops the constraint or does
not exist, while the run branches disown the anchor. The two families then
partition the surviving variable.

**Evidence.** With only `upper = strict = True` predicates, everything passed.
Generating all four decorated comparisons, the two mixed ones failed on a
*single* predicate in two variables: `hill, upper=True, strict=False` gave 4.0
against 2.0 under Lebesgue, and `hill, upper=False, strict=True` gave 0.0
against 2.0. Note this one was wrong **under Lebesgue as well** — it is a
representation defect, not an atom bug, and only its anchor component is
magnitude-specific. It survived because the generators produced one of four
possible decorations.

There is a degenerate sub-case worth a sentence in the text, because it slipped
past a suite that was otherwise looking straight at it: when the staircase has
**no breakpoints at all**, the single piece is `[0, inf)`, whose interior is
`(0, inf)`. A constant bound `f == 0` — the constraint `x_j <= 0` — then gets
read as "every `x` is feasible" for every `x_j` if the anchor is used as the
representative, which invents mass on cells where the term is identically zero.

### R6. The cell family `(pi^-, pi^+]` is left-continuous in `y` and unrepresentable

Equation 42 writes the grid cell as `(pi_i^-(y_i), pi_i^+(y_i)]`, left-open and
right-closed. That is the natural way to write it in prose and the wrong way to
represent it: a family of `(a, b]` cells is left-continuous in `y`, so a
right-continuous step function cannot say which cell `y` falls in. Harmless for
`y > 0`, wrong at `y = 0` where the magnitude atom has mass 1 — the same failure
as R5, now on the `y` side.

**Suggested change.** Use right-open cells `[g_m, g_{m+1})` with the top
boundary placed beyond the domain ceiling, so every cell carries the same
decoration, the last one closes by the ground set rather than by a special case,
and the anchor sits inside the first cell. Both boundary maps are then
right-continuous step functions of `y`. This is a representation choice, not a
change of content: any partition into intervals gives a valid conditioning, and
this is the one the state class can carry.

---

## M1: Theorem 3 bounds the wrong quantity

This is the one item I would call substantive.

Theorem 3 gives `S' <= C_d S` with `C_d` depending only on `d`, and Proposition
1 concedes that "nothing in the current argument proves that a signed sum of
distinct rank-one/one-turn primitives can always be canonically merged back to a
constant number". Both statements are about `C_d`.

**The concession is justified.** Measured, from one incoming term: 450 output
terms at `d = 4` and 78305 at `d = 6`. Canonical merging — adding the
coefficients of identical terms and dropping what cancels, the only
geometry-free merge available — recovers only 4.3x to 4.5x, and nothing at all
in one case. `C_d` really is large and really does not merge away. The remark
that the constants "can be large if the generic label enumeration is implemented
literally" is correct and, on this evidence, understated; the advice to keep a
specialized base case reads as a requirement rather than a convenience.

**But `C_d` is not what governs repeated compression.** Across two generations
on nested grids — `G2` coarser than `G1`, so each `G2` cell is a union of `G1`
cells and the tower property applies — the merged state holds steady, and at
`d = 6` it shrinks:

| d | measure | gen 1 raw → merged | gen 2 raw → merged | merged growth |
|---|---|---|---|---|
| 3 | magnitude | 456 → 274 | 4038 → 468 | 1.7x |
| 4 | magnitude | 3985 → 2805 | 193468 → 11341 | 4.0x |
| 4 | Lebesgue | 1225 → 1075 | 54131 → 114 | 0.1x |
| 6 | Lebesgue | 8330 → 3120 | 602140 → **1296** | **0.42x** |
| 6 | magnitude | 31080 → 7110 | 978120 → **4473** | **0.63x** |

The raw emission reaches **464x** the merged size. So `C_d` bounds what one
compression emits; the state holds something much smaller, and the paper does
not separate the two.

**The reason is structural, and I think it is the missing lemma.** After
compression on `G2`, every function in the state is `G2`-measurable. So the
number of **distinct** terms the state can contain is capped by the grid's
resolution, while the number the construction emits is not. And `G2`, being
coarser than `G1`, admits *fewer* distinct terms — which is why the `d = 6`
state shrinks rather than merely holding steady.

**Suggested change.** Separate the two quantities explicitly. Keep `C_d` as a
bound on what one compression emits, and add the quantity the recursion actually
needs: a bound on the number of distinct `G`-measurable weighted monotone-pair
terms at a given grid resolution. A bound on *that* looks considerably more
tractable than a bound on `C_d`, and it is what Theorem 3 requires for recursive
use. Say also that merging is part of the construction rather than an
implementation detail, and that it must be applied **incrementally** — folding
each contribution in as it is produced — because batching it to the end bounds
memory by the raw count and makes a second generation unreachable in practice.

A corollary worth a remark, because the naive projection is badly wrong (I made
it myself before measuring): **the expensive compression is the first one.** A
generation-1 term at `d = 6` compresses to 27–81 raw terms rather than 78305,
because its predicates are already grid-aligned and coarse. Projecting a second
generation by multiplying the state size by the first generation's `C_d` gives
roughly `3 x 10^9` terms and an estimate of weeks; the measured run is minutes.

None of this touches the asymptotics. `C_d` is independent of instance size,
which is what the recurrence needs, and `T_d(n) = O~_d(n^{d/3})` is unaffected.
M1 is about which constant the theorem should be stated in terms of.

---

## Gaps: claims that measure as true but are not argued

### G1. `B_p`'s independence of instance size is asserted, not proved

The remark after Lemma 2 says the purpose of `B_p` is "only to make explicit
that the branching factor is independent of `n`, `m`, and the number of
staircase breakpoints". The proof establishes that the **label count** depends
only on `p`; it never bounds the number of *emitted terms per label*
independently of the staircase resolution. Those are different quantities, and
the second is what the recurrence consumes.

It does appear to be true. Holding `p` and the predicate graph fixed and varying
only resolution:

| `p = 4` term breakpoints | 6 | 10 | 13 | 22 | 29 | 39 | 53 |
|---|---|---|---|---|---|---|---|
| max branching | 2 | 8 | 10 | 8 | 10 | 9 | 9 |

A nine-fold increase in breakpoints leaves branching at 9–10. At `p = 6` the
maximum also levels off (18, 19, 20) but the mean is still creeping upward at
the largest resolution the oracle allows, so I would describe `p = 6` as
consistent with independence rather than established.

**Suggested change.** Add the missing step: bound emitted terms per label by the
number of branches of the staircase-prefix formula (a constant) times the number
of predicate-list entries per pair (`O_p(1)` after R1), and note that the
feasibility filters are absorbed into unary densities and so emit no terms.

### G2. The `O_p(N)` breakpoint claim is per-elimination; the chain is not addressed

Both lemmas claim output breakpoint complexity `O_p(N)`. Over a chain of depth
`p`, monotone inverses are composed repeatedly, and composition could in
principle refine breakpoints multiplicatively — which would make per-term cost
grow with depth even at bounded branching.

Measured breakpoints **per surviving term** after each elimination:

```
p = 6:  17.0 -> 14.3 -> 11.5 -> 7.2 -> 2.0 -> 0
p = 8:  22.0 -> 20.5 -> 20.2 -> 15.6 -> 12.4 -> 6.5 -> 1.0 -> 0
```

Monotone decrease throughout, so the concern does not materialize. The paper
should say why rather than leaving the reader to worry about it.

---

## Editorial

### E1. Lemma 2 is not needed for closure

For closure purposes the pair lemma is redundant: from `x`'s point of view the
mutual `x`–`y` staircase is just one more bound candidate that happens to depend
on `y`, so two applications of Lemma 1 already give the pair result. The
implementation does exactly this and is exact wherever it applies. Lemma 2 buys
a better constant and a closed form, not extra generality.

Worth saying, because it makes the general theorem rest on the simpler of the
two lemmas — a strengthening of the paper's position, not a weakening.

### E2. The ceiling clamp belongs in the normalization, and its scope matters

`P(A) = integral of Y(f(x)) dmu` diverges when `f == +inf`, which is what a pair
carrying no constraint looks like; the separable split then evaluates as
`inf - inf` while the quantity it represents is finite. In the general setting
the repair is exact, because in the compression integrand every `x_i` is
confined to a grid interval, so a vacuous constraint *is* the cell ceiling.

Two things to add. First, that this is part of the normalization rather than an
implementation detail — the 4-D solver returned NaN on real instances before it
was added. Second, and easy to get backwards: the clamp is valid **only at the
measure step**. Clamping the candidate before the label comparisons is wrong,
because there `+inf` correctly means "weaker than every other bound", and
collapsing it onto the ceiling makes a vacuous bound tie with the ground set and
claim regions it does not own. I tried the stronger clamp and it broke `p = 4`
and `p = 5`.

### E3. The validation-status paragraph is out of date, but should stay narrow

"`6D` and `8D` implementations have not yet been written or experimentally
audited" is no longer accurate. The replacement should be narrower than it is
tempting to make it:

> The elimination calculus is validated against exact oracles to `p = 12`,
> covering the static chain for `d <= 8`. Compression is validated at `d = 4`
> and `d = 6`, through two generations on nested grids. `d = 8` compression,
> which requires `p = 16`, and chains deeper than two generations remain
> untested.

---

## What remains untested

* **`d = 8` compression**, which needs `p = 16`.
* **More than two generations.** Two behave well; nothing here shows what
  happens after ten, and that is where M1's missing bound would earn its keep.
* **Non-monotone exactness above `p = 7`**, and zigzag cost above `p = 6`.
  Above `p = 8` the oracle grid is necessarily coarse enough that staircases
  cannot form a valley, so those rows do not exercise run splitting.
* **A merge cheaper than the canonical one**, which is quadratic in state size
  as written. It turns out to be *sufficient* for keeping the state bounded
  across generations, so the open question is now its cost rather than its
  power.
* **`C_d` as a function of `d`** from more than two data points. 450 at `d = 4`
  and 78305 at `d = 6` is a ratio of 174 across one step — consistent with
  `2^{O(d^3)}`, establishing nothing.

---

## Suggested edit checklist

| # | where | change |
|---|---|---|
| R1 | Definition 1 | admit one-turn boundaries or a predicate list per pair; restate the candidate count as `O_p(1)` |
| R2 | proof of Lemma 1 | inversion returns a union of decorated intervals; split by run index; branches multiply across incident predicates |
| R3 | Definition 1 | encode an empty interval as degenerate `[c, c)`; state that `sup (empty)` is not `0` |
| R4 | proofs of Lemmas 1 and 2 | specify a first-index tie-break in the label partition |
| R5 | Definition 1 | two-component weights (anchor and continuous); anchor branch for left-continuous inversions |
| R6 | equation 42 | right-open cells with the top boundary beyond the ceiling |
| M1 | Theorem 3, Proposition 1 | separate "terms emitted" from "terms held"; add the grid-resolution bound; make merging part of the construction, applied incrementally |
| G1 | remark after Lemma 2 | bound emitted terms per label, not just the label count |
| G2 | Lemmas 1 and 2 | state why breakpoints do not compound along the chain |
| E1 | after Lemma 2 | note that Lemma 1 twice suffices for closure |
| E2 | normalization | ceiling clamp, and that it applies at the measure step only |
| E3 | validation status | replace with the narrower statement above |

---

## Reproducing

```bash
python tests/test_elimination.py    # elimination, p = 2..12
python tests/test_valley.py         # non-monotone boundaries, run decomposition
python tests/test_compression.py    # compression at d = 4 and d = 6
python tests/test_generations.py    # two generations on nested grids
python complexity_probe.py          # branching vs instance size and dimension
```

Supporting detail: `notes/PAIR_ELIMINATION_AUDIT.md` (findings in the order
they were found, including the four bugs the magnitude oracle caught and
Lebesgue did not), `notes/REPORT_PROBLEMS.md` (the same items graded by
severity), `PAIR_ELIMINATION_COMPLEXITY.md` (all measured tables).

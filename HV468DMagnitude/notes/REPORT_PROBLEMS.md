# Problems with the pair-elimination report

Report: *A Product-Measure and Magnitude Simplification of Chan's Hypervolume
Algorithm*, section **General fixed-dimensional pair elimination**
(`paper/magnitude_chan_pair_elimination_report.tex`, `\label{sec:generalpair}`).

Nothing below is a counterexample to the main theorem. I tried to break it and
could not: the constructor now agrees with exact oracles at `p = 2..6, 8, 9, 10`
and `12` under both Lebesgue and magnitude product measures, with zero
abstentions, and with non-monotone boundaries of every shape. The
problems are of three kinds, and they are not equally serious:

* **Blocking** — the statement as written is false or unusable, and a reader
  implementing from the text will get a wrong answer. Five of these.
* **Gap** — the claim appears true (I measured it) but the proof does not
  contain the argument. Three.
* **Editorial** — correct but misleading, understated, or out of date. Three.

Every blocking item was found by implementing, and each one produced a
*silently* wrong number first, not an error. That is the pattern worth noting:
the class is forgiving enough that a wrong implementation looks like a working
one until the dimension is high enough to expose it.

---

## Blocking

### B1. Definition 1's normalization is not closed under the lemmas applied to it

`\label{def:pairterm}` requires every `f_e` to be a **monotone** step map, and
justifies the normal form — "for each unordered variable pair there is at most
one lower and one upper monotone predicate" — with "since same-orientation
predicates are merged by min/max".

That merge exists only when the orientations agree. Elimination routinely
produces, for one pair, an upper predicate from a nondecreasing map *and* one
from a nonincreasing map. Their conjunction is

```
xi_j  <|  min( h_up(xi_i), h_down(xi_i) )
```

which has a turning point and is not monotone. So the class is not closed under
the operation both lemmas perform on it, and the proof step "each labelled
region is itself a weighted monotone-pair term" does not follow from
Definition 1.

This matters beyond tidiness, because the invariant is **load-bearing**: the
`O(p)` candidate count in the proof of Lemma 2 is exactly "one lower and one
upper per incident pair". If a pair carries a list instead, the count becomes
`O_p(1)` with a worse constant rather than `O(p)`, which is still fine for fixed
`p` but is not what the proof says.

**Evidence.** With monotone-only predicates the constructor abstained on 42 of
102 instances — 0% abstention at `p = 2`, 40% at `p = 3`, 63% at `p = 4`, 83% at
`p = 5`, **100% at `p = 6`**. Admitting one-turn predicates took that to zero
abstentions at every `p` up to 8.

**Fix.** Import the four-dimensional section's own repair. `\label{def:latent}`
("latent one-turn normal form") is already exactly `min(h_up, h_down)`; the
general definition should admit one-turn boundaries, or keep a list of
predicates per pair, and restate the candidate count as `O_p(1)`. Either keeps
the asymptotics. This is a repair, not a refutation.

### B2. Decorated inversion is interval-valued, not one-sided

The proof of Lemma 2 (`\label{lem:oneelim}`) says every constraint incident to
the eliminated variable `x` becomes "a lower **or** upper interval bound whose
endpoint is a monotone unary function of one surviving variable".

That disjunction is not exhaustive, for two independent reasons:

1. Even for monotone `f`, the *degenerate* pieces belong to the opposite
   family. Inverting `xi_j <| f(x)` with `f` nondecreasing gives a suffix — a
   lower bound — but on the pieces where every `x` is feasible the natural
   encoding is an upper bound of `+inf`, and on the pieces where none is, a
   lower bound of `+inf`. Forcing one side per predicate misrepresents the
   others.
2. Once B1 is repaired and one-turn maps are in the class, the superlevel set of
   a hill is a genuine middle interval with two nontrivial endpoints.

**Fix.** State the inversion as returning a *finite union of decorated
intervals* whose endpoints are each unary functions of one surviving variable,
plus a `0/1` unary filter recording the pieces where a run is absent. The
feasible set is a single interval for a monotone map and for a hill, but the
superlevel set of a *valley* is two disjoint runs. Since the integral over a
disjoint union is the sum, each run becomes its own signed term; the number of
runs varies per piece, so the split is by run *index*, and it multiplies across
the predicates incident to the eliminated variable (intersecting two unions of
runs gives the union over choices of one run from each). All of that stays
`O_p(1)`.

**Evidence.** Interval-valued inversion moved `p = 6` from 0% coverage to 100%;
run splitting then removed the remaining abstentions at `p = 9..12`, with
`p = 9, 10, 12` exact against a coarse oracle. Measured run counts stay at 2 and
branch products at 96 or below.

### B5. Inverting a non-strict upper or strict lower predicate is left-continuous, and the anchor atom needs its own weight

Definition 1 says the decorated comparison carries "EMPTY/FULL and endpoint
openness explicitly", but neither lemma proof says what that costs. It costs
more than bookkeeping.

Invert `v <| f(x)` for `v` a surviving variable. The feasible `x`-set moves as
`v` crosses the values `f` attains, and the direction of continuity depends on
the decoration:

| predicate | feasible set | constant on | representable right-continuously |
|---|---|---|---|
| `v < f(x)` | `{f > v}` | `[a, b)` | yes |
| `v >= f(x)` | `{f <= v}` | `[a, b)` | yes |
| `v <= f(x)` | `{f >= v}` | `(a, b]` | **no** |
| `v > f(x)` | `{f < v}` | `(a, b]` | **no** |

The two left-continuous cases — exactly those with `upper != strict` — cannot be
written as right-continuous step functions of `v`. A right-continuous
representation is wrong only at the jump points, which is harmless for every
jump except the one at `v = 0`: under magnitude the anchor has mass 1, so that
single point is a mass-1 error, while under Lebesgue it is invisible.

Worse, the natural repair is not available in the stated class. A step function
cannot distinguish `{0}` from `[0, x_0)`, so **neither `[v = 0]` nor `[v > 0]` is
expressible** as a unary step density, and Definition 1 offers nothing else.

**Fix.** Give the weight on each axis two components, as the measure itself has:
a density against the continuous part and a separate weight at the anchor. Then
state in the normalization that a left-continuous inversion emits an extra
branch owning the anchor alone, where the predicate's truth is independent of
`x` — `0 <= f(x)` holds for every `x` and `0 > f(x)` for none, since every value
in this class is non-negative — so that branch either drops the constraint or
does not exist, and the run branches disown the anchor. The two families then
partition the surviving variable.

**Evidence.** With only `upper = strict = True` predicates, everything passed.
Generating all four combinations, the two mixed ones failed on a *single*
predicate in two variables — `hill, upper=True, strict=False` gave 4.0 against
2.0 under Lebesgue, and `hill, upper=False, strict=True` gave 0.0 against 2.0.
The Lebesgue errors came from the interior-sampling half of the bug; the
magnitude errors needed the split weight. After both, all 32 combinations of
shape, side, strictness and measure are exact.

This is the fifth bug in this class that Lebesgue hid and magnitude exposed, and
the first that Lebesgue *also* got wrong — worth noting, because it means the
"magnitude catches endpoint errors" lesson understates the problem: the
representation defect was real under both measures, and only its anchor
component was magnitude-specific.

### B3. EMPTY has no safe numeric encoding, and the text does not supply one

Definition 1 says EMPTY/FULL are "carried explicitly", but both lemma proofs
then write `L_x(z) < x <= U_x(z)` as though the endpoints were numbers. If EMPTY
is encoded with an infinite sentinel, the downstream cumulative produces
`0 * inf = NaN`, and **NaN is not annihilated by multiplying with zero** — so a
region correctly excluded by a label indicator still poisons the whole signed
sum.

**Fix.** Make it part of the normalization that an empty feasible interval is
encoded as a degenerate `[c, c)` — lower closed, upper open at the same point —
which carries mass exactly `0` under Lebesgue and under magnitude alike, so the
anchor atom is neither kept nor double-counted. Say explicitly that `sup ∅` must
not be encoded as `0`: under magnitude the singleton `{0}` has mass 1 while the
empty set has mass 0, and conflating them is a mass-1 error.

**Evidence.** Until this was fixed, *every* instance at every `p` returned NaN.
It went unnoticed for a full test cycle because the natural tolerance check
`abs(got - want) > tol` is `False` when `got` is NaN — a NaN silently passes a
float comparison. Any audit of this report should assert non-NaN separately.

### B4. The label partition is not a partition — ties are unhandled

Both proofs say: "Partition the surviving-variable space according to which
lower endpoint is maximal and which upper endpoint is minimal."

When two candidates coincide on a region of positive measure, both labels hold
and that region is integrated twice. This is not a measure-zero edge case: in
the compression setting distinct predicates routinely invert to the *same*
domain limit, so ties are the common case, not a rarity.

**Fix.** Specify a tie-break and carry it into the label conditions. First
index wins is enough: the candidate at index `k` compares **strictly** against
every earlier candidate and non-strictly against every later one. The report
should state this, because the natural implementation — comparing all others
non-strictly — is the wrong one.

**Evidence.** Before the tie-break, integrals came out exactly 2x and 4x too
large.

---

## Gaps — true as measured, not argued

### G1. `B_p`'s independence of instance size is asserted, not proved

The remark after Lemma 3 says the purpose of `B_p` is "only to make explicit
that the branching factor is independent of `n`, `m`, and the number of
staircase breakpoints". But the proof establishes that the **label count**
depends only on `p`; it never bounds the number of *emitted terms per label*
independently of the staircase resolution. Those are different quantities, and
the second is the one the recurrence consumes.

It does appear to be true. Holding `p` and the predicate graph fixed and varying
only staircase resolution (`complexity_probe.py`):

| p | breakpoints | 6 | 10 | 13 | 22 | 29 | 39 | 53 |
|---|---|---|---|---|---|---|---|---|
| 4 | max B | 2 | 8 | 10 | 8 | 10 | 9 | 9 |
| 4 | mean B | 1.3 | 3.0 | 3.2 | 4.6 | 5.6 | 4.0 | 4.8 |

| p | breakpoints | 11 | 19 | 27 | 42 | 52 |
|---|---|---|---|---|---|---|
| 6 | max B | 3 | 6 | 18 | 19 | 16 |
| 6 | mean B | 1.3 | 2.7 | 9.1 | 10.7 | 12.8 |

Max branching **saturates** — a 9x increase in breakpoints at `p = 4` leaves it
at 9–10. That is the behaviour the lemma needs. But note the mean at `p = 6` is
still creeping upward at the largest resolution tested, so the saturation is
clean at `p = 4` and only probable at `p = 6`.

**Fix.** Add the missing step: bound emitted terms per label by the number of
branches of the staircase-prefix formula (a constant) times the number of
predicate-list entries per pair (`O_p(1)` after B1), and note that the
feasibility filters are absorbed into unary densities and so emit no terms.

### G2. The `O_p(N)` breakpoint claim is per-elimination; the chain is not addressed

Both lemmas claim output breakpoint complexity `O_p(N)`. Over a chain of depth
`p`, monotone inverses are composed repeatedly, and composition could in
principle refine breakpoints multiplicatively — which would make the per-term
cost grow with depth even at bounded branching, and hollow out the
fixed-dimensional claim.

Measured breakpoints **per surviving term** after each elimination:

```
p = 6:  17.0 -> 14.3 -> 11.5 -> 7.2 -> 2.0 -> 0
p = 8:  22.0 -> 20.5 -> 20.2 -> 15.6 -> 12.4 -> 6.5 -> 1.0 -> 0
```

Monotone decrease throughout, so the concern does not materialize. The report
should say why, not leave it to the reader to worry about.

### G3. `C_d` is left implicit where it compounds

`\label{thm:generalsimplification}` gives `S' <= C_d S` with `C_d` a product of
`B_p` along the chain, `p <= 2d`. With the remark's own `B_p = 2^{O(p^2)}`, this
is `C_d = 2^{O(d^3)}`. The report should write that down, because
`\label{prop:complexity}` already concedes that "nothing in the current argument
proves that a signed sum of distinct rank-one/one-turn primitives can always be
canonically merged back to a constant number" — so `C_d` compounds across
compression *generations*, not only within one, and at `2^{O(d^3)}` that caveat
is doing much more work in the general theorem than in the 4-D one.

Peak state size at fixed resolution, as `p` grows: 2, 6, 6, 14, 26, 89, 499 for
`p = 2..8`. The branching constant is modest; the *state* is what grows.

---

## Editorial

### E1. Lemma 3 is not needed for closure

For closure purposes Lemma 3 (`\label{lem:pairelim}`) is redundant: from `x`'s
point of view the mutual `x`–`y` staircase is just one more bound candidate that
happens to depend on `y`, so two applications of Lemma 2 already give the pair
result. The implementation does exactly this and is exact wherever it applies.

Lemma 3 buys a better constant and a closed form (the two-branch formula), not
extra generality. Worth saying, because it makes the general theorem rest on the
simpler of the two lemmas — which is a strengthening of the paper's position,
not a weakening.

### E2. Finding 2 from the 4-D audit is not carried into the general normalization

`P(A) = ∫ Y(f(x)) dμ` diverges when `f ≡ +inf`, which is what a pair carrying no
constraint looks like; the separable split then evaluates as `inf - inf` while
the quantity it represents is finite. In the general setting the repair is
available and exact, because in the compression integrand every `x_i` is
confined to a grid interval, so a vacuous constraint *is* the cell ceiling.

Two things the text should add. First, that the clamp is part of the
normalization rather than an implementation detail — the 4-D solver returned NaN
on real instances before it was added. Second, and this is the part that is easy
to get backwards: the clamp is valid **only at the measure step**. Clamping the
candidate before the label comparisons is wrong, because there `+inf` correctly
means "weaker than every other bound", and collapsing it onto the ceiling makes
a vacuous bound tie with the ground set and claim regions it does not own. I
tried the stronger clamp first and it broke `p = 4` and `p = 5`.

### E3. The validation-status paragraph is out of date

"`6D` and `8D` implementations have not yet been written or experimentally
audited" is no longer accurate — pair elimination is now implemented and exact
against product-grid oracles at `p = 2,3,4,5,6,8` under both measures.

But the replacement sentence should be narrower than it is tempting to make it.
What is validated is **one- and two-variable elimination from a weighted
monotone-pair term**. What is *not* validated is compression in the `2d`
variables `(x, y)` of `\eqref{eq:generalcompression}`, multi-generation closure
above 4-D, or the `S' <= C_d S` bound of the general theorem. The useful
distinction for a reader is "the elimination calculus is tested to `p = 8`; the
compression that uses it is tested only at `d = 4`".

---

## Summary table

| # | kind | item | status in code |
|---|---|---|---|
| B1 | blocking | monotone-only class not closed; invariant is load-bearing | repaired (one-turn admitted) |
| B2 | blocking | inversion is a union of runs, not a one-sided bound | repaired (run splitting) |
| B3 | blocking | EMPTY has no safe encoding; NaN survives zero-multiplication | repaired (degenerate `[c,c)`) |
| B4 | blocking | label partition double-counts ties | repaired (first-winner) |
| B5 | blocking | `upper != strict` inverts left-continuously; anchor needs its own weight | repaired (split weight + anchor branch) |
| G1 | gap | `B_p` independent of `N` asserted, not proved | measured, saturates |
| G2 | gap | breakpoint growth along the chain not addressed | measured, decreases |
| G3 | gap | `C_d = 2^{O(d^3)}` left implicit where it compounds | measured peak state |
| E1 | editorial | Lemma 3 redundant for closure | implementation uses Lemma 2 twice |
| E2 | editorial | ceiling clamp missing from normalization, and its scope | repaired, scope documented |
| E3 | editorial | validation status stale, but should stay narrow | — |

One more observation that belongs in the paper rather than in this list. Of the
seven genuine bugs across these sprints, **five were invisible under Lebesgue
and exposed only by the magnitude measure**, every one an endpoint or atom error
at `0`. Keeping the magnitude oracle in the suite is the single most effective
check on this class of construction, and that is an argument for the
product-measure framing on its own terms, independent of the complexity result.

The counterpart is worth saying too, because it cuts the other way: B5 was
wrong under Lebesgue *as well*, and it went unnoticed because the test
generators only ever produced one of the four decorated comparison types.
A measure-agnostic representation defect hid behind an under-powered generator,
not behind the choice of measure. Both lessons are needed: vary the decorations,
and check against the measure that gives the anchor mass.

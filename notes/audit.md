# Audit of the product-measure / magnitude simplification of Chan's 4-D algorithm

Report audited: *A Product-Measure and Magnitude Simplification of Chan's
Four-Dimensional Hypervolume Algorithm* (2 October 2026).

Everything below was obtained by executing the claims against exact oracles, in
the order the accompanying task specification prescribes. Each identity was
implemented independently of its brute-force twin, and the brute force
enumerates the product grid cell by cell with the anchor atom treated as its own
cell, so it is exact rather than approximate.

Short version: **the mathematics I was able to test holds.** No counterexample
to any stated lemma survived. Three things need correcting or sharpening in the
write-up, and one substantial claim remains untested — see
[What is not verified](#what-is-not-verified), which is the honest limit of this
exercise.

---

## What was verified

| Claim | Where | Evidence |
|---|---|---|
| Magnitude as a product measure; `Mag = Σ_S 2^-\|S\| HV(π_S D)`; the `HV₄` recovery formula | eq. (136), (146), (160) | 8 random 4-D instances, both routes agree |
| Threshold-antiderivative lemma, both orientations | Lemma 4.? (`lem:threshold`) | 480 randomized + 240 exhaustive integer cases, both measures, including plateaus, queries at breakpoints, empty feasible intervals and the anchor atom |
| Two-branch staircase formula | eq. (261) | 12,000 cases over all four measure pairings and both inequality conventions |
| Meet contraction | eq. (225) | 80 cases against explicit pairwise-meet evaluation |
| Boundary class `min(h↑, h↓)` closed under repeated absorption | §"Evaluating one primitive" | 9,600 evaluations, 3 successive masks per trial; orientations preserved |
| **Separability of the source mass on each fixed region** | Thm. `thm:prefix` | 11,250 evaluations against a direct 2-D sweep, both measures |
| Region count is dimension-bounded | Thm. `thm:prefix` | max 14 distinct labels observed against a structural maximum of 18; independent of staircase complexity |
| Four-prefix integral `K(a,b,c,d)` | eq. (337) | 144 queries against exact 4-D enumeration, both measures |
| Tower property / compression closure | eq. (457) | 192 grid cells after absorb-and-coarsen, both measures |
| Model corresponds to real hypervolume | §"Mapping to FastHVChan" | 40 instances agree with the independent `FastHVChan` implementation |

The separability result is the one that matters most. It is the actual content
of the four-prefix normal form: within each fixed region of the target plane the
source mass is built from four *unary* records `X, Y, P, Q` evaluated at
arguments that each depend on at most one target coordinate, and the prefix
arguments `a, b` enter only as anchored rectangles. I tested this pointwise
rather than inferring it, and it held in every case.

---

## Findings that should change the write-up

### 1. The staircase inequality should be strict, and this is not cosmetic

The report writes the easy residual as

> `F(x₁,x₂,x₃,x₄) = Π_{i<j} [x_j ≤ f_ij(x_i)]`  (eq. 182)

and states the first branch of the two-branch formula under `B ≤ f(A)`.

The region actually left uncovered by a union of grounded quadrants
`{x_i ≥ a_t, x_j ≥ b_t}` is `{x_j < f_ij(x_i)}`, with `f(x) = min{b_t : a_t ≤ x}`
— a strict inequality. With the strict form, branch 1 holds under `B < f(A)`.

Under Lebesgue measure the two conventions differ on a null set and nothing
happens. Under the **magnitude** measure — the one the note is built around —
they differ by a positive amount whenever `f` attains the value 0, because the
atom at `x_j = 0` then carries mass. On a two-quadrant instance I measured

```
max |closed convention − strict convention| = 0.500
```

Since the note's thesis is that magnitude is the natural measure here, eq. (182)
and the branch condition in eq. (261) should both be restated strictly. The
report's own "Endpoint convention" paragraph anticipates exactly this hazard; it
simply has not been applied to the two displayed equations.

### 2. `sup ∅` must be an exclusive bound, not 0

`τ(B) = sup{x : f(x) ≥ B}` is used throughout as if it were an inclusive prefix
endpoint. When the set is empty the natural coding `τ = 0` is wrong under
magnitude: it silently contributes the anchor atom, mass 1, to an empty set. The
inverse has to carry explicit "the prefix is `[0, θ)`" semantics, with `θ = 0`
meaning *empty* and therefore distinct from the set `{0}`.

This is not hypothetical. The single implementation error that survived into my
first run was of exactly this shape: I treated the source prefix `[0, 0]` as
empty when it is the anchor atom, which returned `2.0` where the truth is `1.5`
and produced a wrong four-prefix value. It is invisible under Lebesgue and
caught immediately under magnitude. It is pinned as a regression test
(`Phase 6d`).

Recommendation: state the convention once, in the style of
`prefix_bound_ge(f, u) -> θ, meaning {x : f(x) ≥ u} = [0, θ)`, and make every
cumulative evaluation take an explicit open/closed flag.

### 3. The measurability hypothesis in the tower argument is load-bearing

The report justifies discarding fine geometry with

> "Since `B` was defined using coordinates already present in the parent hard
> grid, it is 𝒢-measurable."

That is correct, but the implementation consequence deserves to be stated
explicitly, because it constrains *which* events may be discarded. I built an
instance where the easy state varies at `x₃ = 1.5`, compressed on a grid that
does not contain `1.5`, and then applied a descendant mask cutting at `1.5`:

```
max |compressed − exact| with an off-grid mask = 0.250
```

So the identity fails exactly when the hypothesis does. The practical rule is
stronger than "discard the fine events after compression": the compression grid
must contain every coordinate that **any** descendant could later turn into a
mask — which includes the cut coordinates, not only the surviving hard
orthants' coordinates. In Chan's schedule those are all hard-orthant face
coordinates, so the hypothesis does hold; but an implementation that coarsens
more aggressively than the hard grid would silently lose correctness, and the
failure is not detectable from the compressed state alone.

### 4. Minor: monotone images do introduce new coordinate values

The shared-record argument in Thm. `thm:prefix` says that minima "do not create
new numerical coordinates" because `min(x,y)` is one of `x, y`. True. But the
same paragraph relies on monotone images `ρ₃(u)`, `σ₃(u)`, … , and those *are*
new values — they are breakpoints of the staircases, not hard coordinates. The
`O(m)` bound survives (there are `O(m)` such images), but the sentence as
written proves less than it needs to. Worth one extra clause.

---

## Follow-up: is the remaining gap repairable?

Short answer: **the three findings above are already repaired, and the untested
claim now has strong supporting evidence.** The gap is narrower than the
section below originally stated, so that section is kept but qualified.

The worry was that one compression might multiply the primitive count by a
factor depending on the *instance* (`N`, `m`) rather than on the dimension. Note
first that the bar is low: because a root-to-leaf path sees only `O(1)`
compressions (the recursion unfolds to depth `log_{r^3} N = 1/(3δ)` for
`r = N^δ`), *any* blow-up function depending on the dimension and the incoming
count alone composes to a constant — even a squaring. The construction fails
only if the factor depends on the problem size.

Structurally it cannot. A state entering its second compression carries, per
physical coordinate, threshold constraints from the previous primitive's
arguments `C(x) = min_i c_i(x_i)`, `D(x) = min_i d_i(x_i)`. The admissible source
prefix is then capped by a minimum of a *fixed* list of candidates:

```
x1 :  a, rho3(ut), rho4(vt), theta_{c1}(u), theta_{d1}(v)
x2 :  b, sig3(ut), sig4(vt), theta_{c2}(u), theta_{d2}(v)
```

Five each, whatever the staircases look like, so the label space is bounded by
`5 x 5 x 2 = 50` — dimension-only by construction. What is *not* automatic is
that each labelled region is still separable once the latent coupling is
present; that is the actual content of `lem:parametricclosure`.

`closure_probe.py` measures both, growing instance size by 16x with the probe
density scaled accordingly (a fixed probe count would make flatness an artefact
of undersampling):

| breakpoints/staircase | `N` | max regions, Lebesgue | max regions, magnitude |
|---|---|---|---|
| 2 | 12 | 15 | 19 |
| 4 | 24 | 16 | 17 |
| 8 | 48 | 20 | 26 |
| 16 | 96 | 20 | 24 |
| 32 | 192 | 26 | 22 |

The count stays far below the structural bound. It is *not* flat — pushed over a
64x range in `N` (12 to 768) it creeps 19, 16, 19, 22, 21, 26, 27 — but that is
saturation, not growth: the measured quantity is how many of the finitely many
possible labels a finite sample happens to realise, and both the probe density
and the geometric variety rise with the instance. A 64x increase in `N` buys a
1.4x increase in the count against a hard ceiling of 50. Proportional growth
would have been the falsification; this is not it. (An early version of the
check tested flatness rather than boundedness and flagged saturation as a
failure — worth recording, since the distinction is the whole point.)

**Separability survived the composition in every one of 1,440+ checks**, across
both measures and every instance size.

So the composition does close in the way the lemma needs, and `C_4` is
dimension-only.

### Multi-generation compression (`generations.py`)

The remaining question was whether this survives *iteration*: a state that is
compressed, absorbed into, and compressed again, several times. That experiment
is now done, and it is the one the first draft of this audit listed as missing.

Four generations of (absorb a newly-easy mask → compress onto a strictly coarser
grid), verified at every generation against an exact recomputation from the
**original fine geometry** with all masks so far included:

```
worst |compressed - exact| over 456 cells, Lebesgue : 0.000e+00
worst |compressed - exact| over 456 cells, magnitude: 4.441e-15   (float noise)
```

So an error introduced at generation 1 and masked at generation 2 would have
shown up, and none did. The per-generation blow-up factor is also stable rather
than creeping:

| generation | staircase complexity | regions `R` (Lebesgue) | regions `R` (magnitude) | separability |
|---|---|---|---|---|
| 1 | 4 | 18 | 19 | 288/288 |
| 2 | 6 | 21 | 21 | 288/288 |
| 3 | 8 | 20 | 17 | 288/288 |

The reason `R` is constant across generations is structural, and worth stating
because it is what makes the induction work: **after a compression the state has
the same shape it had before** — physical six-staircase geometry multiplied by
one latent pair from the primitive's arguments `C(x), D(x)`. Generation `g+1`
therefore hands the composer exactly the configuration generation `g` saw, so
the generation map is a fixed point and `p_g ≤ p_0 · R^g` with `R` bounded by
the dimension. Combined with `O(1)` compressions per root-to-leaf path, the
primitive count is constant.

### The constructor, written (`symbolic.py`)

The last step is now implemented: `build_state` emits the primitives as **data**
— explicit unary `Step` objects for `α, β, h↑, h↓` plus cap bounds — rather than
evaluating the normal form on the fly. Queries are answered from that list
alone, and agree with exact 4-D enumeration in 216 of 216 checks across both
measures.

Boundaries are obtained by evaluating each defining predicate at the finitely
many candidate breakpoints, which is exact (everything in sight is piecewise
constant with breakpoints in that set) and avoids hand-deriving compositions of
generalized inverses. The constructor then *asserts* the orientation it expects
of each boundary. That assertion earned its place immediately: the first run
caught three inverted conditions — a predicate passed as its own complement, two
swapped selector negations, and a reversed branch polarity — which would
otherwise have surfaced only as wrong numbers.

The two structural claims measured on the emitted data:

| `m` | primitives (merged) | records / `m` |
|---|---|---|
| 12 | 44 | 66 |
| 18 | 117 | 108 |
| 24 | 113 | 160 |
| 36 | 132 | 172 |
| 48 | 132 | 92 |

The primitive count is bounded by the signed expansion of the 18 region labels
(196 before merging) and **does not scale with the instance** — it saturates
around 130, the same saturation signature as the region count. Records are
`O(primitives × m)`, hence `O(m)` for a constant primitive count. Before
merging the count is exactly 196 at every instance size, which is the cleanest
possible statement of the claim.

One honest qualification: the constant is large. "A dimension-dependent constant
number of primitives" is literally true, but in four dimensions that constant is
of order 10² before anyone tries to optimise it, and the report's complexity
statement would read more usefully with that number in it.

### What is still not proved

Nothing in the normal form itself. What remains is integration: wiring the
compressor into Chan's recursion behind `compression_backend="prefix4"` and
measuring end-to-end, which is Phase 9 of the task specification. That is an
engineering exercise against a now-validated component, not an open question.

## What is not verified

**The parametric closure claim (Lemma `lem:parametricclosure`) is not tested
here, and it is the one I would attack next.** *(Qualified by the section above:
the source-side composition is now tested and closes; the remaining untested
piece is the final re-expression into one latent pair.)*

What I verified splits into two halves:

* the **static** normal form — a six-staircase state has a constant-size
  four-prefix representation (Thm. `thm:prefix`): tested, holds;
* the **semantic** closure — a compressed state may replace the fine geometry
  for all descendants (eq. 457): tested, holds, with the hypothesis shown to be
  necessary.

The untested half is the **representation-level** closure: that integrating a
term which *already contains* a `J`-primitive against a fresh family of easy
masks again yields `O(1)` primitives. This is what makes the scheme recursive
rather than single-shot, and it is the step the report itself identifies as
having previously been implicit ("the simplification claim no longer relies on
an unstated persistence or closure assumption").

The report is candid that one compression may multiply the primitive count by a
constant `C₄ > 1`, and argues that a root-to-leaf path sees only `O(1)`
compressions, so `C₄^{O(1)}` stays constant. That argument is valid *given* the
per-compression constant. What nobody has measured is `C₄` itself. If `C₄`
turned out to depend on the incoming primitive count rather than on the
dimension alone, the state would grow geometrically in the number of
compression generations and the complexity claim would fail while every test in
this repository still passed — because all of my tests compress exactly once.

Concretely, the missing experiment is Phase 8 of the task specification run to
2–4 generations *through the symbolic representation* rather than through the
cell-average function: build the `J`-primitives, absorb, rebuild primitives,
and count them at each generation. That requires a symbolic `PrefixState` with
latent variables, which this prototype does not have — it evaluates the normal
form but does not carry it as data.

Two smaller things are also unverified: no claim about the end-to-end
`Õ(n^{4/3})` bound is tested, since the compressor is not wired into Chan's
recursion; and the magnitude identity `Mag(D) = μ_M(D)` is taken from the cited
companion paper. I verified its *consequence* — that the product measure
expands into `Σ_S 2^-|S| HV(π_S D)` for downward-closed `D`, which needs only
that the slice at 0 equals the projection — but not that this quantity is the
magnitude of the metric space.

---

## Assessment

The construction is in better shape than I expected going in. The `2+2`
factorisation does what the report claims: the target plane splits into a
constant number of fixed regions, the source mass is separable on each, and the
prefix arguments act only as caps. That is a genuine simplification over a
generic symbolic elimination engine, and it is the part most likely to be
doubted, so it is good that it is also the part that tested cleanest.

The exponent claim is explicitly *not* an improvement over Chan, and the report
says so plainly, which is the right posture.

The remaining risk is concentrated in one place — whether the primitive count
stays bounded across compression generations — and it is testable. Until it is
tested, the honest status is: the normal form is verified, the semantics of
compression are verified, and the recursion-level complexity argument rests on a
constant that has not yet been measured.

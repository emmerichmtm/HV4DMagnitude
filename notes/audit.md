# Audit of the product-measure / magnitude simplification of Chan's 4-D algorithm

Report audited: *A Product-Measure and Magnitude Simplification of Chan's
Four-Dimensional Hypervolume Algorithm*.

**Status update.** A revised version of the report
(`paper/magnitude_chan4_report_closed.tex`) now closes the open point and
incorporates every finding below. The audit text is kept as written, with the
resolution recorded at the end, so the trail from finding to fix stays
readable.

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

### The compressor A/B (Phase 7)

The prefix compressor was compared directly against FastHVChan's
`compress_terms` on identical grounded 4-D states. Both claim the same defining
property -- the compressed state integrates to the original mass over every cell
of the hard grid -- so both were checked against exact enumeration on every
cell:

```
32 hard-grid cells: legacy and prefix compressors both agree with exact mass
```

A cell mass from the prefix state is the 16-corner finite difference of the
four-prefix function, arranged as four states (the source caps `a, b` select
which regions are rectangles, so they are fixed at build time) each queried at
four `(c, d)` corners. That is the report's "16 prefix evaluations per grid
cell", realised.

So the prefix compressor is a valid replacement for `compress_terms` on the
grounded four-dimensional path.

### What is still not done, and one reason it matters less than expected

Full recursion integration (Phase 9) -- running the solver with
`compression_backend="prefix4"` end to end -- is not implemented, for a concrete
reason: it needs a constructor that `symbolic.py` does not provide. Compression
at depth must build a state from the newly easy six staircases **times an
incoming J-primitive**, whereas `build_state` handles the generation-0 case
(six staircases alone). The composition is known to close -- that is what
`closure_probe.py` measures -- but emitting its primitives as data is a further
derivation, not a wiring exercise.

There is also an empirical reason the end-to-end A/B would say less than one
would hope. On grounded 4-D workloads the legacy compressor **never fires**:

```
n=20   nodes=55    term_high_water=1  compressions=0
n=40   nodes=155   term_high_water=1  compressions=0
n=80   nodes=419   term_high_water=1  compressions=0
n=160  nodes=1077  term_high_water=1  compressions=0
```

For grounded orthants, absorption only adds or merges conditions and never
multiplies terms, so the representation never outgrows its post-compression
size and the adaptive trigger never pulls. Both the component being replaced and
its replacement are therefore dead code on exactly the workload the report
targets, unless compression is forced on. An end-to-end comparison would be
measuring two things that do not run.

That is worth knowing independently of this prototype: it says the compression
machinery -- the hardest part of Chan's Section 4 and the part this report sets
out to simplify -- is not what costs time on the grounded four-dimensional
hypervolume problem. The cost is in the cutting recursion and in the per-node
integration. A simplification of compression is still worth having for the
general case and for the proof, but it should not be expected to show up in
grounded 4-D hypervolume timings.

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


---

## Resolution in the revised report

The revised report adopts all four findings into the formal statement rather
than into footnotes, and supplies a proof of the step that was open.

| Finding | Where it now lives |
|---|---|
| 1. Strict residual inequality | the residual is now `∏ [x_j < f_ij(x_i)]`, and the staircase prefix formula is stated for the strict residual |
| 2. `sup ∅` must not be the number 0 | the endpoint convention now mandates *decorated interval bounds* — EMPTY/FULL sentinels plus an open/closed flag — and says explicitly that an empty prefix is never encoded by the endpoint 0 |
| 3. Measurability is load-bearing | formalised as a **future-complete** grid: every threshold that can become a descendant mask, including cut coordinates and inherited cell boundaries, must be grid-measurable |
| 4. Monotone images are new values | stated directly — the images "can be new numerical values" — with the `O(m)` bound re-derived on that basis |

**The open point is closed.** The parametric re-expression lemma now carries an
explicit constructor: after an incoming primitive is multiplied by newly easy
masks and integrated, every residual comparison is routed by variable type into
one of four fixed destinations — a physical staircase, a latent `u` cap, a
latent `v` cap, or a one-turn latent boundary — with a crude but dimension-only
bound `C₄ ≤ 16 R₀ K₀⁸ 2^{L₀}`.

Two things about that proof are worth recording from an auditor's seat.

First, it is careful in exactly the place my measurements were ambiguous. It
separates *primitive count* from *breakpoint count*: a same-variable comparison
may alternate arbitrarily often and is stored as one unary step filter rather
than split into primitives. That is precisely why the primitive count cannot
depend on staircase complexity — and it matches what I measured, where the raw
primitive count stayed flat at 196 while breakpoint complexity grew with the
instance.

Second, the complexity claim is now weaker and more defensible than the earlier
draft's. It no longer asserts that the primitive count returns to a constant
after each compression; it asserts only `S' ≤ C₄ S`, and then handles the
accumulation by schedule — `O(1)` generations under a fixed block parameter, or
`O(log log n)` generations and hence `(log n)^{O(1)}` primitives under a locally
recomputed `r_j = N_j^δ`. Both are absorbed by `Õ`, and the report says outright
that no canonical recombination theorem is claimed. That is the right call:
those two situations are not distinguishable by the experiments here.

**What my experiments do and do not say about the new proof.** They are
consistent with it, and are cited as such, but they are not a check of it. I
measured the composition's region labels and separability — the mechanism the
constructor formalises — and multi-generation compression through cell
averages. I did not implement the constructor, so `C₄` itself remains
unmeasured: the 196 figure is the generation-0 count from `symbolic.py`, not
`S₁/S₀`.

**What this unblocks.** The appendix gives the recipe for the composed
constructor whose absence blocked Phase 9, as a normalization table from
comparison type to destination. The remaining engineering is now specified
rather than open. The empirical caveat from the previous section still stands:
compression does not fire on grounded 4-D workloads, so an end-to-end A/B would
measure a path that does not execute unless compression is forced on.


---

## End-to-end: the normal form computes hypervolume inside Chan's recursion

`solver.py` is a grounded 4-D hypervolume solver that keeps Chan's structure
unchanged -- absorption, weighted-median cuts on `(d-3)`-faces, inclusion-
exclusion at the base case -- and replaces only *how the easy residual is
integrated*, routing it through the four-prefix normal form. A box integral is
the signed 16-corner difference of prefix queries.

It agrees with the independent `chan_hypervolume` implementation on every
instance tried: random fronts up to `n=20`, fully non-dominated spherical fronts
up to `n=40`, tie-heavy integer grids, duplicates, nested points, and a point at
the origin.

More tellingly, the recursion tree matches:

| `n` | nodes, prefix solver | nodes, reference d/3 solver |
|---|---|---|
| 5 | 5 | 5 |
| 10 | 19 | 19 |
| 20 | 61 | 61 |
| 40 | 153 | 159 |

Identical up to `n=20`; the small divergence at 40 is tie-breaking in the
weighted median, the same effect already documented between the C and Python
ports of the reference. So the cut decisions really are Chan's, and the
substitution is confined to the integration step, as the report intends.

**No compression is needed for this to work.** Compression bounds symbolic
growth, and on grounded orthants it never triggers; the easy state at any cell
is always exactly six staircases, rebuilt from the orthants active there. That
is why a working solver was reachable without the parametric re-expression
constructor, and it is a second, independent reason the end-to-end A/B
discussed earlier would have measured very little.

### Two bugs found by doing this, both worth recording

**1. The easy orthants must be passed down through a cut.** My first version
recursed with only the hard orthants, on the reasoning that the easy ones had
been absorbed into the staircases. They had -- but into *this* cell's
staircases. Below a cut each child rebuilds its own, so an orthant absorbed at
the parent still covers volume in the children and must be handed down. (Its
active-constraint count can only fall as the cell shrinks, so an easy orthant
stays easy; nothing is lost by re-absorbing.) Results were correct whenever the
base case was reached without cutting, which is why small `n` passed.

**2. Staircases must be clamped to the cell, or the separable split diverges.**
This one is a genuine remark about the normal form, not about my code. Branch 2
of the staircase formula is split as `P(A) + Q(B)` precisely to make the two
factors separable. But

    P(A) = ∫ Y(f(x)) dμ₁

diverges when `f ≡ +∞`, which is exactly what a pair with no quadrant looks
like. The split then evaluates as `∞ − ∞` and returns NaN, while the quantity
it represents is perfectly finite. Inside a bounded cell the repair is exact
rather than a patch: a vacuous constraint *is* the cell's ceiling on that axis,
so clamping `f_ij ← min(f_ij, hi_j)` restores both finiteness and separability.

This is worth a sentence in the report. The separable split is what the whole
normal form is built on, and it is only valid for staircases that are finite on
the cell. Stated as "clamp every staircase to the cell before splitting" it
costs nothing and removes a failure mode that is invisible until a pair happens
to carry no quadrant.


---

## Was there a speed benefit? Was magnitude used at all?

Three things get conflated under "the magnitude simplification". Separating
them, and measuring:

### 1. The magnitude *measure* was never used to compute a hypervolume

`solver.py` runs on `[LEBESGUE] * 4`. The hypervolume **is** the Lebesgue
specialisation of the product measure, so the atom-at-zero measure plays no part
in the computation. Nor would the recovery identity
`HV₄ = 16·Mag − 16 − 8H₁ − 4H₂ − 2H₃` help: it needs `Mag`, which costs what
`HV₄` costs, *plus* fourteen lower-dimensional projections. It is a validation
identity, not a shortcut.

Where magnitude did earn its keep is methodological, and it was not small. The
atom at the anchor makes endpoint errors *visible* that Lebesgue silently
absorbs. Five separate bugs in this work were caught that way and would
otherwise have shipped: the branch-1 short-circuit on the anchor prefix, an
exclusive bound consumed as closed in `cap_u`, a merge key omitting the caps,
and two inverted region conditions. Every one of them was invisible under
Lebesgue. As a test oracle, magnitude is worth keeping even if it never appears
in a production path.

### 2. The product-measure normal form *was* used

The `2+2` prefix machinery does all of the easy-residual integration in
`solver.py`. That part of the proposal is genuinely exercised.

### 3. No speed benefit — and none is claimed

Spherical fronts, 4-D, seconds:

| `n` | prefix (this work) | Chan d/3 | Chan d/2 | dimension sweep | WFG |
|---|---|---|---|---|---|
| 10 | 0.196 | 0.014 | 0.0028 | 0.0005 | 0.0010 |
| 20 | 0.527 | 0.059 | 0.0143 | 0.0030 | 0.0123 |
| 40 | 3.338 | 0.284 | 0.0292 | 0.0058 | 0.0218 |
| 80 | 4.569 | 0.246 | 0.0518 | 0.0138 | 0.0547 |

At `n = 40` the prefix solver is ~240x slower than the dimension sweep, ~12x
slower than the reference Section-4.2 solver and ~100x slower than Section 2.

**How much of that is the prototype rather than the method?** Most of it, and
the profile says where: 99% of the time is in `K_decomposed`, which answers each
prefix query by sweeping the target grid afresh — about 65 target cells per
query, 84,076 integrand evaluations for 1,296 queries. The cost per query also
*grows* with the instance (0.37 → 0.42 → 0.91 ms from `n`=10 to 40), so the
query is not `O(1)`; the report's design says it should be, answered from
precomputed records with a predecessor search at worst.

So the fair reading is: this measurement does **not** show the method is slow,
only that no speed benefit has been demonstrated. A faithful implementation
would build records once per cell and answer queries in `O(1)`, putting
per-node cost in the same class as the reference Section-4.2 solver — which is
itself ~50x slower than the dimension sweep on this workload.

That ordering is the thing to keep in view. The report is explicit that it does
not improve the exponent and offers a simpler simplification step, not a faster
algorithm. Nothing here contradicts that, and nothing here supports a speed
claim either. On grounded 4-D hypervolume the practical ranking remains what the
FastHVChan benchmarks already showed: dimension sweep first, then WFG and
Chan d/2, with the d/3 family well behind on constants despite the better
exponent.

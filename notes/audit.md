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

## What is not verified

**The parametric closure claim (Lemma `lem:parametricclosure`) is not tested
here, and it is the one I would attack next.**

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

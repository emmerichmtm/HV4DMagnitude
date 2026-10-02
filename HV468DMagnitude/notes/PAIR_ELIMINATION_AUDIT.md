# Audit: the general pair-elimination simplification

Report audited: *A Product-Measure and Magnitude Simplification of Chan's
Hypervolume Algorithm* — generalized version with Section "General
fixed-dimensional pair elimination" (`paper/magnitude_chan_pair_elimination_report.tex`).

**Status: mathematically plausible but implementation-incomplete.**

That is the prompt's middle category, and it is chosen deliberately. No
counterexample was found. Everything the constructor covers is exact, under both
measures. But the constructor's coverage falls off sharply with dimension, and
at `p = 6` it covered nothing, so the general-dimensional claim — the entire
point of the generalization — is **not yet tested**. The 4-D base case remains
verified from the previous sprint.

---

## What was tested, and what it says

Elimination was checked against exact product-grid enumeration (anchor atom as
its own cell, so magnitude is exact, not approximated):

| `p` | instances completed | max branching per elimination | exactness |
|---|---|---|---|
| 2 | 30/30 | 2 | exact |
| 3 | 18/30 | 5 | exact |
| 4 | 11/30 | 9 | exact |
| 5 | 1/6 | 2 | exact |
| 6 | **0/6** | — | **untested** |

"Branching" is the number of terms one elimination produces from one incoming
term — the quantity the lemma actually bounds. A cumulative count would
conflate it with chain depth, which is how the first version of this
instrumentation misled me.

Where the constructor completed, the result matched the oracle every time, for
Lebesgue and for magnitude. That is real but narrow evidence.

---

## Finding 1 — the term class is not closed, and the abstentions prove it

Definition 1 requires every `f_e` to be a **monotone** step map, and justifies
"at most one lower and one upper predicate per pair" by saying same-orientation
predicates merge via min/max.

That merge is only available when the orientations agree. Elimination routinely
produces, for one pair, an upper predicate from a nondecreasing map *and* one
from a nonincreasing map. Their conjunction is

```
xi_j  <|  min( h_up(xi_i), h_down(xi_i) )
```

which has a turning point and is **not monotone**. The class as defined does not
contain it.

This is not speculation about the proof. The implementation abstains exactly
when it meets such a predicate, and the abstention rate is the table above:
0% at `p = 2`, 40% at `p = 3`, 63% at `p = 4`, 83% at `p = 5`, 100% at `p = 6`.
The failure mode grows with dimension precisely because more eliminations mean
more chances for two opposite orientations to land on the same pair.

**This is a repair, not a refutation.** The 4-D treatment already contains the
fix: its "latent one-turn normal form" is exactly `min(h_up, h_down)`. Two ways
to close the general definition:

1. allow one-turn boundaries in the class, as 4-D does; or
2. keep a *list* of predicates per pair rather than one per side.

Either keeps the count bounded for fixed `p` — each elimination adds only
`O_p(1)` predicates — so no asymptotic claim changes. But Definition 1 as
written is too tight, and the proof's step "each labelled region is itself a
weighted monotone-pair term" does not follow from it.

The implementation takes route 2 (`pair_state.merge_same_orientation` merges
only within an orientation and keeps the rest). Route 1 is what is needed to
close the remaining abstentions, and is the next piece of work.

## Finding 2 — the separable split needs the staircase clamped

Carried over from the 4-D audit and unchanged by the generalization. Branch 2 of
the staircase prefix formula is split as `P(A) + Q(B)` precisely to make the two
factors separable, but

```
P(A) = ∫ Y(f(x)) dμ
```

diverges when `f ≡ +∞`, which is what a pair carrying no constraint looks like.
The split then evaluates as `∞ − ∞` while the quantity it represents is finite.

In the general setting the repair is available and exact: in the compression
integrand every `x_i` is confined to a grid interval, so a vacuous constraint
*is* the cell ceiling and `f_e` can be clamped to it. This must be stated as
part of the normalization, not left to the implementer — the 4-D solver
returned `NaN` on real instances before it was added.

## Finding 3 — the decorated lower bound inverts in the cumulative difference

Found by implementing, and the single genuine bug of this sprint. The mass of a
feasible interval is

```
Q(U) − Q(L⁻)
```

so on the **lower** side the endpoint flag must be *inverted*: if `L` is
attained, what gets subtracted is the mass strictly below it; if `L` is
excluded, the mass up to and including it. Using the flag directly deletes the
anchor atom whenever `L = 0`.

Under Lebesgue this is invisible — every test passed. Under magnitude it is a
mass-1 error, and it broke every magnitude instance from `p = 2` upward. This is
now the fourth distinct endpoint bug across the two sprints that Lebesgue hid
and magnitude exposed, which is the strongest practical argument for keeping the
magnitude oracle in the test suite.

## Observation — Lemma 3 follows from Lemma 2

The report proves one-variable elimination (Lemma 2) and pair elimination
(Lemma 3) separately. For *closure* purposes Lemma 3 is not needed: from `x`'s
point of view the mutual `x`–`y` staircase is just one more bound candidate that
happens to depend on `y`, so applying Lemma 2 twice already gives the pair
result. The implementation does exactly this and is exact wherever it applies.

Lemma 3 buys a better constant and a closed form (the two-branch formula), not
extra generality. Worth saying in the paper, because it makes the general
argument rest on the simpler of the two lemmas.

---

## The prompt's seven obligations

| # | obligation | status |
|---|---|---|
| 1 | finite active-bound candidate count | holds: domain bound plus one per incident predicate, `O(p)` per side |
| 2 | finite label count | holds: product of the two candidate sets, `O(p²)`; measured max 9 at `p=4` |
| 3 | separable branch formulas | verified wherever the constructor applies, both measures |
| 4 | closure of branch conditions | **fails as stated** — see Finding 1 |
| 5 | complement handling | not reached; blocked behind Finding 1 |
| 6 | breakpoint complexity vs primitive count | partially: branching stayed at 2, 2, 5 while term breakpoints went 2, 6, 9, but the sample is small and biased toward completed cases |
| 7 | no hidden instance-size dependence in `B_p` | consistent with the data, not established — the coverage is too thin to claim it |

---

## What would close this

1. Extend the predicate class to one-turn boundaries (`min(h_up, h_down)`),
   mirroring the 4-D normal form. This should remove most abstentions and is
   the prerequisite for everything below.
2. Re-run `p = 6` and `p = 8` with the extended class, on tiny grids, against
   exact enumeration.
3. Only then measure `B_p` against growing `N` and claim independence.
4. Phase 5 (compression in `2d` variables) and Phase 6 (FastHVChan backend
   behind a flag) have not been started.

Until step 2 produces exact agreement at `p = 6`, the honest statement is that
the generalization is plausible and partially implemented, and that its 4-D
specialization is the only dimension verified end to end.

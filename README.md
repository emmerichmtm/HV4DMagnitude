# HV4DMagnitude

An adversarial validation of the product-measure / magnitude simplification of
Chan's four-dimensional hypervolume algorithm, following the task specification
in `CODEX_CLAUDE_MAGNITUDE_CHAN_PROMPT_REVISED.md`.

The purpose of this folder is **to try to break the proposed lemmas**, not to
ship a fast solver. Every identity is implemented independently of an exact
brute-force oracle and checked against it, including the degenerate cases the
report flags as dangerous: plateaus, queries exactly at breakpoints, empty
feasible intervals, and the magnitude atom at the anchor 0.

The report under audit is in [`paper/`](paper/). Its revised version proves the
step that was open and adopts all four findings into the formal statement.

**Read [`notes/audit.md`](notes/audit.md) first** — it states what was verified,
what needs correcting in the report, and, most importantly, what is *not*
verified.

## Result in one paragraph

No counterexample to any stated lemma survived. The core claim — that the target
plane splits into a constant number of fixed regions on each of which the source
mass is separable — held across 11,250 pointwise checks, and the assembled
four-prefix integral matched exact 4-D enumeration under both measures. Three
items need sharpening in the write-up (a strict-inequality convention, the
`sup ∅` encoding, and the measurability hypothesis in the tower argument), and
one substantial claim was untested: whether the primitive count stays bounded
across *repeated* compressions. That has now been tested too
(`generations.py`): compression is exact over four generations against the
original fine geometry, and the per-generation blow-up factor is stable at
17-21 regions against a dimension-only bound of 50. What remains is to write
the symbolic constructor for the output primitives — engineering, not
discovery.

## Running

```bash
python run_all_tests.py
```

No dependencies beyond the standard library. The FastHVChan cross-check imports
`C:\MyTemp\code\FastHVChan` and is skipped if that is not present.

## Layout

```
core.py          measures, right-continuous step functions, generalized
                 inverses, and the exact grid oracle
identities.py    the three identities: threshold antiderivative, two-branch
                 staircase prefix, meet contraction
prefix4.py       the six-staircase model and the four-prefix normal form,
                 with region labelling and instrumentation
symbolic.py      the normal form emitted as data: explicit primitives with
                 unary Step records, answering queries by arithmetic
closure_probe.py measures whether the blow-up factor grows with the instance
generations.py   multi-generation compress/absorb, verified against fine geometry
tests/           one file per phase group of the validation plan
notes/audit.md   the audit
paper/           the report under audit (revised; the open lemma is proven)
```

## Conventions

These are stated in full at the top of `core.py` and are the whole ballgame for
the magnitude measure; the short form:

* coordinates live on `[0, ∞)` and are never translated — the anchor 0 carries
  an atom under magnitude, so a shifted cell is a different problem;
* measures are `atom0·δ₀ + density·λ`; Lebesgue is `(0, 1)`, magnitude is
  `(1, ½)`;
* step functions are **right-continuous**, which the geometry forces: the
  staircase bounding a union of grounded quadrants drops exactly at a vertex;
* staircase constraints are **strict** (`x_j < f(x_i)`), which is what the
  geometry gives and which differs from the report's `≤` under magnitude;
* generalized inverses return **exclusive** bounds, so `θ = 0` means the empty
  set and not `{0}`.

## Status against the task's acceptance criteria

| Criterion | Status |
|---|---|
| Threshold-antiderivative tests pass, both measures, ties and anchor | ✅ 720 cases |
| Four-prefix queries agree with small-grid enumeration | ✅ 144 queries |
| Compression agrees after absorb + coarsen | ✅ 192 cells, 1 generation |
| Closure over *multiple* generations | ✅ exact over 4 generations (456 cells, error 0 / 4e-15); `R` stable at 17-21 |
| Compressed state is `O(m)` with no Cartesian blowup | ✅ primitives constant in `m` (196 raw, ~130 merged); records `O(m)` |
| Normal form emitted as data | ✅ `symbolic.py`; 216/216 queries from primitives alone |
| New compression agrees with legacy `compress_terms` | ✅ 32 hard-grid cells, both against exact mass |
| End-to-end `prefix4` agrees with the legacy backend | ❌ needs the composed constructor; and compression never fires on grounded 4-D anyway — see the audit |
| Cutting recursion unchanged | ✅ untouched — nothing in `FastHVChan` was modified |
| Documented `Õ(n^{4/3})` overall bound | ❌ not tested here |

`FastHVChan` is used read-only as an independent oracle. Nothing in it was
changed.

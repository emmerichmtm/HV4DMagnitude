# HV468DMagnitude

Verification of the **general pair-elimination** simplification of Chan's
hypervolume algorithm, for fixed dimensions 4, 6 and 8.

The 4-D case is verified end to end in the parent project. This one tests the
generalization: that eliminating variables from a *weighted monotone-pair term*
stays inside the class, with branching bounded by dimension alone.

**Status: mathematically plausible but implementation-incomplete.**
No counterexample. Everything the constructor covers is exact under both
measures. But coverage falls with dimension and reaches 0% at `p = 6`, so the
general-dimensional claim is not yet tested.

Read [`notes/PAIR_ELIMINATION_AUDIT.md`](notes/PAIR_ELIMINATION_AUDIT.md).

## Results so far

| `p` | completed | max branching / elimination | exactness |
|---|---|---|---|
| 2 | 30/30 | 2 | exact |
| 3 | 18/30 | 5 | exact |
| 4 | 11/30 | 9 | exact |
| 5 | 1/6 | 2 | exact |
| 6 | 0/6 | — | untested |

The abstentions are not random: the constructor declines exactly when a pair
acquires two predicates of opposite orientation, whose conjunction is one-turn
and leaves the monotone class. That is Finding 1 of the audit.

## Running

```bash
python tests/test_elimination.py
```

## Layout

```
pair_state.py        weighted monotone-pair terms, decorated bounds
pair_elimination.py  the elimination constructor plus instrumentation
pair_oracles.py      exact product-grid enumeration (atom as its own cell)
tests/               elimination against the oracles
notes/               the audit, and the task prompt
```

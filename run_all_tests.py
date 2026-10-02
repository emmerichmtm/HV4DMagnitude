#!/usr/bin/env python3
"""Run every phase of the validation plan and summarise.

    python run_all_tests.py
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUITES = [
    ("Phases 1-4  identities", "tests/test_identities.py"),
    ("Phases 5-6  four-prefix normal form", "tests/test_prefix4.py"),
    ("Phase 6e-f  normal form as data", "tests/test_symbolic.py"),
    ("Phases 8,10 closure and magnitude oracle", "tests/test_closure.py"),
    ("Phase 7     vs legacy compress_terms", "tests/test_vs_legacy_compressor.py"),
    ("Phase 9     cross-check vs FastHVChan", "tests/test_against_fasthvchan.py"),
    ("Phase 9     end-to-end 4-D hypervolume solver", "tests/test_solver.py"),
    ("Phase 8+    multi-generation compression", "generations.py"),
    ("Phase 11    blow-up factor vs instance size", "closure_probe.py"),
]


def main() -> int:
    failed = []
    for title, path in SUITES:
        print("=" * 72)
        print(title)
        print("=" * 72)
        done = subprocess.run([sys.executable, os.path.join(HERE, path)],
                              cwd=HERE)
        if done.returncode != 0:
            failed.append(title)
        print()
    print("=" * 72)
    if failed:
        print("FAILED: " + "; ".join(failed))
        return 1
    print("All phases passed.")
    print("See notes/audit.md for what this does and does not establish.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

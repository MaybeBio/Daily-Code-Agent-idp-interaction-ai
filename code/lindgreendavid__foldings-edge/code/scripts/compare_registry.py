#!/usr/bin/env python3
"""Compare research registries without requiring cross-platform float byte identity."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

REL_TOLERANCE = 1e-9
ABS_TOLERANCE = 1e-12
CLASSIFIER_BOOTSTRAP_CI_ABS_TOLERANCE = 1e-5


def _absolute_tolerance(path: str) -> float:
    """Allow only seeded classifier bootstrap boundaries to vary by one draw."""
    is_classifier_bootstrap_ci = path.startswith("$.h2_classifier.") and (
        ".f1.ci_95_" in path or ".mcc.ci_95_" in path
    )
    return CLASSIFIER_BOOTSTRAP_CI_ABS_TOLERANCE if is_classifier_bootstrap_ci else ABS_TOLERANCE


def assert_registry_close(expected: Any, actual: Any, path: str = "$") -> None:
    """Require identical structure and tightly equivalent numeric values."""
    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected is not actual:
            raise AssertionError(f"{path}: expected {expected!r}, got {actual!r}")
        return

    if isinstance(expected, dict) and isinstance(actual, dict):
        if expected.keys() != actual.keys():
            missing = sorted(expected.keys() - actual.keys())
            extra = sorted(actual.keys() - expected.keys())
            raise AssertionError(f"{path}: key mismatch; missing={missing}, extra={extra}")
        for key in expected:
            assert_registry_close(expected[key], actual[key], f"{path}.{key}")
        return

    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            raise AssertionError(f"{path}: expected {len(expected)} items, got {len(actual)}")
        for index, (expected_item, actual_item) in enumerate(zip(expected, actual, strict=True)):
            assert_registry_close(
                expected_item,
                actual_item,
                f"{path}[{index}]",
            )
        return

    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        expected_float = float(expected)
        actual_float = float(actual)
        absolute_tolerance = _absolute_tolerance(path)
        if not math.isfinite(expected_float) or not math.isfinite(actual_float):
            raise AssertionError(f"{path}: registry numbers must be finite")
        if not math.isclose(
            expected_float,
            actual_float,
            rel_tol=REL_TOLERANCE,
            abs_tol=absolute_tolerance,
        ):
            raise AssertionError(
                f"{path}: expected {expected!r}, got {actual!r} "
                f"(rel_tol={REL_TOLERANCE}, abs_tol={absolute_tolerance})"
            )
        return

    if type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"{path}: expected {expected!r}, got {actual!r}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expected", type=Path)
    parser.add_argument("actual", type=Path)
    args = parser.parse_args()

    expected = json.loads(args.expected.read_text())
    actual = json.loads(args.actual.read_text())
    assert_registry_close(expected, actual)
    print(
        "registries match structurally and numerically "
        f"(rel_tol={REL_TOLERANCE}, abs_tol={ABS_TOLERANCE})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

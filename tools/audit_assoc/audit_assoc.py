#!/usr/bin/env python3
"""audit_assoc.py -- exact tie-aware rank association for small panels.

This is the dependency-free reference implementation of the rank-inference
core used in the PARP1 GaMD audit.  It takes precomputed descriptor values
(one row per system), a phenotype rank vector, and returns, for every
descriptor column:

  * Spearman rho with midranks (average ranks) for ties;
  * the two-sided exact permutation p-value over all n! orderings;
  * the number of permutations as or more extreme than the observation;
  * the smallest attainable p (p_min) and the smallest |rho| that reaches
    p < alpha for this design.

It does not reweight trajectories, reconstruct PMFs, or choose a production
window: those steps remain study-specific.  The CSV schema is documented in
README.md and illustrated by example/panel_example.csv.

Usage
-----
    python3 tools/audit_assoc/audit_assoc.py \
        --input tools/audit_assoc/example/panel_example.csv \
        --metrics span_minmax,span_q90,span_t95,span_support \
        --out /tmp/audit_assoc.json

    python3 -m unittest discover -s tools/audit_assoc/tests -v
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

MAX_EXACT_N = 8  # 8! = 40,320 permutations; above this, use a Monte Carlo test.


def midrank(values: Sequence[float]) -> List[float]:
    """Average ranks, with tied values sharing the mean of their positions."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j + 2) / 2.0  # 1-based positions i+1 .. j+1
        for k in range(i, j + 1):
            ranks[order[k]] = mean_rank
        i = j + 1
    return ranks


def spearman_midrank(x: Sequence[float], y: Sequence[float]) -> float:
    """Tie-aware Spearman rho; returns 0.0 when either margin is constant."""
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("x and y must have equal length >= 2")
    rx, ry = midrank(x), midrank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0.0 or syy == 0.0:
        return 0.0
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry)) / math.sqrt(sxx * syy)


def exact_test(x: Sequence[float], y: Sequence[float], alpha: float = 0.05) -> Dict[str, object]:
    """Two-sided exact permutation test of |rho| under all n! relabelings."""
    n = len(x)
    if n < 2:
        raise ValueError("need at least two systems")
    if n > MAX_EXACT_N:
        raise ValueError(
            f"n={n} exceeds the exact-enumeration limit ({MAX_EXACT_N}); "
            "this reference implementation intentionally does not approximate."
        )
    observed = spearman_midrank(x, y)
    tol = 1e-12

    # Full permutation distribution of |rho|; the same distribution yields
    # p_min and the smallest |rho| that can reach alpha for this design.
    abs_rhos: List[float] = []
    for perm in itertools.permutations(y):
        abs_rhos.append(abs(spearman_midrank(x, perm)))
    n_perm = len(abs_rhos)
    count = sum(1 for r in abs_rhos if r >= abs(observed) - tol)
    p_value = count / n_perm

    thresholds = sorted({round(r, 12) for r in abs_rhos})
    detectable = None
    p_min = 1.0
    for threshold in thresholds:
        c = sum(1 for r in abs_rhos if r >= threshold - tol)
        p = c / n_perm
        p_min = min(p_min, p)
        if detectable is None and p < alpha:
            detectable = threshold
    return {
        "n": n,
        "rho": round(observed, 6),
        "p_exact": round(p_value, 6),
        "n_permutations": n_perm,
        "n_as_or_more_extreme": count,
        "p_min": round(p_min, 6),
        "detectable_abs_rho_alpha": detectable,
        "alpha": alpha,
    }


def read_panel(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header row")
        rows = [dict(row) for row in reader]
    required = {"system", "rank"}
    missing = required.difference(reader.fieldnames)
    if missing:
        raise ValueError(f"{path} is missing required column(s): {sorted(missing)}")
    if len(rows) < 3:
        raise ValueError("panel needs at least three systems")
    return list(reader.fieldnames), rows


def parse_metrics(fieldnames: Iterable[str], requested: str | None) -> List[str]:
    candidates = [c for c in fieldnames if c not in {"system", "rank"}]
    if requested is None:
        return candidates
    metrics = [m.strip() for m in requested.split(",") if m.strip()]
    unknown = [m for m in metrics if m not in candidates]
    if unknown:
        raise ValueError(f"unknown metric column(s): {unknown}; available: {candidates}")
    return metrics


def run(input_path: Path, metrics: List[str], alpha: float) -> Dict[str, object]:
    _, rows = read_panel(input_path)
    y = [float(row["rank"]) for row in rows]
    systems = [row["system"] for row in rows]
    results: Dict[str, object] = {
        "input": str(input_path),
        "systems": systems,
        "phenotype_rank": y,
        "n": len(rows),
        "alpha": alpha,
        "metrics": {},
    }
    for metric in metrics:
        try:
            x = [float(row[metric]) for row in rows]
        except (TypeError, ValueError) as exc:
            raise ValueError(f"non-numeric value in column {metric!r}") from exc
        results["metrics"][metric] = exact_test(x, y, alpha=alpha)
    return results


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Exact tie-aware Spearman audit for a small system panel."
    )
    parser.add_argument("--input", required=True, type=Path, help="panel CSV")
    parser.add_argument("--metrics", default=None,
                        help="comma-separated descriptor columns; default: all non-system/rank columns")
    parser.add_argument("--alpha", default=0.05, type=float, help="nominal level (default 0.05)")
    parser.add_argument("--out", default=None, type=Path, help="JSON output path (default: stdout)")
    args = parser.parse_args(argv)

    fieldnames, _ = read_panel(args.input)
    metrics = parse_metrics(fieldnames, args.metrics)
    if not metrics:
        parser.error("no descriptor columns found")
    payload = run(args.input, metrics, args.alpha)
    text = json.dumps(payload, indent=2)
    if args.out is None:
        print(text)
    else:
        args.out.write_text(text + "\n")
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

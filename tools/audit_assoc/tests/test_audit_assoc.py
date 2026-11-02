"""Unit tests for the reference exact rank-association CLI.

Run from the repository root:

    python3 -m unittest discover -s tools/audit_assoc/tests -v
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
TOOL_DIR = HERE.parents[1]
sys.path.insert(0, str(TOOL_DIR))

from audit_assoc import exact_test, midrank, spearman_midrank  # noqa: E402


class MidrankTests(unittest.TestCase):
    def test_ties_receive_average_rank(self):
        self.assertEqual(midrank([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])

    def test_all_tied_is_constant(self):
        self.assertEqual(midrank([3, 3, 3]), [2.0, 2.0, 2.0])

    def test_degenerate_margin_returns_zero(self):
        self.assertEqual(spearman_midrank([1, 1, 1], [1, 2, 3]), 0.0)


class ExactPermutationTests(unittest.TestCase):
    def test_monotone_panel_with_one_tie(self):
        x = [14.326, 10.375, 11.231, 8.782, 11.243]
        y = [5.0, 4.0, 2.5, 2.5, 1.0]
        result = exact_test(x, y)
        self.assertAlmostEqual(result["rho"], 0.205, places=3)
        self.assertEqual(result["n_permutations"], 120)
        self.assertEqual(result["n_as_or_more_extreme"], 92)
        self.assertAlmostEqual(result["p_exact"], 0.767, places=3)
        self.assertAlmostEqual(result["p_min"], 0.033, places=3)
        self.assertAlmostEqual(result["detectable_abs_rho_alpha"], 0.975, places=3)

    def test_no_ties_minimum_p_is_two_of_120(self):
        result = exact_test([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
        self.assertAlmostEqual(result["rho"], 1.0, places=12)
        self.assertAlmostEqual(result["p_min"], 2 / 120, places=6)
        self.assertAlmostEqual(result["p_exact"], 2 / 120, places=6)

    def test_permutation_count_is_exact_not_sampled(self):
        result = exact_test([2.0, 1.0, 4.0, 3.0, 5.0], [1, 2, 3, 4, 5])
        self.assertEqual(result["n_permutations"], 120)


class CliTests(unittest.TestCase):
    def test_cli_writes_json_for_example(self):
        example = TOOL_DIR / "example" / "panel_example.csv"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "result.json"
            proc = subprocess.run(
                [sys.executable, str(TOOL_DIR / "audit_assoc.py"),
                 "--input", str(example), "--metrics", "C2_q90", "--out", str(out)],
                capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            payload = json.loads(out.read_text())
            self.assertEqual(payload["n"], 5)
            self.assertAlmostEqual(payload["metrics"]["C2_q90"]["rho"], 0.205, places=3)
            self.assertAlmostEqual(payload["metrics"]["C2_q90"]["p_exact"], 0.767, places=3)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
#
# Purpose:  Prior-sensitivity table: P(rho<0) under flat vs arcsine (Jeffreys) priors.
# Created:  2026-09 (header standardised 2026-10-11)
# Run:      python3 scripts/s36_ordinal_prior_sensitivity.py   (from the repository root)

"""Prior-sensitivity table: P(rho<0) under flat vs arcsine (Jeffreys) priors.
Imports the d4_bayesian_ordinal model. Output: results/analysis/prior_sensitivity.txt
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
spec = importlib.util.spec_from_file_location(
    "d4b", Path(__file__).resolve().parents[0] / "d4_bayesian_ordinal.py")
d4b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d4b)

panels = {
    "nira+1": (d4b.s2cv1, d4b.sd_s2, d4b.s1, d4b.sd_s1, d4b.aai, d4b.sd_aai),
    "ruca+1": (d4b.s2cv1_ruca, d4b.sd_s2_ruca, d4b.s1_ruca, d4b.sd_s1_ruca,
               d4b.aai_ruca, d4b.sd_aai_ruca),
}

out = ["panel  metric  flat_P(rho<0)  arcsine_P(rho<0)  delta  arc_median  arc_CI95"]
for panel, (s2x, s2sd, s1x, s1sd, aaix, aaisd) in panels.items():
    for name, x_obs, sd in [("S2_CV1", s2x, s2sd), ("S1", s1x, s1sd),
                            ("AAI", aaix, aaisd)]:
        s_flat, _ = d4b.bayesian_p_rho_neg(x_obs, sd, d4b.trap_rank, prior="flat")
        s_arc, _ = d4b.bayesian_p_rho_neg(x_obs, sd, d4b.trap_rank, prior="arcsine")
        pf = (s_flat < 0).mean()
        pa = (s_arc < 0).mean()
        med = np.median(s_arc)
        lo, hi = np.percentile(s_arc, [2.5, 97.5])
        line = (f"{panel:7s} {name:6s} {pf:14.3f} {pa:18.3f} {pa-pf:+8.3f} "
                f"{med:11.3f} [{lo:.2f}, {hi:.2f}]")
        out.append(line)
        print(line, flush=True)

Path("results/analysis/prior_sensitivity.txt").write_text("\n".join(out) + "\n")
print("WROTE results/analysis/prior_sensitivity.txt")

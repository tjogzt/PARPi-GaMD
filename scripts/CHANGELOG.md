# scripts/ CHANGELOG

Changes to scripts whose outputs ship with the release are summarised here.
**RULE: any further modification of `s08_consolidate_s2_dbe.py` MUST append an
entry below** — the consolidated table is the manuscript's numerical anchor
and must stay traceable to its generator history.

## s08_consolidate_s2_dbe.py
- **v1** (2026-09-16/17): first consolidation of the S2/S1 well-depth table;
  rebuilt values from the legacy reconstruction and legacy S1 depths
  (superseded — see v2).
- **v2** (2026-10-06, P0-1 code-audit remediation): rewritten to rebuild the
  table from `data/analysis/s2_dbe_well_depths.npy` plus four curated CSVs
  (`nira_prot_dbe`, `rerun_cumulant_wells`, `s1_dbe_unified_wells`,
  `extension_s1_wells`); no hard-coded depths; column-name normalisation;
  pinned asserts against the manuscript anchor values. Verified
  value-equivalent to the published table (wells identical to 0.000000;
  max AAI delta 8.3e-5 vs the 5e-4 tolerance); the downstream guard chain
  (code/13, 18, 21) is green on the rebuilt table.

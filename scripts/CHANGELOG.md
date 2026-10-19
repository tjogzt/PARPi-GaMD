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

- **v3** (2026-10-11, release_v2 reconciliation): rewritten to the corrected-span
  data layer (r6c_s2_spanfix / r6c_s2_cv2_spanfix / r6_panel_values /
  r6c_s1_spanfix / s1_dbe_unified_wells / extension_s1_wells + the archived
  veli_prot R9 summary JSON). The v2 branch was found to recompute the
  pre-correction convention and failed its own stale anchors; v3 carries updated
  frozen-state asserts and passes a cell-for-cell golden check against the
  frozen table (33 rows identical). Changelog discipline retained.

## s55_s10_archive_frame_cvs.py / s56_s10_crossseed_spanfix.py
- **v1** (2026-10-19): added with the cross-seed stability analysis (SI Table S17).
  s55 rebuilds the per-frame records for the two additional seeds (7, 123) of
  both chemical states of the niraparib S2 system. Pairing is step-based, not
  row-based: the seed-123 neutral gamd.log lost ~51 buffered rows when the
  pre-migration process was stopped (steps 16.35-17.60 M; 2.55 ns), while the
  DCD is complete (760 frames = 38 ns). s56 recomputes the table-convention
  C3 spans (occupied-support: occupancy >= 2 bins) and the 2-ns window SDs;
  it reproduces the frozen seed-42 spanfix state exactly (span 26.926294 vs
  26.926038 kcal/mol; window SD 7.4726 vs 7.4725; occupancy bins 44/5/3) and
  asserts the frozen cross-seed summary cell-for-cell (six cells, three seeds
  per state). Outputs: data/analysis/s10/*.

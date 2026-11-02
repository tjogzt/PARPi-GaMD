# `audit_assoc` — exact rank-association audit for small panels

This directory contains the dependency-free reference implementation of the
rank-inference core used in the PARP1 GaMD audit. It answers one narrow
question exactly: given one descriptor value per system and a rank-ordered
phenotype, what is the tie-aware Spearman correlation, and how often does a
permutation of the phenotype produce a correlation at least as extreme?

It is deliberately **not** a full reimplementation of the paper's workflow.
Trajectory loading, production masking, PMF reconstruction by cumulant
expansion, quantile-span definitions, chemical-state re-simulation, and the
panel-specific phenotype vector all remain study-specific. The CLI consumes
precomputed descriptor columns.

## Input schema

A UTF-8 CSV with a header row and one row per system:

| column | required | meaning |
|---|---|---|
| `system` | yes | unique system label |
| `rank` | yes | phenotype rank (ties allowed; use midranks, e.g. `2.5`) |
| any other numeric column | at least one | descriptor value, e.g. `span_minmax`, `span_q90` |

Example (`example/panel_example.csv`, the five-system neutral panel of the
paper under the harmonized 50-ps/first-14-ns convention):

```csv
system,C2_mm,C2_q90,C3_mm,C3_q90,rank
talazoparib,19.431,14.326,30.528,18.942,5
niraparib,12.352,10.375,18.682,15.461,4
olaparib,20.533,11.231,35.973,25.262,2.5
rucaparib,13.477,8.782,31.319,22.745,2.5
veliparib,14.626,11.243,22.853,19.768,1
```

## Usage

```bash
python3 tools/audit_assoc/audit_assoc.py \
  --input tools/audit_assoc/example/panel_example.csv \
  --metrics C2_mm,C2_q90,C3_mm,C3_q90 \
  --out /tmp/audit_assoc.json
```

Output is JSON. For each descriptor it reports `rho`, the two-sided exact
`p_exact`, the number of permutations (`n_permutations`), the number at least
as extreme (`n_as_or_more_extreme`), `p_min` (the smallest attainable p for
this design), and `detectable_abs_rho_alpha` (the smallest |rho| that reaches
`--alpha`). The last two fields make the design's resolution explicit rather
than leaving a non-significant result to be read as evidence of no
association.

## Tests

```bash
python3 -m unittest discover -s tools/audit_assoc/tests -v
```

The tests check midrank handling, the known five-system cell
(`rho = +0.205`, `p = 0.767`, 92/120), the no-tie minimum `p = 2/120`, the
one-tie minimum `p = 4/120 = 0.033`, and the CLI end-to-end path.

## Limits

- Exact enumeration is capped at `n = 8` (40,320 permutations) by design.
  Larger panels need a pre-specified Monte Carlo test, which this module does
  not provide.
- The p-value conditions on the observed descriptor values; it does not model
  uncertainty in those values. Use the paper's window/seed sensitivity and
  ordinal model for that.
- Ties in the phenotype must be encoded as midranks. The module averages tied
  positions, but it does not decide which compounds are tied.

# PARPi-GaMD: Allosteric Encoding of PARP1 Trapping by Inhibitors

Gaussian accelerated molecular dynamics (GaMD) analysis of PARP1--inhibitor systems,
quantifying how structurally diverse PARP inhibitors occupying the same catalytic
pocket encode differential allosteric signals through the helical domain (HD).

## Overview

Dual-system design:
- **S1** — catalytic (CAT) domain alone (~30k atoms); CV = HD--ART COM distance
- **S2** — full-length PARP1 bound to a DNA double-strand break mimic (~286k atoms);
  CV1 = protein--DNA COM distance, CV2 = HD--ART COM distance

The **Allosteric Amplification Index (AAI)** = S1 well depth / S2 CV2 well depth
normalizes the DNA-free allosteric response against the DNA-bound baseline.

### Panel (10 systems)

| Group | Systems |
|---|---|
| Original (7) | talazoparib, olaparib, niraparib, rucaparib, veliparib, AZD5305, APO |
| Extension (3) | fluzoparib, pamiparib, senaparib |

## Protocol

- OpenMM (CUDA) + Amber ff14SB/bsc1/GAFF2/TIP3P, 300 K, 2 fs timestep
- Dual-boost GaMD (lower-dual, σ0 = 6.0 kcal/mol), 10 ns statistics-gathering phase
- Production: S1 19--31 ns (original) / 200 ns (extension); S2 24--26 ns (original) / 22 ns (extension)
- Reweighting: PyReweighting cumulant expansion (C1--C3)
  - S1: dihedral boost energy (DBE) weights
  - S2: dihedral force-weight (DFW) weights
- Ligand preparation: Vina docking → RDKit→TRIPOS mol2 (aromatic-bond corrected) → AM1-BCC charges (GAFF2)
- Zinc-finger model: non-bonded Zn²⁺ (CYM/HID coordinating residues)

## Repository contents

```
figures/pdf/  — publication figures (26, as referenced by the manuscript/SI)
figures/png/  — PNG previews of the same figures
code/         — analysis scripts (R + Python)
code/pipeline/— system building, equilibration, GaMD deployment, CV extraction
data/         — result tables (S1 PMF summary, AAI, S2 well depths, replicate data)
data/new_drugs_s1/ — extension-panel S1 CV + weights (production frames)
data/new_drugs_s2/ — extension/re-run S2 CV1/CV2 + DFW weights
```

## Reproducing

1. Environment: `mamba env create -f code/pipeline/environment_gamd_min.yml`
   (OpenMM ≥8.1 CUDA + ambertools; no rdkit required for the MD pipeline)
2. Build: `python code/pipeline/build_system2_tleap.py --receptor <pdb> --ligand <mol2> --out <dir>`
3. Equilibrate: `python code/pipeline/equilibrate.py --prmtop ... --inpcrd ... --out ... --temp 300`
4. GaMD config + run: `python code/pipeline/gen_gamd_config.py ... --prod-ns 22` then `gamdRunner xml config.xml`
5. Extract CVs: `python code/pipeline/extract_s2_cv.py runs_s2 <drug>`
6. Reweight: `python PyReweighting-1D.py -input cv.dat -T 300 -disc 0.1 -Emax 20 -cutoff 2 -job amdweight_CE -weight weights.dat`

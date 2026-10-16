# AutoDL deployment & run archive (PARP1 GaMD)

GaMD dual-system simulation toolchain for the PARP1 trapping-mechanism study.
OpenMM (CUDA) + gamd-openmm (Miao lab); open-source / commercially usable.

This directory is the archived deployment record as executed on the AutoDL GPU
instances: system builders, per-system deploy scripts, transfer helpers and
environment specifications. The standardised, re-numbered equivalents of the
generic pipeline steps live under `code/pipeline/` (p01-p09); items superseded
by those are listed in `docs/standardisation_report_v2.md`.

## 0. Environment (one-time)
Sync the project to the instance data disk, then create the environment:
    mamba env create -f environment_gamd.yml
plus PyReweighting (Miao lab; vendored under `tools/`) for PMF reweighting.

Recommended instance: RTX 4090 24G x1-2, Ubuntu 22.04 + CUDA 12.1+.

## 1. Standard pipeline (single system)
    build_system -> equilibrate -> gen_gamd_config -> gamdRunner -> analyze_gamd
(standardised copies: `code/pipeline/`; products under `md/runs/<tag>/`).

| step | script | role | product |
|---|---|---|---|
| 1 | build_system.py | pdbfixer loop completion + SystemGenerator parameterisation + solvation | .prmtop/.inpcrd/.pdb |
| 2 | equilibrate.py | min -> NVT heating -> NPT density equilibration | *_equil.rst7 |
| 3 | gen_gamd_config.py | gamdRunner configuration (lower-dual) | config_*.xml |
| 4 | gamdRunner xml config_*.xml | cMD statistics -> GaMD equilibration -> production | gamd_out/*.dcd, gamd.log |
| 5 | analyze_gamd.py | CV + boost dV -> reweighting inputs | *_cv_weights.dat |
| 6 | PyReweighting-1D.py | cumulant-2nd reweighting -> PMF | pmf-*.dat |

## 2. Reaction coordinates (CV) — PARP1 domain boundaries
- system 1 allosteric CV (tests the Zandarashvili allostery model): HD domain vs
  ART core centre-of-mass distance (HD open/closed); or the HD aF helix (~760-780)
  displacement relative to ART.
- system 2 trapping CV (tests the Hopkins catalytic-inhibition model): catalytic
  module vs DNA centre-of-mass distance / contact area (DNA retained <-> dissociated);
  the PMF along this CV = retention free energy.

> Adjust --sel1/--sel2 per system residue numbering (the 4DQY multi-chain setup
> needs segid/chain qualifiers).

## 3. Full system matrix (recommended)
Representative inhibitors x system, each >=300 ns (replicates allowed):
- system 1: veliparib / olaparib / niraparib / talazoparib / saruparib / apo
- system 2: the same representative subset (4DQY-based + docking; see the NOTE in
  build_system.py for the Zn/DNA force-field treatment)

## 4. Wall-clock estimates (RTX 4090, OpenMM mixed precision)
- system 1 (~40-50k atoms): ~120-180 ns/day -> 300 ns ~ 1.7-2.5 days/system
- system 2 (~120-150k atoms incl. DNA): ~50-80 ns/day -> 300 ns ~ 4-6 days/system
- 6x system 1 + 4x system 2 ~ 30 days on one card; ~2 weeks with two cards in parallel.

## Files (curated subset)
- `build_system.py`, `build_system2.sh`, `build_system2_tleap.py` — system builders (S1 / S2 tleap inputs)
- `deploy_s1_only.sh`, `deploy_rerun6.sh`, `deploy_{nira,ruca,veli}_prot.sh`, `deploy_veli_prot_resume.sh` — per-system deployment record
- `upload_{nira,ruca,veli}_prot.sh`, `pull_ruca_prot.sh`, `fetch_veli_dcds.sh` — transfer helpers (SSH host/port configured locally, never stored here)
- `replicate_runs.sh`, `run_review_replicates.sh` — replicate batches (REV-1/REV-2)
- `environment_gamd.yml`, `environment_gamd_min.yml` — environment specifications
- `verify_env.py`, `extract_rerun_starts.py` — utilities
- The veliparib S2 CV extraction originally deployed here is archived as
  `code/pipeline/p09_extract_veli_prot_cv.py`.

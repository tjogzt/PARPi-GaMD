# ③ Multi-seed extension — run recipe (prepared locally, 2026-10-20)

Everything that does **not** need a GPU is done. The remaining step is to run the
two jobs on an AutoDL/Linux GPU box.

## What is already prepared

| file | content |
|---|---|
| `tala_s2_equil12ns.rst7` | talazoparib S2, frame 12000 of `sys2_talazoparib/output.dcd` (step 6,000,000 = 12.0 ns, equilibration end), 285,796 atoms, CV1 63.9 Å |
| `veli_s2_equil12ns.rst7` | veliparib S2, frame 1200 of `sys2_veliparib/output.dcd` (10-ps cadence; 12.0 ns), 285,788 atoms, CV1 63.6 Å |
| `restart_provenance.json` | source paths, frame indices, steps, box vectors, prmtop SHA-256 |
| `config_seed7/*.xml`, `config_seed123/*.xml` | run configurations cloned from the veliparib-protonated run (`veli_prot`), with `random-seed` set to 7 / 123 and paths pointing at `runs_multiseed/` |

The protocol is identical to the deposited re-simulations: conventional MD 1,000,000
steps + GaMD equilibration 5,000,000 steps + GaMD production 13,000,000 steps
(production starts at step 6,000,000 = 12 ns; 26 ns production). Fresh velocities
are drawn at 300 K, so the two new seeds share the starting structure but not the
random history.

## Steps on the GPU box

```bash
# 1. upload (≈20 MB rst7 + 2×48 MB prmtop + configs)
mkdir -p /root/autodl-tmp/PARPi_design/runs_multiseed/{tala_s2,veli_s2}
# talazoparib.prmtop / veliparib.prmtop: copy the deposited prmtop into each dir,
# named tala.prmtop / veli.prmtop (paths already point there)
cp <local>/data/restarts/tala_s2_equil12ns.rst7  .../runs_multiseed/tala_s2/
cp <local>/data/restarts/veli_s2_equil12ns.rst7  .../runs_multiseed/veli_s2/

# 2. run (same launcher used for veli_prot; 2 jobs in parallel, one per GPU)
python gamd.py tala_config.xml   # seed 7   or seed 123
python gamd.py veli_config.xml   # seed 7   or seed 123
```

Measured throughput for this system size on the archived veliparib-protonated run:
**736 steps/s ≈ 127 ns/day** (`veli_prot/runner.log`, 19,000,000 steps in 31.2 h).
A conservative plan is 8–15 h per job; two jobs fit on 2×RTX 4090 24 GB in parallel.

## After the runs

1. Confirm `production-start-step.txt` equals `6000000` in each output directory.
2. Recompute CV1/CV2 spans with the archived pipeline (`s55/s56` logic, 50-ps
   cadence, first-14-ns → full 26-ns window) and add the seeds to Table S17.
3. Re-run `scripts/s57_manuscript_number_audit.py` after any text change; update
   the seed-summary variant in Table S18 (three-seed means) and the SI note that
   the talazoparib/veliparib S2 systems now carry three independent seeds.

Provenance note for the manuscript: the two new seeds share the 12-ns equilibrated
structure with the archived trajectories; only the stochastic history is independent.

#!/usr/bin/env python3
"""
p05_gen_gamd_config.py — generate a gamd-openmm (gamdRunner) configuration XML

Purpose:   Emit the GaMD config XML from the official schema (lower-dual
           variant), parameterised for protein-ligand / protein-DNA-ligand
           systems: boost type, cMD/equilibration/production step counts,
           temperature and sigma0.
Author:    Tao Zhu (tjogzt@gmail.com)
Created:   2026-09-15
Inputs:    --prmtop/--rst7  Amber topology + equilibrated coordinates (p04)
Outputs:   the config XML at --out; gamdRunner then writes to --outdir
Depends:   standard library only (schema: gamd-openmm tools/gamd-openmm/)
Run:       python3 code/pipeline/p05_gen_gamd_config.py --prmtop <prmtop> \
               --rst7 <rst7> --outdir <dir> --prod-ns 300 --out <config.xml>
"""
from __future__ import annotations
import argparse
from pathlib import Path

DT_PS = 0.002  # integration timestep: 2 fs


def ns_to_steps(ns: float) -> int:
    return int(round(ns * 1000 / DT_PS))


TEMPLATE = """<?xml version="1.0" ?>
<gamd>
    <temperature>{temp}</temperature>

    <system>
        <nonbonded-method>PME</nonbonded-method>
        <nonbonded-cutoff>1.0</nonbonded-cutoff>
        <constraints>HBonds</constraints>
    </system>

    <barostat>
        <pressure>1.0</pressure>
        <frequency>25</frequency>
    </barostat>

    <run-minimization>{run_min}</run-minimization>

    <integrator>
        <algorithm>langevin</algorithm>
        <boost-type>{boost_type}</boost-type>
        <sigma0>
            <primary>{sigma0}</primary>
            <secondary>{sigma0}</secondary>
        </sigma0>
        <random-seed>{seed}</random-seed>
        <dt>{dt}</dt>
        <friction-coefficient>1.0</friction-coefficient>
        <number-of-steps>
            <conventional-md-prep>{cmd_prep}</conventional-md-prep>
            <conventional-md>{cmd}</conventional-md>
            <gamd-equilibration-prep>{eq_prep}</gamd-equilibration-prep>
            <gamd-equilibration>{eq}</gamd-equilibration>
            <gamd-production>{prod}</gamd-production>
            <averaging-window-interval>{avg_win}</averaging-window-interval>
        </number-of-steps>
    </integrator>

    <input-files>
        <amber>
            <topology>{prmtop}</topology>
            <coordinates type="rst7">{rst7}</coordinates>
        </amber>
    </input-files>

    <outputs>
        <directory>{outdir}</directory>
        <overwrite-output>True</overwrite-output>
        <reporting>
            <energy><interval>{rep_int}</interval></energy>
            <coordinates><file-type>DCD</file-type></coordinates>
            <statistics><interval>{rep_int}</interval></statistics>
        </reporting>
    </outputs>
</gamd>
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--rst7", required=True)
    ap.add_argument("--outdir", required=True, help="gamdRunner output directory")
    ap.add_argument("--out", required=True, help="path of the generated config XML")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--boost-type", default="lower-dual",
                    choices=["lower-dual", "upper-dual", "lower-total",
                             "upper-total", "lower-dihedral", "upper-dihedral"])
    ap.add_argument("--sigma0", type=float, default=6.0, help="kcal/mol")
    ap.add_argument("--seed", type=int, default=4321)
    ap.add_argument("--cmd-ns", type=float, default=2.0, help="conventional-MD statistics collection (ns)")
    ap.add_argument("--equil-ns", type=float, default=10.0, help="GaMD equilibration / boost update (ns)")
    ap.add_argument("--prod-ns", type=float, default=300.0, help="GaMD production sampling (ns)")
    ap.add_argument("--no-min", action="store_true", help="skip minimisation (system already pre-equilibrated)")
    ap.add_argument("--report-ps", type=float, default=10.0, help="reporting interval (ps)")
    args = ap.parse_args()

    rep_int = max(1, ns_to_steps(args.report_ps / 1000))
    cfg = TEMPLATE.format(
        temp=args.temp, run_min="False" if args.no_min else "True",
        boost_type=args.boost_type, sigma0=args.sigma0, seed=args.seed, dt=DT_PS,
        cmd_prep=ns_to_steps(0.2), cmd=ns_to_steps(args.cmd_ns),
        eq_prep=ns_to_steps(0.2), eq=ns_to_steps(args.equil_ns),
        prod=ns_to_steps(args.prod_ns), avg_win=max(1, ns_to_steps(args.cmd_ns) // 40),
        prmtop=args.prmtop, rst7=args.rst7, outdir=args.outdir, rep_int=rep_int,
    )
    Path(args.out).write_text(cfg)
    print(f"[OK] -> {args.out}")
    print(f"   boost={args.boost_type} T={args.temp}K  "
          f"cMD={args.cmd_ns}ns equil={args.equil_ns}ns prod={args.prod_ns}ns")
    print(f"   total steps ~{ns_to_steps(args.cmd_ns+args.equil_ns+args.prod_ns):,} (@2fs)")
    print(f"   run: gamdRunner xml {args.out}")


if __name__ == "__main__":
    main()

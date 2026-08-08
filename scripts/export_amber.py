#!/usr/bin/env python3
"""
export_amber.py — 从已建好的 OpenMM system.xml + pdb 导出 AMBER prmtop/inpcrd

build_system.py 的 parmed.openmm.load_topology 报 'NoneType' object has no attribute 'used'。
此脚本独立调试导出: 打完整 traceback + 多策略尝试, 找出能产出有效 prmtop/rst7 的方法。
"""
from __future__ import annotations
import sys
import traceback
from pathlib import Path


def main() -> int:
    xml = sys.argv[1]   # sys1_talazoparib.system.xml
    pdb = sys.argv[2]   # sys1_talazoparib.pdb
    out = Path(sys.argv[3])  # 输出前缀

    from openmm import XmlSerializer
    from openmm.app import PDBFile
    import parmed
    print("parmed version:", parmed.__version__)

    system = XmlSerializer.deserialize(Path(xml).read_text())
    pf = PDBFile(pdb)
    print(f"system forces: {[f.__class__.__name__ for f in system.getForces()]}")
    print(f"atoms: {pf.topology.getNumAtoms()}")

    # 策略1: 直接 load_topology(topology, system, xyz)
    try:
        print("\n[策略1] parmed.openmm.load_topology(topology, system, xyz=positions)")
        st = parmed.openmm.load_topology(pf.topology, system, xyz=pf.positions)
        st.save(str(out.with_suffix(".prmtop")), overwrite=True)
        st.save(str(out.with_suffix(".inpcrd")), format="rst7", overwrite=True)
        print(f"[OK-策略1] -> {out.with_suffix('.prmtop').name} + .inpcrd")
        return 0
    except Exception:
        print("[策略1 失败] 完整 traceback:")
        traceback.print_exc()

    # 策略2: 先 load_topology(topology, system), 再单独赋坐标/盒子
    try:
        print("\n[策略2] load_topology(topology, system) 后单独赋坐标")
        st = parmed.openmm.load_topology(pf.topology, system)
        st.coordinates = pf.positions
        if pf.topology.getPeriodicBoxVectors() is not None:
            import numpy as np
            from openmm import unit
            bv = pf.topology.getPeriodicBoxVectors().value_in_unit(unit.angstrom)
            bv = np.array(bv)
            lengths = np.linalg.norm(bv, axis=1)
            st.box = [lengths[0], lengths[1], lengths[2], 90.0, 90.0, 90.0]
        st.save(str(out.with_suffix(".prmtop")), overwrite=True)
        st.save(str(out.with_suffix(".inpcrd")), format="rst7", overwrite=True)
        print(f"[OK-策略2] -> {out.with_suffix('.prmtop').name} + .inpcrd")
        return 0
    except Exception:
        print("[策略2 失败] 完整 traceback:")
        traceback.print_exc()

    print("\n[全部策略失败]")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

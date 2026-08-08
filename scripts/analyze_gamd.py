#!/usr/bin/env python3
"""
analyze_gamd.py — GaMD 轨迹 CV 计算 + 重加权准备

从 GaMD 输出 (DCD 轨迹 + gamd.log 的 boost ΔV) 计算反应坐标(CV)并生成
PyReweighting 输入, 用于重加权得到 PMF/自由能 (cumulant expansion 2nd order)。

预设 CV (MDAnalysis selection, 可命令行覆盖):
  - 体系1 变构: HD 域 vs ART 催化中心 质心距离 (HD 开/闭)
  - 体系2 trapping: PARP1(催化模块) vs DNA 质心距离 (DNA 保留/解离)

用法:
  python analyze_gamd.py --top sys.prmtop --traj gamd_out/.../output.dcd \
      --gamdlog gamd_out/.../gamd.log \
      --sel1 "resid 790-810 and name CA" --sel2 "resid 888-908 and name CA" \
      --label HD_open --out analysis/sys1_talazoparib
然后 PyReweighting:
  python PyReweighting-1D.py -input cv_weights.dat -Emax 8 -cutoff 10 \
      -disc 0.1 -T 300 -job amdweight_CE
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path


def log(m: str) -> None:
    print(f"[analyze] {m}", flush=True)


def read_gamd_dv(gamdlog: Path) -> list[float]:
    """从 gamd.log 提取每帧总 boost ΔV (kcal/mol)。
    gamd-openmm 日志列含 'Boost-Energy-Potential' 等; 取所有 boost 列之和。
    修复: 从 # 开头表头行解析 boost 列索引，再应用到 Tab 分隔数据行。"""
    all_lines = [l for l in gamdlog.read_text().splitlines() if l.strip()]
    if not all_lines:
        return []

    # 1. 从 # 表头行找 boost 列
    boost_cols = None
    for l in all_lines:
        if l.startswith("#") and re.search(r"boost", l, re.I):
            parts = re.split(r"[,\s]+", l.strip())
            boost_cols = [i for i, c in enumerate(parts)
                          if re.search(r"(?i)(?<!un)boost", c)]
            break

    # 2. 解析非 # 数据行，取 boost 列之和
    dv = []
    for l in all_lines:
        if l.startswith("#"):
            continue
        # 数据行以 Tab 开头: "\t1\t5000\t..." → parts[0] == '', 跳过
        parts = re.split(r"[,\s]+", l.strip())
        try:
            vals = [float(x) for x in parts if x != '']
        except ValueError:
            continue
        if boost_cols and len(vals) > max(boost_cols):
            # boost_cols 基于 header parts (含 '#' 在 index 0), 数据 vals 已经去空
            # header parts: ['#','ntwx','nstep',...,'Total-Boost','Dihedral-Boost',...]
            # 数据 vals: [1, 5000, ..., boost_val1, boost_val2, ...]
            # → boost_cols 偏移 -1 (去掉了 '#')
            actual_cols = [c - 1 for c in boost_cols]
            dv.append(sum(vals[c] for c in actual_cols if c < len(vals)))
        else:
            # fallback: 如果数据列数与表头 boost 列不匹配，取所有疑似 boost 列
            # (不依赖表头解析的情况)
            dv.append(vals[-1])  # 兜底
    return dv


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--gamdlog", required=True)
    ap.add_argument("--sel1", required=True, help="CV 原子组1 (MDAnalysis 选择)")
    ap.add_argument("--sel2", required=True, help="CV 原子组2")
    ap.add_argument("--label", default="CV")
    ap.add_argument("--out", required=True, help="输出前缀")
    ap.add_argument("--temp", type=float, default=300.0)
    args = ap.parse_args()

    try:
        import MDAnalysis as mda
        import numpy as np
        from MDAnalysis.analysis.distances import distance_array
    except ImportError as e:
        sys.exit(f"需 MDAnalysis/numpy: {e}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    log(f"载入 {args.top} + {args.traj}")
    u = mda.Universe(args.top, args.traj)
    g1, g2 = u.select_atoms(args.sel1), u.select_atoms(args.sel2)
    if len(g1) == 0 or len(g2) == 0:
        sys.exit(f"选择为空: sel1={len(g1)} sel2={len(g2)}")
    log(f"CV={args.label}: |sel1|={len(g1)} |sel2|={len(g2)} 质心距离")

    cv = []
    for ts in u.trajectory:
        d = np.linalg.norm(g1.center_of_mass() - g2.center_of_mass())
        cv.append(d)
    cv = np.array(cv)
    log(f"CV: n={len(cv)} 范围=[{cv.min():.2f},{cv.max():.2f}] Å 均值={cv.mean():.2f}")

    dv = read_gamd_dv(Path(args.gamdlog))
    log(f"boost ΔV: n={len(dv)}")
    n = min(len(cv), len(dv)) if dv else len(cv)
    if dv:
        dv = np.array(dv[:n]); cv = cv[:n]
    else:
        log("  [warn] 未解析到 ΔV, 仅输出 CV (重加权需手动核对 gamd.log 列)")
        dv = np.zeros(n)

    # PyReweighting 1D 输入: CV 与对应 ΔV(kcal/mol)
    cvw = out.with_name(out.name + "_cv_weights.dat")
    np.savetxt(cvw, np.column_stack([cv, dv]), fmt="%.4f",
               header=f"{args.label}(A)  dV(kcal/mol)")
    np.save(out.with_name(out.name + "_cv.npy"), cv)
    log(f"[OK] -> {cvw.name} (喂 PyReweighting-1D.py -job amdweight_CE)")
    print("\n后续重加权 (clone github.com/MiaoLab20/PyReweighting):")
    print(f"  python PyReweighting-1D.py -input {cvw.name} -T {args.temp} "
          f"-Emax 8 -cutoff 10 -disc 0.1 -job amdweight_CE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

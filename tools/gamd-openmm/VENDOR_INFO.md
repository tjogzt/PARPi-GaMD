# gamd-openmm (vendored GaMD engine)

- Upstream: gamd-openmm package by the Miao lab (University of Kansas);
  https://github.com/MiaoLab/gamd-openmm  (see AUTHORS / LICENSE, MIT).
- Vendored copy = the engine-of-record used for all S2 production runs.
- Base commit: 949e877 ("Merge pull request #39 ... pytest-compatibility").
- Local modifications vs upstream base (do not drop when updating):
  * gamd/gamdSimulation.py — inpcrd file handling refactor
    (positions/boxVectors read via an explicit inpcrd object).
  * gamdRunner — shebang python -> python3.
- Deploy scripts invoke: tools/gamd-openmm/gamdRunner xml <config.xml>
  (optionally -p CUDA/CPU; default CUDA).
- Run without install: PYTHONPATH=<this dir> python3 <this dir>/gamdRunner xml config.xml

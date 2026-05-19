# Sodium-Ion Battery Electrolyte Dataset for MLIP Benchmarking

Equilibrated atomistic configurations of Na-ion (and reference Li⁺ / K⁺) battery electrolytes in glyme and carbonate solvents, paired with experimental densities at 298.2 K. Intended as a benchmark for evaluating machine-learning interatomic potentials (MLIPs) — most directly, NPT-predicted density vs. experiment.

Configurations were equilibrated under NPT with OPLS-AA and the OMol25-trained UMA potential ([FAIRChem](https://github.com/facebookresearch/fairchem)).

> **Prediction and Experimental Verification of Electrolyte Solvation Structure from an OMol25-Trained Interatomic Potential**
> N. Kumar, J. Lai, C. S. Mezerkor, J. Wang, K. M. Wiaderek, J. D. Bazak, S. M. Blau, E. J. Crumlin (2026).
> arXiv: [2603.20183](https://arxiv.org/abs/2603.20183)

**If you use this dataset, please [cite the paper](#citation).**

---

## Contents

```
sodium-ion-battery-electrolyte-dataset/
├── README.md
├── LICENSE                           # CC-BY-4.0
├── CITATION.cff
├── CITATION.bib
├── .gitignore
├── data/
│   ├── *.xyz                         # 25 equilibrated configurations (extended XYZ)
│   └── na_ion_salt_concentrations_and_dens.csv
└── examples/
    └── load_with_ase.py
```

### Systems (25 configurations, all at 298.2 K)

| Category | Systems |
|---|---|
| **Pure solvents** (6) | DME (C₄H₁₀O₂), DEGDME (C₆H₁₄O₃), TEGDME (C₁₀H₂₂O₅), PC (C₄H₆O₃), DEG (C₄H₁₀O₃), DMC (C₃H₆O₃) |
| **0.1 M electrolytes** (8) | NaPF₆, NaOTf in DME, DEGDME, TEGDME, PC |
| **0.5 M electrolyte** (1) | NaPF₆ in DME |
| **1.0 M electrolytes** (10) | NaPF₆, NaOTf in DME / DEGDME / TEGDME / PC; NaTFSI in DME; KPF₆ in DME |

---

## File formats

### Configurations (`data/*.xyz`)

[Extended XYZ](https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz), ASE-compatible:

```
2136
Lattice="28.5962 0 0 0 28.5962 0 0 0 28.5962" Properties=species:S:1:pos:R:3 pbc="T T T"
H  0.89  24.89  11.55
...
```

Each file is a single equilibrated snapshot. Some files additionally carry per-atom momenta and forces, and per-frame energy/stress in the header; `ase.io.read` handles both transparently.

### Metadata (`data/na_ion_salt_concentrations_and_dens.csv`)

| Column | Description |
|---|---|
| `filename` | Configuration filename in `data/` |
| `concentration_M` | Effective concentration (mol L⁻¹) from n_salt / cell volume |
| `n_salt` | Number of salt formula units in the cell |
| `box_a`, `box_b`, `box_c` | Cell vectors (Å, orthorhombic) |
| `volume_ang3` | Cell volume (Å³) |
| `exp_density_g_cm3` | Experimental density at 298.2 K (g cm⁻³), from Table S1 of the paper |

---

## Quick start

```bash
pip install ase numpy pandas
```

```python
import pandas as pd
from ase.io import read

df = pd.read_csv("data/na_ion_salt_concentrations_and_dens.csv")
for _, row in df.iterrows():
    atoms = read(f"data/{row.filename}")
    # atoms.calc = MyMLIPCalculator(...)
    # run NPT at 298.2 K, average <density>, compare to row.exp_density_g_cm3
```

See `examples/load_with_ase.py` for a fuller example.

---

## Citation

```bibtex
@misc{kumar2026predictionexperimentalverificationelectrolyte,
  title         = {Prediction and Experimental Verification of Electrolyte Solvation Structure from an OMol25-Trained Interatomic Potential},
  author        = {Nitesh Kumar and Jianwei Lai and Casey S. Mezerkor and Jiaqi Wang and Kamila M. Wiaderek and J. David Bazak and Samuel M. Blau and Ethan J. Crumlin},
  year          = {2026},
  eprint        = {2603.20183},
  archivePrefix = {arXiv},
  primaryClass  = {physics.chem-ph},
  url           = {https://arxiv.org/abs/2603.20183}
}
```

---

## License

[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). See [`LICENSE`](LICENSE) for the full text.

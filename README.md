# Sodium-Ion Battery Electrolyte Dataset for MLIP Benchmarking

A curated collection of equilibrated atomistic configurations of sodium-ion (and reference Li⁺ / K⁺) battery electrolyte solutions in glyme and carbonate solvents, paired with **experimental density measurements at 298.2 K**. The dataset is intended as a **benchmark for evaluating machine-learning interatomic potentials (MLIPs)** on liquid-phase electrolyte properties — most directly, NPT-predicted density vs. experiment, but also as starting points for structural (RDF, coordination, S(Q)) and dynamic benchmarks.

The configurations themselves were equilibrated using the OMol25-trained UMA potential ([FAIRChem](https://github.com/facebookresearch/fairchem)) under NPT conditions and accompany the paper:

> **Prediction and Experimental Verification of Electrolyte Solvation Structure from an OMol25-Trained Interatomic Potential**
> Nitesh Kumar, Jianwei Lai, Casey S. Mezerkor, Jiaqi Wang, Kamila M. Wiaderek, J. David Bazak, Samuel M. Blau, Ethan J. Crumlin (2026).
> arXiv: [2603.20183](https://arxiv.org/abs/2603.20183)

**If you use this dataset, please [cite the paper](#citation).**

---

## What this dataset is for

The primary intended use is **MLIP evaluation and benchmarking**:

- Run NPT MD with a candidate MLIP (UMA, MACE, Orb, SevenNet, ORB, ESEN, etc.) starting from the supplied configurations, equilibrate at 298.2 K, and compare the predicted density against `exp_density_g_cm3` in the metadata CSV.
- Use the equilibrated configurations as inputs for structural analyses (RDF, coordination-number distributions, structure factors S(Q)) and compare against the experimental observables reported in the accompanying paper.
- Use the cell vectors and atom counts in `data/na_ion_salt_concentrations_and_dens.csv` to set up matched simulations at the correct effective concentration.

A reference baseline (UMA-OMol vs. experiment, with SevenNet-OMat and Orb-OMat comparisons for a subset) is reported in Table S1 of the paper.

---

## Contents

```
sodium-ion-battery-electrolyte-dataset/
├── README.md
├── LICENSE                           # CC-BY-4.0
├── CITATION.cff                      # GitHub citation metadata
├── CITATION.bib                      # BibTeX entry
├── .gitignore
├── data/
│   ├── *.xyz                         # 25 equilibrated configurations (extended XYZ)
│   └── na_ion_salt_concentrations_and_dens.csv
└── examples/
    └── load_with_ase.py              # Minimal loading + benchmarking example
```

### Systems included (25 configurations)

| Category | Systems |
|---|---|
| **Pure solvents** (6) | DME (C₄H₁₀O₂), diglyme / DEGDME (C₆H₁₄O₃), **tetraglyme / TEGDME (C₁₀H₂₂O₅)**, propylene carbonate / PC (C₄H₆O₃), diethylene glycol / DEG (C₄H₁₀O₃), dimethyl carbonate / DMC (C₃H₆O₃) |
| **0.1 M electrolytes** (8) | NaPF₆, NaOTf in DME, DEGDME, TEGDME, PC |
| **0.5 M electrolyte** (1) | NaPF₆ in DME |
| **1.0 M electrolytes** (10) | NaPF₆, NaOTf in DME / DEGDME / TEGDME / PC; NaTFSI in DME; KPF₆ in DME |

All configurations: cubic / orthorhombic periodic cells, equilibrated under NPT at 298.2 K.

> ⚠️ **Solvent abbreviation note.** "TEGDME" in this dataset refers to **tetraglyme (G4, C₁₀H₂₂O₅)** — that is, CH₃O(CH₂CH₂O)₄CH₃. The abbreviation is genuinely ambiguous in the literature (some authors use it for triglyme/G3); when comparing to other published values, ensure the experimental reference is also for the tetraglyme variant. The filename prefix `tgdme_` likewise refers to tetraglyme.

---

## File formats

### Configurations (`data/*.xyz`)

[Extended XYZ format](https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz), ASE-compatible. The second line encodes the periodic lattice and column schema:

```
2136
Lattice="28.5962 0 0 0 28.5962 0 0 0 28.5962" Properties=species:S:1:pos:R:3 pbc="T T T"
H  0.89  24.89  11.55
...
```

Each file is a **single equilibrated snapshot** (not a multi-frame trajectory). Some files additionally carry per-atom momenta and forces, and per-frame energy/stress in the header (extended `Properties` schema). ASE's `read("file.xyz")` handles both variants transparently.

### Metadata (`data/na_ion_salt_concentrations_and_dens.csv`)

| Column | Description |
|---|---|
| `filename` | Configuration filename in `data/` |
| `concentration_M` | Effective salt concentration (mol L⁻¹), computed from cell volume and salt count |
| `n_salt` | Number of salt formula units in the cell |
| `box_a`, `box_b`, `box_c` | Cell vectors in Å (orthorhombic) |
| `volume_ang3` | Cell volume in Å³ |
| `exp_density_g_cm3` | Experimental density at 298.2 K in g cm⁻³ |

---

## Quick start

### Install

```bash
pip install ase numpy pandas
```

### Benchmark loop (sketch)

```python
import pandas as pd
from ase.io import read

df = pd.read_csv("data/na_ion_salt_concentrations_and_dens.csv")
for _, row in df.iterrows():
    atoms = read(f"data/{row.filename}")
    # 1) attach your MLIP calculator:  atoms.calc = MyMLIPCalculator(...)
    # 2) run NPT at 298.2 K to equilibrium
    # 3) record <density> averaged over the NPT trajectory
    # 4) compare against row.exp_density_g_cm3
```

A more complete example (composition counts, density estimate from the starting cell, simple benchmarking scaffold) is in `examples/load_with_ase.py`.

---

## Important caveats — please read before benchmarking

1. **TEGDME = tetraglyme (G4, C₁₀H₂₂O₅)**, not triglyme. See the note in the systems table above. The pure-tetraglyme box contains 62 solvent molecules; corresponding electrolyte boxes contain 62 tetraglyme molecules plus salt.

2. **Single equilibrated snapshots, not trajectories.** Each `.xyz` is one configuration after NPT equilibration with the OMol25-trained UMA potential. They are intended as starting points for the *user's own* MLIP NPT runs; the density benchmark requires the user to re-equilibrate with their MLIP and time-average. For full trajectory data, please contact the corresponding author.

3. **Reference experimental densities (`exp_density_g_cm3`)** are the values reported in Table S1 of the accompanying paper.

4. **0.5 M NaPF₆ / DME** has no matched experimental density in our tabulation; the CSV value (0.93452 g cm⁻³) is a carry-over from an earlier internal estimate and should be treated cautiously when used as a benchmark target.

5. **Concentrations are effective** — they are computed from the actual number of formula units divided by the cell volume, so the 0.1 M / 1 M labels in filenames are nominal targets. Always use `concentration_M` from the CSV for quantitative comparisons.

---

## Citation

If you use this dataset in any form, please cite:

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

A `CITATION.cff` file is provided so GitHub's "Cite this repository" widget renders the entry automatically, and `CITATION.bib` contains the same BibTeX for convenience.

---

## License

The dataset (configurations, metadata CSV, derived files) and accompanying scripts are released under [**Creative Commons Attribution 4.0 International (CC-BY-4.0)**](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the material for any purpose, including commercially, as long as you give appropriate credit by citing the paper above. See [`LICENSE`](LICENSE) for the full text.

---

## Contact

- **Nitesh Kumar** — Materials Sciences Division, Lawrence Berkeley National Laboratory
- Website: [kmnitesh05.github.io](https://kmnitesh05.github.io)
- GitHub: [@KMNitesh05](https://github.com/KMNitesh05)

For questions, benchmark contributions, bug reports, or full trajectory requests, please open an [issue](https://github.com/KMNitesh05/sodium-ion-battery-electrolyte-dataset/issues) or contact the corresponding author.

## Acknowledgments

This work was carried out at Lawrence Berkeley National Laboratory. Simulations used resources of the National Energy Research Scientific Computing Center (NERSC) under allocation `m4292`. The OMol25-trained UMA potential is developed and maintained by the [FAIRChem](https://github.com/facebookresearch/fairchem) team.

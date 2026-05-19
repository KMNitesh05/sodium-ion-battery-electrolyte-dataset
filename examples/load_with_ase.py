#!/usr/bin/env python3
"""
Minimal example: load all configurations from the dataset with ASE, verify
that compositions and starting-cell densities make sense, and sketch the
shape of an MLIP benchmarking loop.

Requires:  ase, numpy, pandas

Run from the repo root:
    python examples/load_with_ase.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from ase.units import mol  # for atomic mass -> g/mol bookkeeping

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
CSV_PATH = DATA_DIR / "na_ion_salt_concentrations_and_dens.csv"


def relabel_mislabeled_sodium(atoms, filename):
    """For napf6_dme_1M.xyz and naotf_dme_1M.xyz, Na atoms were written as N.

    Remap them so MLIPs and analysis tools see the correct species. See the
    README ('Important caveats') for context.
    """
    if filename in {"napf6_dme_1M.xyz", "naotf_dme_1M.xyz"}:
        atoms.symbols = ["Na" if s == "N" else s for s in atoms.symbols]
    return atoms


def starting_density_g_cm3(atoms):
    """Density of the starting (pre-MLIP-NPT) configuration, in g/cm^3."""
    total_mass_amu = atoms.get_masses().sum()
    volume_ang3 = atoms.get_volume()
    # 1 amu = 1/N_A g; 1 Å^3 = 1e-24 cm^3
    return (total_mass_amu / mol) / (volume_ang3 * 1e-24)


def main():
    df = pd.read_csv(CSV_PATH)

    print(f"{'filename':<32}{'n_atoms':>9}{'rho_start':>12}{'rho_exp':>10}{'Δ':>10}")
    print("-" * 73)

    rows = []
    for _, row in df.iterrows():
        atoms = read(DATA_DIR / row.filename)
        atoms = relabel_mislabeled_sodium(atoms, row.filename)

        rho_start = starting_density_g_cm3(atoms)
        delta = rho_start - row.exp_density_g_cm3
        rows.append({"filename": row.filename, "rho_start": rho_start,
                     "rho_exp": row.exp_density_g_cm3, "delta": delta})

        print(f"{row.filename:<32}{len(atoms):>9}"
              f"{rho_start:>12.4f}{row.exp_density_g_cm3:>10.4f}{delta:>+10.4f}")

    out = pd.DataFrame(rows)
    mae = out["delta"].abs().mean()
    print(f"\nMAE(starting density − experiment) = {mae:.4f} g/cm^3")
    print("\nThe `rho_start` values above are NOT MLIP-equilibrated densities — they")
    print("are the densities of the supplied starting configurations. To benchmark")
    print("your MLIP, run NPT at 298.2 K from each configuration until the density")
    print("converges, then compare the time-averaged density against `rho_exp`.")
    print("\nSketch of a benchmarking loop:")
    print("""
    from ase.md.npt import NPT
    from ase import units

    for _, row in df.iterrows():
        atoms = read(DATA_DIR / row.filename)
        atoms = relabel_mislabeled_sodium(atoms, row.filename)
        atoms.calc = MyMLIPCalculator(...)                # <-- your MLIP here

        dyn = NPT(atoms,
                  timestep=0.5*units.fs,
                  temperature_K=298.2,
                  externalstress=0.0,
                  ttime=25*units.fs,
                  pfactor=(75*units.fs)**2 * units.GPa)
        # equilibrate, then sample density over a production segment
        ...
    """)


if __name__ == "__main__":
    main()

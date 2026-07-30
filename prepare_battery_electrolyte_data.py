#!/usr/bin/env python3
"""
Package the Na-ion battery electrolyte dataset for ML-PEG.

ML-PEG's molecular dynamics benchmarks read reference values straight off the
structure files (``atoms.info``), rather than from a side-car table. This script
therefore takes the configurations and the metadata CSV from

    https://github.com/KMNitesh05/sodium-ion-battery-electrolyte-dataset

and writes a single zip in which every extended-XYZ file carries its own
``exp_density`` and ``exp_temperature``, plus provenance metadata used for plot
hover text.

Run from a local clone of the dataset repository::

    python prepare_battery_electrolyte_data.py \
        --dataset-dir /path/to/sodium-ion-battery-electrolyte-dataset \
        --out-dir ./staging

The resulting ``battery_electrolyte_densities.zip`` should be committed to
``data/`` in the dataset repository, where ``download_github_data`` can fetch it.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
import zipfile

from ase import Atoms, units
from ase.io import read, write

# Every configuration in the dataset was equilibrated at this temperature, and
# the experimental densities in the CSV are quoted at the same temperature.
EXP_TEMPERATURE_K = 298.2

# amu / A^3 -> g / cm^3
AU_TO_G_CM3 = 1e24 / units.mol

# Name of the directory created inside the zip. ``download_github_data`` extracts
# into a flat cache directory, so the benchmark addresses its data as
# <cache>/BENCHMARK_NAME/STRUCT_SUBDIR.
BENCHMARK_NAME = "battery_electrolyte_densities"
STRUCT_SUBDIR = "structures"

# Per-atom arrays and info keys left over from the OPLS-AA / UMA equilibration
# runs. These are dropped so the benchmark inputs cannot be confused with
# reference values produced by another potential. Momenta are deliberately kept:
# they give the NPT run a sensible starting velocity distribution.
STALE_ARRAYS = ("forces", "energies", "stresses")
STALE_INFO = ("energy", "free_energy", "stress", "virial", "dipole", "magmom")

# Volume agreement tolerance between the CSV bookkeeping and the actual cell.
VOLUME_RTOL = 1e-3


def density_g_cm3(atoms: Atoms) -> float:
    """
    Get the density of a periodic structure in g/cm^3.

    Parameters
    ----------
    atoms
        ASE Atoms object of the periodic system.

    Returns
    -------
    float
        Density in g/cm^3.
    """
    return AU_TO_G_CM3 * atoms.get_masses().sum() / atoms.get_volume()


def classify(row: dict[str, str]) -> str:
    """
    Assign a system to a coarse category used for plot symbols.

    Parameters
    ----------
    row
        One row of the metadata CSV.

    Returns
    -------
    str
        Category label: "pure solvent", "0.1 M", "0.5 M" or "1.0 M".
    """
    n_salt = int(row["n_salt"])
    if n_salt == 0:
        return "pure solvent"
    concentration = float(row["concentration_M"])
    if concentration < 0.3:
        return "0.1 M"
    if concentration < 0.8:
        return "0.5 M"
    return "1.0 M"


def split_name(filename: str) -> tuple[str, str]:
    """
    Split a dataset filename into salt and solvent labels.

    Filenames follow ``<salt>_<solvent>_<concentration>.xyz`` for electrolytes
    and ``<solvent>_pure.xyz`` for neat solvents.

    Parameters
    ----------
    filename
        Configuration filename, e.g. "napf6_diglyme_1M.xyz".

    Returns
    -------
    tuple[str, str]
        Salt label ("none" for pure solvents) and solvent label, lower-cased.
    """
    stem = Path(filename).stem.lower()
    parts = stem.split("_")
    if parts[-1] == "pure":
        return "none", "_".join(parts[:-1])
    return parts[0], "_".join(parts[1:-1])


def read_metadata(csv_path: Path) -> list[dict[str, str]]:
    """
    Read the dataset metadata CSV.

    Parameters
    ----------
    csv_path
        Path to na_ion_salt_concentrations_and_dens.csv.

    Returns
    -------
    list[dict[str, str]]
        One dictionary per system, ordered as in the file.
    """
    with csv_path.open(encoding="utf8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    required = {"filename", "exp_density_g_cm3", "concentration_M", "n_salt"}
    missing = required - set(rows[0])
    if missing:
        raise KeyError(f"{csv_path} is missing columns: {sorted(missing)}")

    return rows


def annotate(
    dataset_dir: Path, out_dir: Path, *, strict: bool = True
) -> list[dict[str, object]]:
    """
    Write annotated extended-XYZ files and return a manifest.

    Parameters
    ----------
    dataset_dir
        Root of a clone of the dataset repository.
    out_dir
        Directory to stage the annotated structures and manifest in.
    strict
        Whether to raise if a cell volume disagrees with the CSV, or if the
        density of the supplied snapshot is implausibly far from experiment.
        Default is True.

    Returns
    -------
    list[dict[str, object]]
        Manifest rows, ordered by ``system_id``.
    """
    data_dir = dataset_dir / "data"
    rows = read_metadata(data_dir / "na_ion_salt_concentrations_and_dens.csv")

    struct_dir = out_dir / BENCHMARK_NAME / STRUCT_SUBDIR
    struct_dir.mkdir(parents=True, exist_ok=True)

    by_filename = {row["filename"]: row for row in rows}

    # ``sorted`` fixes the mapping from --system-id to system. It is applied to
    # the CSV filenames here and to the extracted files in the calc script, so
    # the two orderings are guaranteed to agree.
    manifest: list[dict[str, object]] = []
    problems: list[str] = []

    for system_id, filename in enumerate(sorted(by_filename)):
        row = by_filename[filename]
        source = data_dir / filename
        if not source.exists():
            raise FileNotFoundError(f"{source} listed in CSV but not found on disk")

        atoms = read(source)
        atoms.calc = None

        for key in STALE_ARRAYS:
            if key in atoms.arrays:
                del atoms.arrays[key]
        for key in STALE_INFO:
            atoms.info.pop(key, None)

        if not all(atoms.pbc):
            raise ValueError(f"{filename} is not fully periodic")

        csv_volume = float(row["volume_ang3"])
        cell_volume = atoms.get_volume()
        if abs(cell_volume - csv_volume) > VOLUME_RTOL * csv_volume:
            problems.append(
                f"{filename}: cell volume {cell_volume:.2f} A^3 disagrees with "
                f"CSV value {csv_volume:.2f} A^3"
            )

        exp_density = float(row["exp_density_g_cm3"])
        initial_density = density_g_cm3(atoms)
        # The snapshots are already NPT-equilibrated, so a large gap here means
        # a structure/row mismatch rather than a model error.
        if abs(initial_density - exp_density) > 0.15:
            problems.append(
                f"{filename}: snapshot density {initial_density:.3f} g/cm^3 is far "
                f"from experiment {exp_density:.3f} g/cm^3"
            )

        salt, solvent = split_name(filename)
        category = classify(row)

        atoms.info.update(
            {
                "system": Path(filename).stem,
                "exp_density": exp_density,
                "exp_temperature": EXP_TEMPERATURE_K,
                "concentration_M": float(row["concentration_M"]),
                "n_salt": int(row["n_salt"]),
                "salt": salt,
                "solvent": solvent,
                "category": category,
                # Overall charge-neutral: every salt contributes a cation and an
                # anion. Set explicitly for charge/spin-aware calculators.
                "charge": 0,
                "spin": 1,
                "reference": "arXiv:2603.20183",
            }
        )

        write(struct_dir / filename, atoms, format="extxyz")

        manifest.append(
            {
                "system_id": system_id,
                "filename": filename,
                "system": Path(filename).stem,
                "n_atoms": len(atoms),
                "salt": salt,
                "solvent": solvent,
                "category": category,
                "concentration_M": float(row["concentration_M"]),
                "exp_density_g_cm3": exp_density,
                "snapshot_density_g_cm3": round(initial_density, 5),
                "exp_temperature_K": EXP_TEMPERATURE_K,
                "has_momenta": bool(atoms.has("momenta")),
            }
        )

    if problems:
        message = "Consistency problems:\n  " + "\n  ".join(problems)
        if strict:
            raise ValueError(message)
        print(f"WARNING: {message}", file=sys.stderr)

    manifest_path = out_dir / BENCHMARK_NAME / "manifest.csv"
    with manifest_path.open("w", encoding="utf8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    return manifest


def make_zip(out_dir: Path) -> Path:
    """
    Zip the staged benchmark directory.

    Parameters
    ----------
    out_dir
        Directory containing the staged BENCHMARK_NAME directory.

    Returns
    -------
    Path
        Path to the created zip file.
    """
    root = out_dir / BENCHMARK_NAME
    zip_path = out_dir / f"{BENCHMARK_NAME}.zip"

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=str(path.relative_to(out_dir)))

    return zip_path


def main() -> None:
    """Package the dataset and report a system_id table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        required=True,
        help="Local clone of sodium-ion-battery-electrolyte-dataset.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("staging"),
        help="Directory to stage annotated structures and the zip in.",
    )
    parser.add_argument(
        "--no-strict",
        action="store_true",
        help="Warn instead of failing on consistency problems.",
    )
    args = parser.parse_args()

    manifest = annotate(args.dataset_dir, args.out_dir, strict=not args.no_strict)
    zip_path = make_zip(args.out_dir)

    header = f"{'id':>3}  {'system':<28} {'atoms':>6} {'exp':>7} {'snapshot':>9}"
    print(header)
    print("-" * len(header))
    for row in manifest:
        print(
            f"{row['system_id']:>3}  {row['system']:<28} {row['n_atoms']:>6} "
            f"{row['exp_density_g_cm3']:>7.3f} {row['snapshot_density_g_cm3']:>9.3f}"
        )

    print(f"\n{len(manifest)} systems packaged -> {zip_path}")
    print(
        "Set N_SYSTEMS in calc_battery_electrolyte_densities.py to "
        f"{len(manifest)} if it differs."
    )


if __name__ == "__main__":
    main()

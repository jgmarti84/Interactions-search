# Baseline Snapshot

Recorded state of `Interactions_search.py` before any refactoring begins.
This document is a regression anchor — use it to verify that refactoring steps preserve behavior.

---

## Version info

- File: `Interactions_search.py` — 1 214 lines
- Commit: see `git log -1 Interactions_search.py`
- Python: 3.11+
- Key dependencies: biopython, rdkit, pandas, numpy, pyyaml

---

## CLI interface

```
python Interactions_search.py [-h] (-x COMPLEX.pdb | -r RECEPTOR.pdb)
                               -c CHAIN
                               [-l LIGAND.pdb [LIGAND.pdb ...]]
                               [-n LIG_NAME]
                               [-f FORCE_LIG [FORCE_LIG ...]]
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `-x / --complex` | mutually exclusive with `-r` | Combined PDB (protein + ligand); auto-split |
| `-r / --receptor_pdb` | mutually exclusive with `-x` | Pre-separated receptor PDB |
| `-l / --ligand_input` | required with `-r` | One or more ligand PDB files |
| `-c / --chain_receptor` | always required | Protein chain letter (e.g. `A`) |
| `-n / --lig_name` | optional | HETATM residue name to select when complex has multiple groups |
| `-f / --force_ligand` | optional | Residue names stored as ATOM to treat as ligand |

### Usage modes

**Mode 1 — Separate PDB files**
```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A
```

**Mode 2 — Batch (multiple ligands)**
```bash
python Interactions_search.py -r protein.pdb -l lig1.pdb lig2.pdb lig3.pdb -c A
```

**Mode 3 — Complex PDB**
```bash
python Interactions_search.py -x complex.pdb -c A             # auto-select if single HETATM group
python Interactions_search.py -x complex.pdb -c A -n TF3      # select named group
python Interactions_search.py -x complex.pdb -c A -f TF3      # ligand stored as ATOM records
```

---

## Analysis pipeline (function call order)

1. `main()` — parses args, resolves pairs, calls `carga_variables()`, loops over `analyze_pair()`
2. `split_pdb()` — if `--complex`: splits into protein PDB + one HETATM PDB per residue name
3. `validate_inputs()` — checks files exist and chain is present in receptor
4. `analyze_pair()` — orchestrates steps 5–14 below
5. `remove_bias()` — strips CM atoms from ligand PDB; saves `<ligand>_old.pdb`
6. `extract_coords_from_pdb()` — reads ligand PDB into list of tuples + centre of mass
7. `active_site_residues()` — BioPython: residues within `centroid_distance` Å of ligand CM
8. `search_hot_points()` — RDKit SMARTS: ligand acceptors and donors; generates PNGs if `ligand_plot: Yes`
9. `search_rings()` — RDKit ring detection filtered by planarity RMSD
10. `Coordenadas_interes_receptor()` — maps acceptors/donors/aromatic centroids from YAML onto active site atoms
11. Contact searches (vectorised numpy distances):
    - H-bonds: inline in `analyze_pair` — `receptor_points × acceptor/donor atoms`
    - Aromatic: `residuos_contacto()` — ring centroid distances + `Interaccion_Aromatica()`
    - Hydrophobic: `search_hydrophobic()` — apolar C atoms, deduplicated to 1 row per receptor residue
    - Salt bridges: `search_salt_bridges()`
    - Pi-cation: `search_pi_cation()`
12. Angle validation — `angle_three_points()` for H-bonds; `aromatic_angle()` for aromatics
13. CSV output — three CSV files (all / threshold / true)
14. VMD output — `scripting_vmd()` + `scripting_vmd_hydrophobic()` if `vmd_output: Yes`
15. Cumulative append — `_append_cumulative_csv()` if `cumulative_output: Yes`
16. `shutil.rmtree` — deletes temp dir created by `split_pdb` (complex mode only)

---

## Output files (per receptor–ligand pair)

All outputs are written to `<receptor_stem>_<ligand_stem>/`.

| File | Description |
|---|---|
| `<receptor>.pdb` | Copy of receptor PDB |
| `<ligand>.pdb` | Copy of ligand PDB (post-bias-removal) |
| `<ligand>_old.pdb` | Pre-bias-removal ligand copy |
| `Interaction_*_all.csv` | All contacts found (no filter) |
| `Interaction_*_threshold.csv` | Distance-filtered contacts |
| `Interaction_*_true.csv` | Distance + angle validated contacts |
| `summary.csv` | Count per interaction type |
| `CM.csv` | Ligand centre of mass |
| `vmd_*.tcl` | VMD visualisation script (if `vmd_output: Yes`) |
| `*_acceptors.png` | Ligand acceptors highlighted (if `ligand_plot: Yes`) |
| `*_donors.png` | Ligand donors highlighted (if `ligand_plot: Yes`) |
| `*_aromatic.png` | Ligand aromatic rings highlighted (if `ligand_plot: Yes`) |

**Cumulative files** (in the working directory, not per-pair):
- `Interactions_close.csv` — one row per pair (if `cumulative_output: Yes`)
- `CM_all.csv` — centre of mass per pair (if `cumulative_output: Yes`)

---

## CSV column schema

| Column | Type | Description |
|---|---|---|
| `Pos R` | int | Receptor residue sequence number |
| `Res` | str | Residue name (e.g. `SER`, `TYR`) |
| `Atom` | str | Receptor atom name |
| `Dist` | float | Distance in Å |
| `Lig` | str | Ligand atom name or ring label (e.g. `O1`, `R1`) |
| `Type` | str | `acceptor`, `donor`, `aromatic`, `hydrophobic`, `salt_bridge`, `pi_cation` |
| `Angle` | float | Validation angle in degrees (0 for types without angle check) |
| `Interaction` | str | `'Yes'` if distance + angle criteria met, else `'No'` |

---

## Configuration

Loaded from `Interacciones_variables.yml` by `carga_variables()`.
Returns a 14-element tuple (no type safety). Key values:

| Key | Default | Used by |
|---|---|---|
| `ligand_plot` | `'Yes'` | PNG generation |
| `vmd_output` | `'Yes'` | TCL generation |
| `cumulative_output` | `'Yes'` | Append to shared CSVs |
| `Distances_Hidrogen_Bonds` | `3.2` | H-bond distance filter |
| `Distances_Aromatic` | `5.5` | Aromatic centroid distance |
| `Distances_Hidrofobica` | `4.0` | Hydrophobic C-C distance |
| `centroid_distance` | `12.0` | Active site search radius |
| `Angle_Hidrogen_Bonds_Min` | `100` | H-bond angle minimum |
| `Angle_Hidrogen_Bonds_Max` | `180` | H-bond angle maximum |
| `Ring_Planarity_RMSD_Max` | `0.15` | Aromatic planarity filter |

---

## Known limitations (pre-refactor)

- `carga_variables()` returns an untyped 14-element tuple — easy to pass in wrong order
- `cfg` dict uses mixed-language keys (`Distancia_Hidrofobica`, `Aceptores_Prot`, etc.)
- `analyze_pair()` is ~200 lines; H-bond logic is inline, not extractable without touching the function
- No type annotations anywhere in the file
- Spanish identifiers throughout: `carga_variables`, `Coordenadas_interes_receptor`, `Busqueda_Antecesor_Lig`, `Interaccion_Aromatica`, `residuos_contacto`, etc.
- `src/interactions_search/` subpackage structure exists but `.py` sources are missing
- `tests/test_smoke.py` source is missing; needs to be recreated from scratch
- `pyproject.toml` references a non-existent `Geometry.py` module
- No CI configured

---

## Regression check command

After each refactoring step, run this to verify parity:

```bash
python Interactions_search.py -r tests/fixtures/receptor_mini.pdb \
                               -l tests/fixtures/ligand_mini.pdb \
                               -c A
# then compare output CSV to golden files in tests/fixtures/golden/
```

Or via pytest:
```bash
pytest tests/test_smoke.py -v
```

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

| Argument | Required | Effect |
|---|---|---|
| `-x / --complex` | mutually exclusive with `-r` | Combined PDB (protein + ligand); auto-split into protein and one HETATM file per residue name. Water excluded. Temp files deleted after run. |
| `-r / --receptor_pdb` | mutually exclusive with `-x` | Pre-separated receptor PDB passed directly to `analyze_pair` |
| `-l / --ligand_input` | required with `-r` | One or more ligand PDB files; each creates a separate analysis pair |
| `-c / --chain_receptor` | always required | Protein chain letter (e.g. `A`); used by BioPython to scope the active-site search |
| `-n / --lig_name` | optional, `--complex` only | HETATM residue name to select when the complex contains multiple ligand groups; without it the script exits if more than one group is found |
| `-f / --force_ligand` | optional, `--complex` only | Residue name(s) stored as `ATOM` records to re-label as `HETATM` during splitting (e.g. CHARMM-prepped structures) |

### Usage modes

**Mode 1 — Separate PDB files**
```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A
```
One analysis pair. Receptor and ligand are already separate files.

**Mode 2 — Batch (multiple ligands, one receptor)**
```bash
python Interactions_search.py -r protein.pdb -l lig1.pdb lig2.pdb lig3.pdb -c A
```
`analyze_pair` is called once per ligand. Each pair writes to its own output folder.
If `cumulative_output: 'Yes'`, `Interactions_close.csv` and `CM_all.csv` in the working
directory get one row appended per pair (tagged by Receptor/Ligand, so successive runs
accumulate without overwriting).

**Mode 3 — Complex PDB**
```bash
# Single HETATM group — selected automatically
python Interactions_search.py -x complex.pdb -c A

# Multiple HETATM groups — must name one with -n
python Interactions_search.py -x complex.pdb -c A -n TF3

# Ligand saved as ATOM records, not HETATM
python Interactions_search.py -x complex.pdb -c A -f TF3
```
`split_pdb()` separates the complex into a protein PDB and one file per HETATM residue name
(water excluded). Files go to a `tempfile.mkdtemp` directory that is deleted at exit.
If multiple HETATM groups exist and `-n` is not given, the script prints the available
names and exits with code 1 without running any analysis.

---

## Analysis pipeline (function call order)

1. `main()` — parses args, resolves the list of `(receptor, ligand)` pairs, calls
   `carga_variables()` once, then loops over each pair.
2. `split_pdb()` — (`--complex` mode only) separates protein from HETATM groups;
   distributes `CONECT` records to the matching HETATM file by atom serial.
3. `validate_inputs()` — checks that both files exist, the receptor has `ATOM` records,
   and the requested chain is present. Errors are printed and the pair is skipped.
4. `analyze_pair()` — orchestrates steps 5–16 below for one receptor–ligand pair.
5. `remove_bias()` — strips lines containing ` CM ` from the ligand PDB; saves a
   `<ligand>_old.pdb` copy before modification.
6. `extract_coords_from_pdb()` — reads ligand PDB into a list of 8-tuples
   `(atom_id, atom_name, res_name, chain_id, res_seq, x, y, z)` plus the geometric
   centre of mass (mean of all atom positions).
7. `search_hot_points()` — RDKit SMARTS matching on the ligand molecule to identify
   H-bond acceptor and donor atom indices. Generates `*_acceptors.png` and
   `*_donors.png` if `ligand_plot: 'Yes'`.
   - Acceptor SMARTS: `[O;H1]`, `[O;H0]`, `[N;H1]`, `[N;H0]`, `[n]`, `[o]`, `[N+]`
   - Donor SMARTS: `[O;H]`, `[N;H2]`, `[N;H]`, `[S;H]`, `[nH]`
8. `search_rings()` — `mol.GetRingInfo().AtomRings()` filtered to rings with more than 5
   atoms and a best-fit-plane RMSD ≤ `Ring_Planarity_RMSD_Max`. RDKit's own aromaticity
   flag is not used (unreliable for PDB-sourced molecules without explicit bond orders).
   Generates `*_aromatic.png` if `ligand_plot: 'Yes'`.
9. `active_site_residues()` — BioPython: iterates over all residues in the specified
   chain; keeps those whose gravitic centre of mass is within `centroid_distance` Å of
   the ligand CM. Excludes `HOH` and the ligand residue name itself.
10. `Coordenadas_interes_receptor()` — maps YAML acceptor/donor tables onto the active
    site atoms to produce a `receptor_points` DataFrame. Adds aromatic centroids for
    TYR, PHE, TRP using `get_aromatic_coord()` → `center_aromatic_ring()`.
11. Contact searches (vectorised numpy distances):
    - **H-bond** — `residuos_contacto()` called twice: receptor donors vs. ligand
      acceptors; receptor acceptors vs. ligand donors. Proximity threshold: 4 Å (hard-coded).
    - **Aromatic** — inline loop in `analyze_pair`: each ligand ring centroid vs. each
      receptor aromatic centroid, threshold `Distances_Aromatic`.
    - **Hydrophobic** — `search_hydrophobic()`: apolar C atoms matched by SMARTS
      `[c,C;!$([C,c]~[#7,#8,#16,#15,#9,#17,#35,#53])]` vs. `_HYDROPHOBIC_ATOMS` dict.
      Multiple contacts to the same receptor residue are collapsed to one row per
      (ligand atom, residue) pair; `Atom` lists all receptor atoms comma-separated and
      `Dist` is their mean.
    - **Salt bridge** — `search_salt_bridges()`: charged SMARTS on ligand vs.
      `_SALT_POS_ATOMS`/`_SALT_NEG_ATOMS` dicts on receptor. Hard-coded threshold: 4.0 Å.
    - **π-cation** — `search_pi_cation()`: ligand ring centroids vs. ARG/LYS/HIS cation
      atoms on receptor. Hard-coded threshold: 5.0 Å.
12. Angle validation — iterates over `DF_Interacciones` by row index:
    - `acceptor` rows: donor = receptor atom; acceptor = ligand atom; antecedent =
      nearest non-H ligand atom (`Busqueda_Antecesor_Lig`). Angle via `angle_three_points`.
    - `donor` rows: donor = ligand atom; acceptor = receptor atom; antecedent looked up
      from `acceptors_antecedent` YAML table, falling back to backbone `C` if missing.
    - `aromatic` rows: `Interaccion_Aromatica()` → `aromatic_angle()` computes the angle
      between ring-plane normal vectors.
    - `hydrophobic`, `salt_bridge`, `pi_cation` rows: no angle computed (left at 0.0).
13. Final classification — sets `Interaction` column to `'Yes'` or `'No'`:
    - `acceptor`/`donor`: `dist < Distances_Hidrogen_Bonds` AND
      `Angle_Hidrogen_Bonds_Min < angle <= Angle_Hidrogen_Bonds_Max` (default 100°–180°).
    - `aromatic`: `dist < Distances_Aromatic` AND (`angle < 30°` OR `angle > 60°`).
    - `hydrophobic`, `salt_bridge`, `pi_cation`: `Interaction` is set to `'Yes'` at
      creation time in their search functions; not re-evaluated here.
14. CSV output — three files written (all interactions / distance-filtered / validated).
    Internal `LigID` column is dropped before writing.
15. VMD output — `scripting_vmd()` + `scripting_vmd_hydrophobic()` if `vmd_output: 'Yes'`.
16. Cumulative append — `_append_cumulative_csv()` for `Interactions_close.csv` and
    `CM_all.csv` if `cumulative_output: 'Yes'`.
17. Temp-dir cleanup — `shutil.rmtree(tmp_dir)` at the end of `main()` (complex mode only).

---

## Output files (per receptor–ligand pair)

All per-pair outputs are written to `<receptor_stem>_<ligand_stem>/`.

| File | Condition | Description |
|---|---|---|
| `<receptor>.pdb` | always | Copy of receptor PDB |
| `<ligand>.pdb` | always | Copy of ligand PDB (post-bias-removal) |
| `<ligand>_old.pdb` | always | Pre-bias-removal ligand copy |
| `Interaction_*_all.csv` | always | All contacts found (no filter) |
| `Interaction_*_threshold.csv` | always | Contacts where `Dist < Distances_Aromatic` |
| `Interaction_*_true.csv` | always | Contacts where `Interaction == 'Yes'` |
| `summary.csv` | always | Count per interaction type for all/threshold/true |
| `CM.csv` | always | Ligand geometric centre of mass |
| `vmd_*.tcl` | `vmd_output: 'Yes'` | VMD script for H-bond / aromatic interactions |
| `vmd_hydrophobic_*.tcl` | `vmd_output: 'Yes'` | VMD script for hydrophobic contacts |
| `*_acceptors.png` | `ligand_plot: 'Yes'` | Ligand with acceptor atoms highlighted |
| `*_donors.png` | `ligand_plot: 'Yes'` | Ligand with donor atoms highlighted |
| `*_aromatic.png` | `ligand_plot: 'Yes'` | Ligand with aromatic rings highlighted |

**Cumulative files** (working directory, not per-pair):

| File | Condition | Description |
|---|---|---|
| `Interactions_close.csv` | `cumulative_output: 'Yes'` | One row per pair: total counts + per-type counts for dist and true subsets |
| `CM_all.csv` | `cumulative_output: 'Yes'` | Ligand centre of mass, one row per pair |

---

## CSV column schema

Applies to `_all.csv`, `_threshold.csv`, and `_true.csv`.

| Column | Type | Description |
|---|---|---|
| `Pos R` | int | Receptor residue sequence number |
| `Res` | str | Residue name (e.g. `SER`, `TYR`) |
| `Atom` | str | Receptor atom name; for `hydrophobic`, comma-separated when multiple atoms of the same residue are involved |
| `Dist` | float | Distance in Å; for collapsed hydrophobic rows, the mean distance |
| `Lig` | str | Ligand atom name, ring label (e.g. `aromatic 1 (#6)`), or SMARTS-matched atom name |
| `Type` | str | `acceptor`, `donor`, `aromatic`, `hydrophobic`, `salt_bridge`, `pi_cation` |
| `Angle` | float | Validation angle in degrees; 0.0 for types without angle validation |
| `Interaction` | str | `'Yes'` if all distance and angle criteria are met, else `'No'` |

---

## Configuration (`Interacciones_variables.yml`)

Loaded by `carga_variables()`, which returns a 14-element positional tuple. `analyze_pair`
converts it to a plain `dict` keyed by the variable names below.

| Key | Default | Effect |
|---|---|---|
| `ligand_plot` | `'Yes'` | Generates PNG images of ligand acceptors, donors, and rings |
| `vmd_output` | `'Yes'` | Generates TCL scripts for VMD |
| `cumulative_output` | `'Yes'` | Appends each pair to `Interactions_close.csv` / `CM_all.csv` |
| `Distances_Hidrogen_Bonds` | `3.2` | H-bond distance cutoff (Å) for final classification |
| `Distances_Aromatic` | `5.5` | Aromatic centroid-to-centroid cutoff (Å); also used as threshold filter |
| `Distances_Hidrofobica` | `4.0` | Hydrophobic C-C distance cutoff (Å) |
| `centroid_distance` | `12.0` | Active site search radius (Å) from ligand CM |
| `Angle_Hidrogen_Bonds_Min` | `100` | H-bond angle minimum (°) |
| `Angle_Hidrogen_Bonds_Max` | `180` | H-bond angle maximum (°); `np.arccos` ceiling is 180° |
| `Ring_Planarity_RMSD_Max` | `0.15` | Maximum RMSD (Å) to best-fit plane for a ring to be considered aromatic |

The YAML also contains per-residue `acceptors`, `donors`, `acceptors_antecedent`, and
`special` tables used by `Coordenadas_interes_receptor()` and the angle validation step.

---

## Known limitations (pre-refactor)

- `carga_variables()` returns a 14-element positional tuple — callers must unpack in
  exact order; wrong-position bugs are silent.
- `cfg` dict passed to `analyze_pair` uses mixed-language keys (`Distancia_Hidrofobica`,
  `Aceptores_Prot`, `Dadores_Prot`, `Aceptot_antecedent`).
- `analyze_pair()` is ~200 lines; H-bond contact search and aromatic contact search are
  inline (no dedicated function), making them hard to test or replace independently.
- No type annotations anywhere in `Interactions_search.py`.
- Spanish identifiers throughout: `carga_variables`, `Coordenadas_interes_receptor`,
  `Busqueda_Antecesor_Lig`, `Interaccion_Aromatica`, `residuos_contacto`,
  `center_aromatic_ring` (the `x/y/z` variable names inside it), etc.
- `center_of_mass` is defined directly in `Interactions_search.py`; the `src/`
  package structure exists but contains no importable source modules yet.
- Salt bridge and π-cation thresholds (4.0 Å and 5.0 Å) are hard-coded in
  `_SALT_DIST` and `_PI_CATION_DIST` constants, not read from the YAML.
- The H-bond proximity pre-filter in `residuos_contacto` is hard-coded at 4 Å
  (`threshold_PH = 4` in `analyze_pair`), independent of `Distances_Hidrogen_Bonds`.
- No CI configured.

---

## Current test coverage

### Smoke tests — `tests/test_smoke.py`

Four tests run the full script as a subprocess against minimal fixture PDBs
(`tests/fixtures/receptor_mini.pdb` + `tests/fixtures/ligand_mini.pdb` for the normal
case; `receptor_noint.pdb` + `ligand_noint.pdb` for the no-interaction case).

| Test | Fixtures used | What it checks |
|---|---|---|
| `test_exit_code_zero` | both pairs | `returncode == 0` for both runs |
| `test_output_folder_and_all_csv_exist` | normal pair | Output folder created; at least one `Interaction_*_all.csv` present |
| `test_all_csv_has_expected_columns` | normal pair | `_all.csv` has exactly `["Pos R", "Res", "Atom", "Dist", "Lig", "Type", "Angle", "Interaction"]` |
| `test_true_csv_interactions_all_yes` | normal pair | `_true.csv` is non-empty; every row has `Interaction == 'Yes'` |

Cleanup is handled by a single `autouse=True` module-scoped fixture (`clean_module`)
that wipes both output folders and the cumulative CSVs before and after the module.
The two run fixtures (`run_output`, `run_output_noint`) perform no cleanup themselves.

**Not covered by smoke tests:** intermediate DataFrames, angle values, per-interaction-type
counts, hydrophobic/salt-bridge/π-cation specific rows, VMD script content, PNG generation,
cumulative CSV append behavior, batch mode, complex-PDB splitting.

---

## Regression check command

After each refactoring step, run the full test suite:

```bash
pytest tests/ -v
```

For a manual end-to-end check against the fixture data:

```bash
python Interactions_search.py \
    -r tests/fixtures/receptor_mini.pdb \
    -l tests/fixtures/ligand_mini.pdb \
    -c A
```

Expected: exit code 0, folder `receptor_mini_ligand_mini/` created, `_true.csv` non-empty
with all rows having `Interaction == 'Yes'`, columns matching the schema above.

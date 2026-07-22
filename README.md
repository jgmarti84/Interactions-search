# Interactions Search

Detects and classifies non-covalent interactions between a ligand and a protein from PDB files. Identifies hydrogen bonds, aromatic interactions (π-π and T-shaped), hydrophobic contacts, salt bridges, and π-cation interactions. Each contact is validated by distance and angle criteria. Outputs CSV files and TCL scripts for visualization in VMD.

---

## Files

| File | Description |
|---|---|
| `Interactions_search.py` | Main script |
| `Interacciones_variables.yml` | Distance thresholds, acceptors and donors per residue |

---

## Dependencies

```
biopython
rdkit
pandas
numpy
pyyaml
```

---

## Usage

### Mode 1 — Separate PDB files

```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A
```

### Mode 2 — Batch (multiple ligands, one receptor)

```bash
python Interactions_search.py -r protein.pdb -l lig1.pdb lig2.pdb lig3.pdb -c A
```

Each pair generates its own output folder. If `cumulative_output: 'Yes'` in the config, `Interactions_close.csv` and `CM_all.csv` (in the current directory) get one row appended per pair, tagged by `Receptor`/`Ligand`, so successive or batch runs accumulate without overwriting each other.

### Mode 3 — Complex PDB

A single PDB containing protein + ligand(s). The script splits it automatically.

```bash
# Single HETATM group (selected automatically)
python Interactions_search.py -x complex.pdb -c A

# Multiple HETATM groups (select with -n)
python Interactions_search.py -x complex.pdb -c A -n LIG
```

If multiple HETATM groups are present and `-n` is not specified, the script lists the available names and exits without analysing.

If the complex PDB has no `HETATM` records at all (ligand saved as `ATOM`, e.g. some CHARMM/AMBER-prepped structures), use `-f` to tell `split_pdb` which residue name(s) to treat as ligand:

```bash
# Ligand TF3 is stored as ATOM records, not HETATM
python Interactions_search.py -x complex.pdb -c A -f TF3

# Multiple forced ligand names; still use -n to pick one for this run
python Interactions_search.py -x complex.pdb -c A -f TF3 7FW -n TF3
```

### Arguments

| Argument | Description |
|---|---|
| `-x / --complex` | Complex PDB. Alternative to `-r`. |
| `-r / --receptor_pdb` | Receptor PDB (already separated). |
| `-l / --ligand_input` | Ligand PDB(s). Accepts one or several (batch). |
| `-c / --chain_receptor` | Protein chain (e.g. `A`). |
| `-n / --lig_name` | HETATM name when multiple groups exist in `--complex`. |
| `-f / --force_ligand` | Residue name(s) to treat as ligand even if stored as `ATOM` instead of `HETATM` in `--complex`. |

---

## Analysis Pipeline

```
Complex PDB (optional)
        │
        ▼
[1] split_pdb()
    ├── <stem>_protein.pdb      ← ATOM records
    └── <stem>_RESNAME.pdb      ← HETATM records per residue (water excluded)
        │
        ▼
[2] Ligand cleanup
    └── remove_bias()           ← removes CM atoms from ligand PDB
        │
        ▼
[3] Ligand hot-points  (RDKit + SMARTS)
    ├── H-bond acceptors:  [O;H1], [O;H0], [N;H1], [N;H0], [n], [o], [N+]
    ├── H-bond donors:     [O;H], [N;H2], [N;H], [S;H], [nH]
    └── Aromatic rings:    detected via ring_info, filtered by size > 5 and by
                           planarity (best-fit-plane RMSD ≤ Ring_Planarity_RMSD_Max)
        │
        ▼
[4] Receptor active site  (BioPython)
    └── active_site_residues()
        Residues whose centre of mass is within 12 Å of the ligand CM
        HOH and the ligand itself are excluded
        │
        ▼
[5] Receptor points of interest
    └── Coordenadas_interes_receptor()
        ├── Receptor acceptors  (from YAML table)
        ├── Receptor donors     (from YAML table)
        └── Aromatic ring centroids (TYR, PHE, TRP)
        │
        ▼
[6] Contact search  (numpy, vectorised distances)
    ├── H-bond:        lig acceptor ↔ rec donor      (threshold: Distances_Hidrogen_Bonds)
    │                  lig donor    ↔ rec acceptor    (threshold: Distances_Hidrogen_Bonds)
    ├── Aromatic:      centroid ↔ centroid             (threshold: Distances_Aromatic)
    ├── Hydrophobic:   apolar C lig ↔ apolar C rec     (threshold: Distances_Hidrofobica)
    ├── Salt bridge:   ± group lig ↔ ∓ group rec       (threshold: 4.0 Å)
    └── π-cation:      lig ring ↔ ARG/LYS/HIS rec      (threshold: 5.0 Å)
        │
        ▼
[7] Angle validation
    ├── H-bond:    D-A···Antecedent angle between 100° and 200°
    ├── Aromatic:  angle between ring planes
    │               0°–30°  → π-π (parallel / sandwich)
    │               60°–90° → T-shaped (perpendicular)
    └── Hyd. / salt / π-cat:  validated by distance only
        │
        ▼
[8] Hydrophobic pocket detection  (search_hydrophobic_pockets, independent of step 7)
    ├── Group ligand hydrophobic atoms into fragments by bond connectivity
    ├── Per fragment, collect distinct contacting receptor residues (≥ Pocket_Min_Residues)
    └── Score spatial coverage around the fragment (Coverage_R, see "Hydrophobic Pockets" below)
        │
        ▼
[9] Outputs
    ├── CSV with all raw interactions
    ├── CSV filtered by distance
    ├── CSV filtered by distance + angle  (validated interactions)
    ├── CSV of hydrophobic pocket candidates (Pockets_<rec>_<lig>.csv)
    ├── Console summary  (table of validated interactions + pocket count)
    ├── Cumulative summary Interactions_close.csv
    ├── Centre of mass    CM_all.csv
    └── TCL scripts for VMD (if vmd_output: Yes)
```

---

## PDB Splitting (`split_pdb`)

When `--complex` is used, `split_pdb()` pre-processes the PDB before analysis:

- Separates `ATOM` records (protein) from `HETATM` records (ligands/cofactors)
- Groups HETATM records by residue name (`resName`)
- **Excludes water** automatically: HOH, WAT, TIP, TIP3, SOL, DOD
- Distributes `CONECT` records to the corresponding HETATM file (by atom serial)
- Saves files to a temporary directory that is deleted after the run

Example with a PDB containing protein + ligand LIG + HEM group:

```
<tmpdir>/
├── complex_protein.pdb   ← full protein
├── complex_LIG.pdb       ← organic ligand + its CONECT records
└── complex_HEM.pdb       ← haem group
```

---

## Configuration (`Interacciones_variables.yml`)

```yaml
options:
  ligand_plot: 'Yes'         # generates PNG images of ligand acceptors, donors and rings
  vmd_output:  'Yes'         # generates TCL script for VMD visualisation
  cumulative_output: 'Yes'   # appends each pair to Interactions_close.csv / CM_all.csv

distancias:
  Distances_Hidrogen_Bonds: 3.2   # Å — H-bond threshold
  Distances_Aromatic:       5.5   # Å — centre-to-centre aromatic threshold
  Distances_Hidrofobica:    4.0   # Å — hydrophobic threshold (aligned with PLIP)
  centroid_distance:       12.0   # Å — active site search radius
  Distances_C_Simple:       1.54  # Å — C-C single bond (reference)
  Distances_C_Doble:        2.56  # Å — C=C double bond (reference)

angulos:
  Angle_Hidrogen_Bonds_Min: 100    # ° — minimum Donor-Acceptor-Antecedent angle
  Angle_Hidrogen_Bonds_Max: 180    # ° — maximum angle (180° is the geometric ceiling)

aromaticidad:
  Ring_Planarity_RMSD_Max:  0.15   # Å — max RMSD to the ring's best-fit plane to
                                   # be considered aromatic (real rings ~0.01, chair ~0.25)

acceptors:             # acceptor atoms per residue
  ALA: [O]
  TYR: [O, OH]
  ASP: [O, OD1, OD2]
  ...

donors:                # donor atoms per residue
  ALA: [N]
  ARG: [N, HNE, HH11, HH12, HH21, HH22, HE, NE, NH1, NH2]
  ...

acceptors_antecedent:  # antecedent atom of each acceptor (for angle calculation)
  TYR: {OH: CZ}
  ASP: {OD1: CG, OD2: CG}
  ...

special:               # special cases (e.g. haem group)
  HEM: [FE, 1.59]

pockets:
  min_residues:        3     # minimum distinct residues contacting the same ligand fragment
  coverage_threshold:  0.5   # max Coverage_R (0-1) to qualify as an enclosing pocket
```

---

## Hydrophobic Pockets

A single hydrophobic contact (one residue, one ligand atom) doesn't tell you whether the
ligand sits in a real, enclosing binding pocket, or just brushes past a residue on one
side. `search_hydrophobic_pockets()` (in `Interactions_search.py`) answers that question
with two independent criteria, both must hold:

1. **Multiple residues on the same ligand fragment.** Ligand hydrophobic atoms (matched by
   the same SMARTS used for `hydrophobic` contacts, `_HPHO_LIG_SMARTS`) are grouped into
   *fragments* by **bond connectivity** (RDKit's bond graph), not spatial proximity — a ring
   or a contiguous aliphatic chain that is contacted counts as one fragment. For each
   fragment, every receptor residue with at least one apolar atom within
   `Distances_Hidrofobica` Å of any atom in that fragment is collected. The fragment
   qualifies only if it has **≥ `min_residues`** distinct contacting residues (default 3).

2. **Spatial coverage around the fragment.** Having 3+ residues touching the same fragment
   isn't enough on its own — they could all be sitting on the same face of the ligand
   (a flat, superficial contact) rather than wrapping around it. `Coverage_R` measures this:
   for each contacting residue, take the unit vector from the fragment's centroid to that
   residue's centroid, then compute the magnitude of the **average of those unit vectors**.

   ```
   Coverage_R = | Σ unit_vectors | / n_residues        (0 ≤ Coverage_R ≤ 1)
   ```

   - **Coverage_R ≈ 0** — the vectors point in different directions and cancel out:
     residues surround the fragment from multiple sides → a real, enclosing pocket,
     consistent with well-defined binding sites seen in crystal structures.
   - **Coverage_R ≈ 1** — the vectors mostly point the same way: all residues are on
     the same side → a superficial contact, not an enclosing pocket, even with 3+ residues.

   A fragment is marked `Is_Pocket = Yes` only if `n_residues ≥ min_residues` **and**
   `Coverage_R < coverage_threshold` (default 0.5).

This runs independently of the per-contact `Interaction == Yes` validation in step 7 — a
fragment can have several individually-validated hydrophobic contacts and still fail the
pocket criteria (e.g. only 2 residues), or vice versa.

---

## Outputs

### CSV files

| File | Content |
|---|---|
| `<folder>/Interaction_<rec>_<lig>_all.csv` | All interactions found (no filters) |
| `<folder>/Interaction_<rec>_<lig>_threshold.csv` | Filtered by distance |
| `<folder>/Interaction_<rec>_<lig>_true.csv` | Validated by distance and angle |
| `<folder>/Pockets_<rec>_<lig>.csv` | Hydrophobic pocket candidates (see "Hydrophobic Pockets" above), one row per ligand fragment |
| `Interactions_close.csv` | Cumulative run summary, one row per pair (same content as `summary.csv`, including per-type counts) — requires `cumulative_output: 'Yes'` |
| `CM_all.csv` | Ligand centre of mass, one row per pair — requires `cumulative_output: 'Yes'` |

`Pockets_<rec>_<lig>.csv` columns:

| Column | Description |
|---|---|
| `Pocket` | Fragment id (arbitrary, stable within the run) |
| `Fragment_Atoms` | Ligand atom names in the fragment, comma-separated |
| `N_Ligand_Atoms` | Number of ligand atoms in the fragment actually in contact |
| `Residues` | Contacting receptor residues, e.g. `LEU63,VAL67,TYR129` |
| `N_Residues` | Distinct contacting residue count |
| `Coverage_R` | Spatial coverage score, 0–1 (see above); lower = more enclosing |
| `Is_Pocket` | `Yes` / `No` — whether both criteria (`N_Residues` and `Coverage_R`) are met |

Interaction CSV columns:

| Column | Description |
|---|---|
| `Pos R` | Receptor residue number |
| `Res` | Residue name (e.g. SER, TYR) |
| `Atom` | Receptor atom involved. For `hydrophobic`, when the same ligand atom contacts several atoms of the same residue, they are collapsed into one row and listed comma-separated (e.g. `CD1,CD2,CG`), with `Dist` averaged across them |
| `Dist` | Distance in Å |
| `Lig` | Ligand atom or ring involved |
| `Type` | Type: `acceptor`, `donor`, `aromatic`, `hydrophobic`, `salt_bridge`, `pi_cation` |
| `Angle` | Validation angle in degrees |
| `Interaction` | `Yes` / `No` — whether distance and angle criteria are met |

### VMD scripts (if `vmd_output: 'Yes'`)

Three independent `.tcl` scripts are generated per pair, each self-contained (they load
the PDB copies already saved in the same output folder, so the folder can be moved or run
on a different machine without editing paths):

| File | Content |
|---|---|
| `vmd_<rec>_<lig>.tcl` | Full protein + active site residues (Licorice) + ligand (Licorice); dashed lines with distance labels for each validated H-bond/aromatic interaction — **white** aromatic, **red** ligand-acceptor, **yellow** ligand-donor |
| `vmd_hydrophobic_<rec>_<lig>.tcl` | Same base scene; dashed **orange** lines for each validated hydrophobic contact |
| `vmd_pockets_<rec>_<lig>.tcl` | Same base scene; one `Surf` (MSMS) representation per qualifying pocket (`Is_Pocket == Yes`), colour-rotated per pocket, covering that pocket's contacting residues |

The full protein is rendered with **Lines** (not `NewCartoon`/`Tube`/`Trace`): some VMD
builds — notably early `2.0.0` alpha releases — silently truncate spline-based backbone
representations to the first ~40 residues regardless of selection, a confirmed VMD bug
unrelated to the input PDB. `Lines` draws bond-by-bond and is unaffected, so it's used as
the reliable default; switch to `NewCartoon` manually in VMD's *Graphics > Representations*
if your VMD build renders it correctly.

`Surf`/`MSMS` in VMD doesn't expose a scriptable "Wireframe" draw style (only probe radius
and resolution are settable via `mol modstyle`) — to see the pocket surface as a mesh
instead of solid, change it manually: *Graphics > Representations* → select the pocket's
`Surf` rep → *Draw style* → Wireframe/Points.

### Ligand PNG images (if `ligand_plot: Yes`)

| File | Content |
|---|---|
| `<lig>_acceptors.png` | Ligand with acceptor atoms highlighted and labelled by PDB atom name |
| `<lig>_donors.png` | Ligand with donor atoms highlighted and labelled by PDB atom name |
| `<lig>_aromatic.png` | Ligand with aromatic rings highlighted; each ring has a distinct colour and is labelled R1, R2, … |

The atom names in the PNG images match the `Lig` column in the CSV files directly.

---

## Output folder structure

Everything is stored inside a single folder per pair `<receptor>_<ligand>/`:

```
<receptor>_<ligand>/
├── <receptor>.pdb             ← copy of the receptor PDB
├── <ligand>.pdb               ← copy of the ligand PDB
├── <ligand>_old.pdb           ← pre-cleanup copy (remove_bias)
├── Interaction_*_all.csv      ← all interactions, no filter
├── Interaction_*_threshold.csv← filtered by distance
├── Interaction_*_true.csv     ← validated by distance + angle
├── Pockets_*.csv              ← hydrophobic pocket candidates
├── summary.csv                ← interaction count by type
├── CM.csv                     ← ligand centre of mass
├── vmd_*.tcl                  ← main VMD script (H-bonds/aromatic) (if vmd_output: Yes)
├── vmd_hydrophobic_*.tcl      ← hydrophobic contacts VMD script (if vmd_output: Yes)
├── vmd_pockets_*.tcl          ← hydrophobic pockets VMD script (if vmd_output: Yes)
├── *_acceptors.png            ← ligand with acceptors highlighted
├── *_donors.png               ← ligand with donors highlighted
└── *_aromatic.png             ← ligand with aromatic rings highlighted
```

In batch mode each pair generates its own independent folder.

---

## Notes

- The script should be run from the directory containing the PDB files, or use absolute paths.
- For batch analysis of multiple ligands, the script can be called in a shell loop; with `cumulative_output: 'Yes'`, `Interactions_close.csv` and `CM_all.csv` are appended automatically across runs (set to `'No'` to disable).
- Non-standard residues not listed in `acceptors` / `donors` in the YAML are silently skipped.
- Aromatic rings in the ligand must contain more than 5 atoms and be planar (RMSD to the best-fit plane ≤ `Ring_Planarity_RMSD_Max`) to be considered. Planarity, not RDKit's aromaticity flag, is used because `Chem.MolFromPDBFile` does not reliably perceive aromaticity from PDB files without explicit bond orders.

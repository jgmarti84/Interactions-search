# Interactions Search

Detects and classifies non-covalent interactions between a ligand and a protein from PDB files. Identifies hydrogen bonds, aromatic interactions (π-π and T-shaped), hydrophobic contacts, salt bridges, and π-cation interactions. Each contact is validated by distance and angle criteria. Outputs CSV files and TCL scripts for visualization in VMD.

---

## Files

| File | Description |
|---|---|
| `Interactions_search.py` | Main script |
| `Interacciones_variables.yml` | Distance thresholds, acceptors and donors per residue |
| `Geometry.py` | Auxiliary geometry module (reference only, not imported directly) |

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

Each pair generates its own output folder. Cumulative CSVs (`Interactions_close.csv`, `CM_all.csv`) are appended automatically.

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
    └── Aromatic rings:    detected via ring_info, filtered by size > 5
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
[8] Outputs
    ├── CSV with all raw interactions
    ├── CSV filtered by distance
    ├── CSV filtered by distance + angle  (validated interactions)
    ├── Console summary  (table of validated interactions)
    ├── Cumulative summary Interactions_close.csv
    ├── Centre of mass    CM_all.csv
    └── TCL script for VMD (if vmd_output: Yes)
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
  ligand_plot: 'Yes'    # generates PNG images of ligand acceptors, donors and rings
  vmd_output:  'Yes'    # generates TCL script for VMD visualisation

distancias:
  Distances_Hidrogen_Bonds: 3.2   # Å — H-bond threshold
  Distances_Aromatic:       5.5   # Å — centre-to-centre aromatic threshold
  Distances_Hidrofobica:    4.0   # Å — hydrophobic threshold (aligned with PLIP)
  centroid_distance:        9.0   # Å — active site search radius (reference)
  Distances_C_Simple:       1.54  # Å — C-C single bond (reference)
  Distances_C_Doble:        2.56  # Å — C=C double bond (reference)

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
```

---

## Outputs

### CSV files

| File | Content |
|---|---|
| `<folder>/Interaction_<rec>_<lig>_all.csv` | All interactions found (no filters) |
| `<folder>/Interaction_<rec>_<lig>_threshold.csv` | Filtered by distance |
| `<folder>/Interaction_<rec>_<lig>_true.csv` | Validated by distance and angle |
| `Interactions_close.csv` | Cumulative run summary |
| `Interactions_all_count.csv` | Count by type (acceptor / donor / aromatic) |
| `CM_all.csv` | Ligand centre of mass per run |

Interaction CSV columns:

| Column | Description |
|---|---|
| `Pos R` | Receptor residue number |
| `Res` | Residue name (e.g. SER, TYR) |
| `Atom` | Receptor atom involved |
| `Dist` | Distance in Å |
| `Lig` | Ligand atom or ring involved |
| `Type` | Type: `acceptor`, `donor`, `aromatic`, `hydrophobic`, `salt_bridge`, `pi_cation` |
| `Angle` | Validation angle in degrees |
| `Interaction` | `Yes` / `No` — whether distance and angle criteria are met |

### VMD script (`vmd_<receptor>_<ligand>.tcl`)

Generates a visualisation ready to load in VMD:
- Full protein in transparent NewCartoon
- Active site residues in Licorice
- Ligand in Licorice
- Dashed lines for each validated interaction:
  - **White** — aromatic
  - **Red** — ligand acceptor (receptor donor)
  - **Yellow** — ligand donor (receptor acceptor)
- Distance label in Ångströms over each line

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
├── summary.csv                ← interaction count by type
├── CM.csv                     ← ligand centre of mass
├── vmd_*.tcl                  ← VMD script (if vmd_output: Yes)
├── *_acceptors.png            ← ligand with acceptors highlighted
├── *_donors.png               ← ligand with donors highlighted
└── *_aromatic.png             ← ligand with aromatic rings highlighted
```

In batch mode each pair generates its own independent folder.

---

## Notes

- The script should be run from the directory containing the PDB files, or use absolute paths.
- For batch analysis of multiple ligands, the script can be called in a shell loop; the cumulative CSVs (`Interactions_close.csv`, etc.) are appended automatically.
- Non-standard residues not listed in `acceptors` / `donors` in the YAML are silently skipped.
- Aromatic rings in the ligand must contain more than 5 atoms to be considered (filters out cyclopentane and similar non-aromatic rings).

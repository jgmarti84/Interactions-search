# Agentic Roadmap — interactions-search

Fine-grained roadmap for converting `Interactions_search.py` into a proper Python package.
Each step is designed for **one focused agent session** with a clear scope and testable acceptance criteria.

## Decisions in force

| Decision | Choice | Note |
|---|---|---|
| Product shape | CLI + importable library | Both first-class |
| Python baseline | 3.11+ | |
| Config format | YAML default, CLI overrides | Pydantic v2 validates |
| Yes/No option type | `Literal['Yes', 'No']` | Deferred to Phase 14 for `bool` coercion |
| Module layout | Nested subpackages | `io/`, `ligand/`, `receptor/`, `interactions/` |
| Function names | Freeze Spanish now, rename in Phase 12 | Smaller diffs per step |
| Language | English for new code and docs | Spanish preserved until rename step |

---

## How to use this roadmap

1. Open a new agent session.
2. Paste the **Agent prompt** block for that step verbatim.
3. The agent stops after its acceptance criteria are met.
4. Review the diff, then move to the next step.

Keep each step independent. If a step feels large, stop at the first acceptance criterion and open a new session.

---

## Current state (before Phase 0)

- `Interactions_search.py`: 1 214-line monolith, 35 functions, `main()` extracted
- `Interacciones_variables.yml`: config with `options`, `distancias`, `angulos`, `aromaticidad`, `acceptors`, `donors`, `acceptors_antecedent`, `special`
- `pyproject.toml`: skeleton — references non-existent `Geometry.py` module
- `src/interactions_search/`: directory exists but all `.py` sources are missing (only stale `.pyc` caches)
- `tests/`: test `.py` source file missing; only compiled `.pyc` cache remains
- `docs/BASELINE.md`: does not exist yet
- `data/PDK4/5m4m_ligands.pdb`: full complex PDB — suitable as integration fixture but too large for unit tests; no minimal fixture pair exists

---

## Phase 0 — Foundation cleanup

**Goal**: Establish a clean, working baseline before any refactoring begins. No logic changes.

### Step 0.1 — Fix pyproject.toml and add pydantic

**What to change**
- Remove `py-modules = ["Interactions_search", "Geometry"]` from `[tool.setuptools]`
- Add `packages = {find = {where = ["src"]}}` under `[tool.setuptools]`
- Add `"pydantic>=2.7"` to `[project].dependencies`
- Remove stale `__pycache__` directories under `src/` and `tests/`

**Acceptance criteria**
- `pip install -e .[dev]` completes without errors
- `python Interactions_search.py --help` still works

**Agent prompt**
```text
Fix pyproject.toml for the interactions-search project:
1. Remove the py-modules line that references Geometry.py (does not exist).
2. Add setuptools find_packages pointing to src/.
3. Add pydantic>=2.7 to project dependencies.
4. Delete any __pycache__ directories under src/ and tests/ (stale .pyc from missing sources).
5. Verify pip install -e .[dev] succeeds and python Interactions_search.py --help works.
Make no other changes.
```

---

### Step 0.2 — Minimal test fixtures

**What to create**
- `tests/fixtures/receptor_mini.pdb` — small receptor (≤ 50 ATOM lines, single chain A)
- `tests/fixtures/ligand_mini.pdb` — small ligand (≤ 15 HETATM lines)
- The pair must produce at least one interaction when run through the script

**Acceptance criteria**
- `python Interactions_search.py -r tests/fixtures/receptor_mini.pdb -l tests/fixtures/ligand_mini.pdb -c A` exits 0 and creates an output folder

**Agent prompt**
```text
Create minimal PDB fixture files for the interactions-search test suite.
- tests/fixtures/receptor_mini.pdb: a small protein fragment with at least one H-bond donor (e.g. backbone N of ALA) and one hydrophobic residue (e.g. LEU side chain), chain A.
- tests/fixtures/ligand_mini.pdb: a small HETATM ligand with at least one acceptor oxygen and one carbon-only region.
Keep both files under 60 lines. Verify that running:
  python Interactions_search.py -r tests/fixtures/receptor_mini.pdb -l tests/fixtures/ligand_mini.pdb -c A
exits successfully and produces an output folder. Document any limitations in comments inside each PDB file.
Do not modify Interactions_search.py.
```

---

### Step 0.3 — Smoke tests

**What to create**
- `tests/test_smoke.py` — at minimum 4 tests:
  1. Script exits 0 on fixture pair
  2. Output folder is created
  3. `_all.csv` has expected columns
  4. `_true.csv` rows all have `Interaction == 'Yes'`

**Acceptance criteria**
- `pytest tests/test_smoke.py` passes

**Agent prompt**
```text
Create tests/test_smoke.py for interactions-search.
Use tests/fixtures/receptor_mini.pdb and tests/fixtures/ligand_mini.pdb as input.
Test at minimum:
1. Running python Interactions_search.py -r ... -l ... -c A exits with code 0.
2. Output folder is created and contains Interaction_*_all.csv.
3. The CSV has the expected columns: ['Pos R', 'Res', 'Atom', 'Dist', 'Lig', 'Type', 'Angle', 'Interaction'].
4. All rows in Interaction_*_true.csv have Interaction == 'Yes'.
Use subprocess.run to invoke the script. Clean up output folders in a fixture teardown.
Run pytest to confirm all pass before finishing.
```

---

### Step 0.4 — Baseline documentation

**What to create**
- `docs/BASELINE.md` — documents exactly what the current script does and produces

**Acceptance criteria**
- File exists and accurately reflects current CLI interface, output files, and known behavior

**Agent prompt**
```text
Create docs/BASELINE.md for interactions-search.
Read Interactions_search.py and README.md. Document:
- All CLI arguments and their effects
- The three usage modes (-r/-l, batch, -x complex)
- The analysis pipeline step-by-step
- Output files produced per run
- Known limitations (e.g. Spanish function names, monolithic structure, no type safety)
- Current test coverage (what tests/test_smoke.py checks)
Keep it factual. This document is a snapshot for regression comparison during refactoring.
```

---

## Phase 1 — Pydantic config model

**Goal**: Replace the untyped 14-tuple returned by `carga_variables()` and the raw `cfg` dict with a validated Pydantic model. The rest of the script remains unchanged.

### Step 1.1 — Define Pydantic models

**File to create**: `src/interactions_search/config.py`

**Models**
```python
YesNo = Literal['Yes', 'No']   # flag: convert to bool in Phase 14

class Options(BaseModel):
    ligand_plot: YesNo
    vmd_output: YesNo
    cumulative_output: YesNo

class Distances(BaseModel):
    Distances_Hidrogen_Bonds: float   # gt=0
    Distances_Aromatic: float         # gt=0
    Distances_Hidrofobica: float      # gt=0
    centroid_distance: float          # gt=0
    Distances_C_Simple: float         # gt=0
    Distances_C_Doble: float          # gt=0

class Angles(BaseModel):
    Angle_Hidrogen_Bonds_Min: float   # ge=0, le=360
    Angle_Hidrogen_Bonds_Max: float   # ge=0, le=360

class Aromaticity(BaseModel):
    Ring_Planarity_RMSD_Max: float    # gt=0

class InteractionConfig(BaseModel):
    options: Options
    distancias: Distances
    angulos: Angles
    aromaticidad: Aromaticity
    acceptors: dict[str, list[str]]
    donors: dict[str, list[str]]
    acceptors_antecedent: dict[str, dict[str, str]]
    special: dict[str, list]
```

Also add `load_config(path: Path | None = None) -> InteractionConfig`.

**Acceptance criteria**
- `from interactions_search.config import load_config, InteractionConfig` works
- Loading `Interacciones_variables.yml` returns a valid `InteractionConfig`
- Setting `Distances_Hidrogen_Bonds: -1.0` raises `ValidationError`

**Agent prompt**
```text
Create src/interactions_search/config.py for interactions-search.
Implement Pydantic v2 models as described in docs/AGENTIC_ROADMAP.md Step 1.1.
Key constraints:
- YesNo = Literal['Yes', 'No'] — do NOT coerce to bool yet (flag: Phase 14)
- All distance fields must be gt=0
- Angle fields must be ge=0 and le=360
- load_config(path) reads the YAML at path (default: Interacciones_variables.yml next to the package root) and returns InteractionConfig
- ValidationError messages must be human-readable (use Field descriptions)
Also create tests/config/test_config_loading.py with at least:
1. test_valid_yaml_loads — loads Interacciones_variables.yml without error
2. test_negative_distance_raises — sets a distance < 0, expects ValidationError
3. test_invalid_yes_no_raises — sets ligand_plot: 'maybe', expects ValidationError
Run pytest before finishing.
```

---

### Step 1.2 — Wire config into carga_variables and analyze_pair

**What to change in `Interactions_search.py`**
- `carga_variables()` internally uses `load_config()`, then unpacks the model into the existing tuple return for backward compatibility
- `main()` builds `cfg` dict from the model's fields (names unchanged for now)
- `analyze_pair` signature unchanged; internally it still uses the `cfg` dict keys

**Acceptance criteria**
- All smoke tests still pass
- Invalid YAML raises a clear error before analysis starts

**Agent prompt**
```text
Wire src/interactions_search/config.py into Interactions_search.py.
Changes:
1. carga_variables() — replace the bare yaml.load call with load_config(); unpack into the same 14-element tuple it currently returns.
2. main() — no signature change; the cfg dict keys must remain identical.
Do NOT change analyze_pair or any other function.
Run tests/test_smoke.py after the change to confirm parity.
```

---

### Step 1.3 — CLI distance overrides

**New argparse flags** (add to `main()`)
```
--hbond-dist FLOAT      override Distances_Hidrogen_Bonds
--aromatic-dist FLOAT   override Distances_Aromatic
--hydrophobic-dist FLOAT override Distances_Hidrofobica
```

Precedence: CLI flag > YAML value.

**Acceptance criteria**
- `--hbond-dist 2.5` sets `cfg['Distances_Hidrogen_Bonds']` to 2.5 regardless of YAML
- Invalid float value prints argparse error and exits non-zero

**Agent prompt**
```text
Add CLI distance override flags to Interactions_search.py main():
  --hbond-dist, --aromatic-dist, --hydrophobic-dist
Each is optional float. When provided, override the corresponding key in cfg after carga_variables() loads the YAML.
Add tests/config/test_cli_overrides.py: invoke the script with --hbond-dist 2.5, capture cfg value via a monkeypatched analyze_pair, assert the override was applied.
Run all tests before finishing.
```

---

### Step 1.4 — Config docs

**File to create**: `docs/CONFIG.md`

**Acceptance criteria**
- Covers all YAML keys, types, valid ranges, and CLI override flags
- Notes the `YesNo` flag for future bool coercion

**Agent prompt**
```text
Create docs/CONFIG.md for interactions-search.
Document:
- Every YAML key in Interacciones_variables.yml with its type, valid range, and default value
- CLI override flags from Step 1.3 and their precedence rule
- The YesNo Literal['Yes','No'] convention and the note that Phase 14 will coerce these to bool
- How to use a custom config file (if load_config supports a path argument)
Keep it concise.
```

---

## Phase 2 — Package skeleton

**Goal**: Create the nested subpackage directory structure with stub `__init__.py` files. No logic moves yet.

### Step 2.1 — Create directory structure

**Directories to create**
```
src/interactions_search/
├── __init__.py
├── cli.py             (stub — from Interactions_search import main)
├── config.py          (already exists from Phase 1)
├── geometry.py        (stub)
├── io/
│   ├── __init__.py
│   ├── pdb.py         (stub)
│   └── output.py      (stub)
├── ligand/
│   ├── __init__.py
│   ├── hotpoints.py   (stub)
│   └── rings.py       (stub)
├── receptor/
│   ├── __init__.py
│   ├── active_site.py (stub)
│   └── points.py      (stub)
├── interactions/
│   ├── __init__.py
│   ├── hbonds.py      (stub)
│   ├── hydrophobic.py (stub)
│   ├── aromatic.py    (stub)
│   ├── salt_bridges.py (stub)
│   └── pi_cation.py   (stub)
├── vmd.py             (stub)
└── pipeline.py        (stub)
```

**Acceptance criteria**
- `pip install -e .` succeeds
- `python -m interactions_search.cli --help` works (delegates to Interactions_search.main)
- `pytest tests/test_smoke.py` still passes

**Agent prompt**
```text
Create the nested subpackage skeleton for src/interactions_search/ as described in docs/AGENTIC_ROADMAP.md Phase 2.
Each file is a stub with one-line docstrings only — no logic yet.
cli.py must do: from Interactions_search import main (so --help works via python -m interactions_search.cli).
Update pyproject.toml so find_packages discovers all subpackages under src/.
Verify pip install -e ., python -m interactions_search.cli --help, and pytest tests/test_smoke.py all pass.
```

---

## Phase 3 — Geometry module

**Goal**: Move pure geometry functions to `src/interactions_search/geometry.py`.

**Functions to move**
| Original name | Keep name? |
|---|---|
| `get_atom_coords` | yes |
| `center_of_mass` | yes |
| `angle_three_points` | yes |
| `center_aromatic_ring` | yes |
| `_ring_planarity_rmsd` | yes |

**Acceptance criteria**
- Functions importable from `interactions_search.geometry`
- `Interactions_search.py` imports them from there (no duplicate definitions)
- Unit tests for `angle_three_points` and `center_of_mass` pass

**Agent prompt**
```text
Move the five geometry functions from Interactions_search.py to src/interactions_search/geometry.py.
Functions: get_atom_coords, center_of_mass, angle_three_points, center_aromatic_ring, _ring_planarity_rmsd.
In Interactions_search.py replace each definition with: from interactions_search.geometry import <name>.
Create tests/geometry/test_geometry.py with deterministic unit tests:
- angle_three_points: right angle (90°), collinear (180°)
- center_of_mass: two atoms equidistant from origin → midpoint
Do NOT rename functions (Phase 12). Run all tests.
```

---

## Phase 4 — I/O layer

**Goal**: Move all PDB file I/O and output functions to dedicated modules.

### Step 4.1 — PDB I/O module (`io/pdb.py`)

**Functions to move**
- `extract_coords_from_pdb`
- `split_pdb`
- `remove_bias`
- `validate_inputs`

**Acceptance criteria**
- `from interactions_search.io.pdb import split_pdb` works
- Complex-mode smoke test (`-x`) still passes

**Agent prompt**
```text
Move extract_coords_from_pdb, split_pdb, remove_bias, validate_inputs from Interactions_search.py to src/interactions_search/io/pdb.py.
In Interactions_search.py replace definitions with imports from interactions_search.io.pdb.
Run pytest tests/test_smoke.py to confirm parity.
Do not rename functions.
```

---

### Step 4.2 — Output module (`io/output.py`)

**Functions to move**
- `_append_cumulative_csv`
- `print_summary`

**Acceptance criteria**
- `from interactions_search.io.output import print_summary` works
- Cumulative output smoke test (if any) still passes

**Agent prompt**
```text
Move _append_cumulative_csv and print_summary from Interactions_search.py to src/interactions_search/io/output.py.
In Interactions_search.py replace definitions with imports.
Run all smoke tests.
```

---

## Phase 5 — Ligand analysis modules

### Step 5.1 — Hotpoints module (`ligand/hotpoints.py`)

**Functions to move**
- `search_hot_points`
- `_draw_mol_labeled`
- `generate_df_ligand`
- `Busqueda_Antecesor_Lig`

**Acceptance criteria**
- `from interactions_search.ligand.hotpoints import search_hot_points` works
- Smoke test passes (ligand PNG generation still works if `ligand_plot: Yes`)

**Agent prompt**
```text
Move search_hot_points, _draw_mol_labeled, generate_df_ligand, Busqueda_Antecesor_Lig from Interactions_search.py to src/interactions_search/ligand/hotpoints.py.
Replace definitions with imports in Interactions_search.py.
Add tests/ligand/test_hotpoints.py: use a known SMARTS pattern (e.g. [O;H0]) on a simple RDKit mol, confirm the atom index is returned.
Run all tests.
```

---

### Step 5.2 — Rings module (`ligand/rings.py`)

**Functions to move**
- `search_rings`
- `visualize_rings`

**Acceptance criteria**
- `from interactions_search.ligand.rings import search_rings` works
- Smoke test passes

**Agent prompt**
```text
Move search_rings and visualize_rings from Interactions_search.py to src/interactions_search/ligand/rings.py.
Replace with imports. Run all smoke tests.
```

---

## Phase 6 — Receptor analysis modules

### Step 6.1 — Active site module (`receptor/active_site.py`)

**Functions to move**
- `active_site_residues`

**Agent prompt**
```text
Move active_site_residues from Interactions_search.py to src/interactions_search/receptor/active_site.py.
Replace with import. Run all tests.
```

---

### Step 6.2 — Receptor points module (`receptor/points.py`)

**Functions to move**
- `Coordenadas_interes_receptor`
- `get_aromatic_coord`

**Agent prompt**
```text
Move Coordenadas_interes_receptor and get_aromatic_coord from Interactions_search.py to src/interactions_search/receptor/points.py.
Replace with imports. Run all tests.
```

---

## Phase 7 — Interaction modules

Move each interaction type to its own file. Each step follows the same pattern: move, import, run tests.

### Step 7.1 — H-bond module (`interactions/hbonds.py`)

The H-bond logic is currently inline inside `analyze_pair` (loops over receptor_points × ligand acceptors/donors). Extract it into a standalone function.

**New function signature**
```python
def search_hbonds(
    receptor_points: pd.DataFrame,
    acceptor_atoms: list[int],
    donor_atoms: list[int],
    pdb_coords: list,
    DF_Lig: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
```

**Acceptance criteria**
- `analyze_pair` calls `search_hbonds` instead of inline logic
- Smoke test CSV output is identical

**Agent prompt**
```text
Extract the H-bond search logic from analyze_pair in Interactions_search.py into a new function search_hbonds in src/interactions_search/interactions/hbonds.py.
The function signature must match what is described in docs/AGENTIC_ROADMAP.md Step 7.1.
Update analyze_pair to call it.
Add tests/interactions/test_hbonds.py: two atoms within threshold → returns row; two atoms outside → empty DataFrame.
Run all tests.
```

---

### Step 7.2 — Hydrophobic module (`interactions/hydrophobic.py`)

**Functions to move**
- `search_hydrophobic`
- `_collapse_same_residue_contacts`

**Agent prompt**
```text
Move search_hydrophobic and _collapse_same_residue_contacts from Interactions_search.py to src/interactions_search/interactions/hydrophobic.py.
Replace with imports. Run all tests.
```

---

### Step 7.3 — Aromatic module (`interactions/aromatic.py`)

**Functions to move**
- `Interaccion_Aromatica`
- `aromatic_angle`
- `residuos_contacto`

**Agent prompt**
```text
Move Interaccion_Aromatica, aromatic_angle, residuos_contacto from Interactions_search.py to src/interactions_search/interactions/aromatic.py.
Replace with imports. Run all tests.
```

---

### Step 7.4 — Salt bridges module (`interactions/salt_bridges.py`)

**Function to move**: `search_salt_bridges`

**Agent prompt**
```text
Move search_salt_bridges from Interactions_search.py to src/interactions_search/interactions/salt_bridges.py.
Replace with import. Run all tests.
```

---

### Step 7.5 — Pi-cation module (`interactions/pi_cation.py`)

**Function to move**: `search_pi_cation`

**Agent prompt**
```text
Move search_pi_cation from Interactions_search.py to src/interactions_search/interactions/pi_cation.py.
Replace with import. Run all tests.
```

---

## Phase 8 — VMD output module

### Step 8.1 — VMD module (`vmd.py`)

**Functions to move**
- `_vmd_write_interaction`
- `scripting_vmd`
- `scripting_vmd_hydrophobic`

**Acceptance criteria**
- Generated TCL files use `[file join [file dirname [info script]] …]` — no absolute `/tmp` paths
- Test verifies TCL content does not contain `/tmp`

**Agent prompt**
```text
Move _vmd_write_interaction, scripting_vmd, scripting_vmd_hydrophobic from Interactions_search.py to src/interactions_search/vmd.py.
Replace with imports. 
Add tests/test_vmd.py: run scripting_vmd on fixture data and assert the resulting .tcl file does not contain any '/tmp' substring.
Run all tests.
```

---

## Phase 9 — Pipeline orchestrator

**Goal**: `analyze_pair` is the main orchestration function. Move it to `pipeline.py` and change it to accept `InteractionConfig` instead of the raw `cfg` dict.

### Step 9.1 — Move analyze_pair

**What changes**
- `analyze_pair(receptor_pdb, ligand_input, chain, cfg: dict)` → moves to `pipeline.py`
- Signature changes to `analyze_pair(receptor_pdb, ligand_input, chain, config: InteractionConfig)`
- All `cfg['key']` accesses become `config.distancias.Distances_Hidrogen_Bonds` etc.
- `Interactions_search.py` imports it; `main()` passes `InteractionConfig` instead of the dict

**Acceptance criteria**
- Full integration test on fixture produces identical CSV output
- `python Interactions_search.py -r ... -l ... -c A` still works

**Agent prompt**
```text
Move analyze_pair from Interactions_search.py to src/interactions_search/pipeline.py.
Change its third argument from cfg: dict to config: InteractionConfig (from interactions_search.config).
Replace all cfg['key'] dict accesses with the corresponding Pydantic model attributes.
In Interactions_search.py, update main() to call load_config() and pass the InteractionConfig object directly.
Remove the old carga_variables() 14-tuple unpacking and the cfg dict construction.
Run all tests and confirm the fixture CSV output is byte-for-byte identical to pre-refactor output (save a copy before running).
```

---

### Step 9.2 — Thin Interactions_search.py

After Phase 9.1, `Interactions_search.py` should contain only:
```python
from interactions_search.cli import main
if __name__ == '__main__':
    main()
```

Everything else is now in the package.

**Agent prompt**
```text
After Phase 9.1, verify that Interactions_search.py can be reduced to just the from…import main and if __name__ block.
Move the argparse setup and main() function body to src/interactions_search/cli.py.
Confirm python Interactions_search.py --help and python -m interactions_search.cli --help both work.
Run all tests.
```

---

## Phase 10 — Public API surface

**Goal**: Define what users can import directly from `interactions_search`.

**Public exports** (minimum viable)
```python
from interactions_search import analyze_pair, load_config, InteractionConfig
```

**Files to create/update**
- `src/interactions_search/__init__.py` — selective exports
- `docs/API.md` — usage examples for both CLI and library

**Agent prompt**
```text
Define the public API for interactions_search.
In src/interactions_search/__init__.py, export: analyze_pair, load_config, InteractionConfig.
Create docs/API.md with:
- CLI usage examples for all three modes (-r/-l, batch, -x)
- Python import examples: load_config, analyze_pair
- Notes on what is and is not part of the public API
Run all tests.
```

---

## Phase 11 — Integration test expansion

**Goal**: Add integration tests that run the full pipeline on the fixture pair and compare CSV outputs to a saved snapshot.

**Agent prompt**
```text
Add tests/test_integration.py for interactions-search.
Run the full pipeline on tests/fixtures/receptor_mini.pdb + ligand_mini.pdb.
Save the output CSVs as golden files in tests/fixtures/golden/.
On subsequent runs, compare current output to golden files column-by-column (not byte-for-byte, to tolerate float rounding).
Add a --update-golden flag (via an environment variable INTERACTIONS_UPDATE_GOLDEN=1) that regenerates the golden files.
Run all tests.
```

---

## Phase 12 — Rename Spanish identifiers (isolated step)

**Goal**: Rename all Spanish function names, variable names, and comments to English. This is a single dedicated step with no logic changes.

**Rename map** (partial — agent must complete from full source)
| Old name | New name |
|---|---|
| `carga_variables` | `load_config_legacy` → removed (replaced by `load_config`) |
| `Coordenadas_interes_receptor` | `receptor_points_of_interest` |
| `Busqueda_Antecesor_Lig` | `find_ligand_antecedent` |
| `Interaccion_Aromatica` | `aromatic_interaction` |
| `residuos_contacto` | `contact_residues` |
| `Aceptores_Prot` / `Dadores_Prot` | `acceptors` / `donors` |
| `Distancia_Hidrofobica` | `hydrophobic_dist` |
| `Distancia_Centro_Activo` | `active_site_radius` |

**Acceptance criteria**
- No Spanish identifiers remain in `src/` (grep verifiable)
- All tests pass after rename

**Agent prompt**
```text
Rename all Spanish identifiers in the interactions-search package (src/ directory only).
For each function or variable with a Spanish name: rename it and update all call sites.
Do not change any logic, only names.
After renaming, run grep -r "[A-Za-z]*[áéíóúñÁÉÍÓÚÑ][A-Za-z]*" src/ to verify no accented characters remain.
Also grep for common Spanish words in identifiers (Distancia, Aceptor, Dador, Recept, Prot, Coordenadas, Busqueda, Antecesor, Residuo) and rename any remaining ones.
Run all tests.
```

---

## Phase 13 — Quality gates + CI

**Agent prompt**
```text
Add GitHub Actions CI for interactions-search.
Create .github/workflows/ci.yml:
- Trigger: push and pull_request on main and feature/*
- Matrix: Python 3.11 only (expand later)
- Steps: checkout, setup-python, pip install -e .[dev], ruff check src/ tests/, pytest tests/
Also fix any ruff violations in src/ before committing the workflow.
```

---

## Phase 14 — Yes/No → bool coercion (deferred)

**Context**: Options in `Interacciones_variables.yml` currently use `'Yes'`/`'No'` strings and are typed as `Literal['Yes', 'No']` in the Pydantic model. All downstream code uses `== 'Yes'` comparisons.

**What this step does**
- Change `YesNo = Literal['Yes', 'No']` to `bool` in `config.py`
- Add YAML validator: accept both `'Yes'`/`'No'` (for backward compat) and `true`/`false`
- Replace every `== 'Yes'` and `== 'No'` in `src/` with direct boolean checks
- Update `Interacciones_variables.yml` to use `true`/`false`
- Update `docs/CONFIG.md`

**Agent prompt**
```text
Convert the Yes/No string options in interactions-search to proper booleans.
In src/interactions_search/config.py:
- Change Options fields from Literal['Yes','No'] to bool
- Add a field_validator that maps 'Yes' → True, 'No' → False for backward compatibility
In src/ replace every == 'Yes' and == 'No' with direct bool checks.
Update Interacciones_variables.yml to use true/false instead of 'Yes'/'No'.
Update docs/CONFIG.md to reflect the change and note the backward-compat validator.
Run all tests.
```

---

## Guardrails (applies to every step)

- One step = one agent session = one PR
- Do not mix logic changes with refactoring in the same step
- Run tests after every change; a failing test is a blocker, not a warning
- Preserve scientific thresholds exactly — no default value changes unless explicitly requested
- Update docs together with code
- Spanish identifiers in the monolith are tolerated until Phase 12; do not rename them in passing

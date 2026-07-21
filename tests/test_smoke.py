"""Smoke tests — freeze current behavior of Interactions_search.py.

Each test uses a module-scoped fixture that runs the script once against the
minimal fixture PDBs and cleans up all output afterward.
"""
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent
PYTHON = REPO_ROOT / ".venv" / "bin" / "python"
SCRIPT = REPO_ROOT / "Interactions_search.py"
RECEPTOR = REPO_ROOT / "tests" / "fixtures" / "receptor_mini.pdb"
LIGAND = REPO_ROOT / "tests" / "fixtures" / "ligand_mini.pdb"

OUTPUT_FOLDER = REPO_ROOT / "receptor_mini_ligand_mini"
CUMULATIVE_FILES = [
    REPO_ROOT / "Interactions_close.csv",
    REPO_ROOT / "CM_all.csv",
]

EXPECTED_COLUMNS = ["Pos R", "Res", "Atom", "Dist", "Lig", "Type", "Angle", "Interaction"]


def _cleanup():
    if OUTPUT_FOLDER.exists():
        shutil.rmtree(OUTPUT_FOLDER)
    for f in CUMULATIVE_FILES:
        f.unlink(missing_ok=True)


@pytest.fixture(scope="module")
def run_output():
    _cleanup()
    result = subprocess.run(
        [
            str(PYTHON), str(SCRIPT),
            "-r", str(RECEPTOR),
            "-l", str(LIGAND),
            "-c", "A",
        ],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    yield result
    _cleanup()


def test_exit_code_zero(run_output):
    assert run_output.returncode == 0, (
        f"Script exited {run_output.returncode}\n--- stderr ---\n{run_output.stderr}"
    )


def test_output_folder_and_all_csv_exist(run_output):
    assert OUTPUT_FOLDER.is_dir(), "Output folder was not created"
    csv_files = list(OUTPUT_FOLDER.glob("Interaction_*_all.csv"))
    assert csv_files, "No Interaction_*_all.csv found in output folder"


def test_all_csv_has_expected_columns(run_output):
    csv = next(OUTPUT_FOLDER.glob("Interaction_*_all.csv"))
    df = pd.read_csv(csv, index_col=0)
    assert list(df.columns) == EXPECTED_COLUMNS, (
        f"Column mismatch.\n  got : {list(df.columns)}\n  want: {EXPECTED_COLUMNS}"
    )


def test_true_csv_interactions_all_yes(run_output):
    true_csv = next(OUTPUT_FOLDER.glob("Interaction_*_true.csv"))
    df = pd.read_csv(true_csv, index_col=0)
    assert not df.empty, "_true.csv is empty — expected at least one validated interaction"
    bad = df[df["Interaction"] != "Yes"]
    assert bad.empty, (
        f"{len(bad)} row(s) in _true.csv have Interaction != 'Yes':\n{bad}"
    )

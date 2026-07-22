import numpy as np
import pytest
from Bio.PDB.Atom import Atom

from interactions_search.geometry import angle_three_points, center_of_mass, center_aromatic_ring, _ring_planarity_rmsd


def _pdb_coords(atoms):
    """Build a minimal pdb_coords dict from {key: (x, y, z)} mapping.
    Slots 0-4 are unused by _ring_planarity_rmsd; x/y/z go to indices 5/6/7.
    """
    return {key: (None, None, None, None, None, x, y, z) for key, (x, y, z) in atoms.items()}


def test_angle_right_angle():
    # Vectors Aceptor→Donor and Aceptor→Antecedent are perpendicular → 90°
    result = angle_three_points([1, 0, 0], [0, 0, 0], [0, 1, 0])
    assert abs(result - 90.0) < 1e-10


def test_angle_collinear():
    # Donor, Aceptor, Antecedent are collinear and on opposite sides → 180°
    result = angle_three_points([-1, 0, 0], [0, 0, 0], [1, 0, 0])
    assert abs(result - 180.0) < 1e-10


def test_center_of_mass_geometric_midpoint():
    # Two atoms at (-1,0,0) and (1,0,0); geometric center must be the origin.
    atom1 = Atom('N', np.array([-1.0, 0.0, 0.0]), 1.0, 1.0, ' ', ' N  ', 1, 'N')
    atom2 = Atom('C', np.array([1.0, 0.0, 0.0]), 1.0, 1.0, ' ', ' C  ', 2, 'C')
    result = center_of_mass([atom1, atom2], geometric=True)
    assert np.allclose(result, [0.0, 0.0, 0.0])


def test_center_aromatic_ring_regular_hexagon():
    """PHE/TYR atom ordering: CG=0, CD1=1, CD2=2, CE1=3, CE2=4, CZ=5."""
    ring = [
        [ 0.0,  1.4, 0.0],   # CG
        [-1.2,  0.7, 0.0],   # CD1
        [ 1.2,  0.7, 0.0],   # CD2
        [-1.2, -0.7, 0.0],   # CE1
        [ 1.2, -0.7, 0.0],   # CE2
        [ 0.0, -1.4, 0.0],   # CZ
    ]
    cx, cy, cz = center_aromatic_ring(ring)
    assert abs(cx) < 1e-6
    assert abs(cy) < 1e-6
    assert abs(cz) < 1e-9


# --- _ring_planarity_rmsd ---

def test_planarity_rmsd_perfect_planar_ring():
    """Six atoms on a regular hexagon in the XY plane → RMSD must be (near) zero."""
    r = 1.4  # benzene C-C bond length approximation
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    atoms = {i: (r * np.cos(a), r * np.sin(a), 0.0) for i, a in enumerate(angles)}
    coords = _pdb_coords(atoms)
    ring = list(atoms.keys())
    assert _ring_planarity_rmsd(coords, ring) < 1e-10


def test_planarity_rmsd_nonplanar_chair_ring():
    """Chair cyclohexane: atoms alternate ±0.25 Å out of the mean plane → RMSD = 0.25 Å."""
    r = 1.25
    d = 0.25  # half-amplitude of chair puckering
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    atoms = {i: (r * np.cos(a), r * np.sin(a), d * ((-1) ** i)) for i, a in enumerate(angles)}
    coords = _pdb_coords(atoms)
    ring = list(atoms.keys())
    assert _ring_planarity_rmsd(coords, ring) == pytest.approx(d, abs=1e-10)


def test_planarity_rmsd_triangle_always_zero():
    """Any three non-collinear points span a plane exactly → RMSD must be zero."""
    atoms = {0: (0.0, 0.0, 0.0), 1: (1.0, 0.0, 0.5), 2: (0.0, 1.0, -0.3)}
    coords = _pdb_coords(atoms)
    ring = [0, 1, 2]
    assert _ring_planarity_rmsd(coords, ring) < 1e-10
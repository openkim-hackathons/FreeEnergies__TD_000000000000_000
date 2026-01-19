"""Utility functions for structure manipulation and supercell generation."""

from typing import Tuple

import numpy as np


def compute_supercell_reps_for_cutoff(cell: np.ndarray, r: float) -> Tuple[int, int, int]:
    """Compute supercell repetitions to ensure minimum image distance >= 2*r.

    This creates a supercell where each lattice direction has enough repetitions
    to ensure the minimum image distance is at least 2*r (i.e., the supercell
    can contain a sphere of radius r without self-interaction across periodic
    boundaries). This is particularly useful for non-cubic cells where uniform
    expansion would be suboptimal.

    Args:
        cell: 3x3 lattice vectors of the primitive cell (Å)
        r: Target radius in Å

    Returns:
        Tuple of repetitions (n_a, n_b, n_c) along each lattice vector
    """
    cell = np.asarray(cell)
    a, b, c = cell[0], cell[1], cell[2]

    # Volume of parallelepiped: v_p = |a · (b × c)|
    v_p = abs(np.dot(a, np.cross(b, c)))

    # Heights (perpendicular distances from origin to opposite face)
    # h_a = v_p / |b × c| (height along a direction)
    h_a = v_p / np.linalg.norm(np.cross(b, c))
    h_b = v_p / np.linalg.norm(np.cross(a, c))
    h_c = v_p / np.linalg.norm(np.cross(a, b))

    # Number of repeats needed: n = ceil(2*r / h)
    n_a = max(1, int(np.ceil(2 * r / h_a)))
    n_b = max(1, int(np.ceil(2 * r / h_b)))
    n_c = max(1, int(np.ceil(2 * r / h_c)))

    return (n_a, n_b, n_c)


def compute_supercell_reps_uniform_cubic(
    n_atoms_primitive: int, target_atoms: int
) -> Tuple[int, int, int]:
    """Compute uniform cubic supercell repetitions to reach a target atom count.

    Calculates the smallest integer `n` such that a supercell of `(n, n, n)`
    repetitions contains at least `target_atoms`. This method assumes a cubic
    expansion.

    Args:
        n_atoms_primitive: Number of atoms in the primitive unit cell.
        target_atoms: The desired minimum number of atoms in the supercell.

    Returns:
        Tuple of repetitions (n, n, n) for uniform cubic expansion.
    """
    n = max(1, int(np.ceil(np.cbrt(target_atoms / n_atoms_primitive))))
    return (n, n, n)

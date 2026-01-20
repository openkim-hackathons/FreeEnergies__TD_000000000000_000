"""Utility functions for structure manipulation and supercell generation."""

from typing import Tuple

import numpy as np
from ase import Atoms


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


def compute_supercell_for_target_size(
    atoms: Atoms,
    target_size: int = 10000,
    initial_radius: float = 20.0,
    radius_step: float = 0.5,
    min_radius: float = 1.0
) -> Tuple[Atoms, Tuple[int, int, int]]:
    """Compute supercell using cutoff-based approach with target size constraint.

    Creates a supercell by iteratively adjusting the cutoff radius until the
    resulting supercell has fewer atoms than the target size. Starts with
    initial_radius and decreases it recursively until the target is met.

    This approach is particularly useful for non-cubic cells where uniform
    expansion would be suboptimal, but you still want to control the total
    number of atoms.

    Args:
        atoms: Primitive unit cell as ASE Atoms object.
        target_size: Maximum number of atoms in the supercell. Default: 10000.
        initial_radius: Starting cutoff radius in Å. Default: 8.0.
        radius_step: Amount to decrease radius by in each iteration (Å). Default: 0.5.
        min_radius: Minimum radius to try before giving up (Å). Default: 1.0.

    Returns:
        Tuple of (supercell Atoms object, repeat tuple (n_a, n_b, n_c)).
        The supercell has natoms < target_size.

    Raises:
        ValueError: If unable to find a supercell below target_size even at min_radius.
    """
    # Compute repetitions for current radius
    repeat = compute_supercell_reps_for_cutoff(atoms.get_cell(), initial_radius)
    
    # Create supercell
    supercell = atoms.repeat(repeat)
    natoms = len(supercell)
    
    # Base case: if we're below target size, return the supercell and repeat
    if natoms < target_size:
        return supercell, repeat
    
    # If we've hit the minimum radius, raise an error
    if initial_radius <= min_radius:
        raise ValueError(
            f"Unable to find supercell below target_size={target_size} atoms. "
            f"Even at min_radius={min_radius}Å, got {natoms} atoms."
        )
    
    # Recursive case: decrease radius and try again
    new_radius = max(min_radius, initial_radius - radius_step)
    return compute_supercell_for_target_size(
        atoms, target_size, new_radius, radius_step, min_radius
    )

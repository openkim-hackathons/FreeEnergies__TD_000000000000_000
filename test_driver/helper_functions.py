import numpy as np
from ase import Atoms
from typing import Dict, Iterable, List, Tuple
from ase import Atoms

def reduce_and_avg(atoms: Atoms, repeat: Tuple[int, int, int]) -> Atoms:
    """
    Function to reduce all atoms to the original unit cell position.
    """
    new_atoms = atoms.copy()

    cell = new_atoms.get_cell()

    # Divide each unit vector by its number of repeats.
    # See https://stackoverflow.com/questions/19602187/numpy-divide-each-row-by-a-vector-element.
    cell = cell / np.array(repeat)[:, None]

    # Decrease size of cell in the atoms object.
    new_atoms.set_cell(cell)
    new_atoms.set_pbc((True, True, True))

    # Set averaging factor
    M = np.prod(repeat)

    # Wrap back the repeated atoms on top of the reference atoms in the original unit cell.
    positions = new_atoms.get_positions(wrap=True)

    number_atoms = len(new_atoms)
    original_number_atoms = number_atoms // M
    assert number_atoms == original_number_atoms * M
    positions_in_prim_cell = np.zeros((original_number_atoms, 3))

    # Start from end of the atoms because we will remove all atoms except the reference ones.
    for i in reversed(range(number_atoms)):
        if i >= original_number_atoms:
            # Get the distance to the reference atom in the original unit cell with the
            # minimum image convention.
            distance = new_atoms.get_distance(i % original_number_atoms, i,
                                              mic=True, vector=True)
            # Get the position that has the closest distance to the reference atom in the
            # original unit cell.
            position_i = positions[i % original_number_atoms] + distance
            # Remove atom from atoms object.
            new_atoms.pop()
        else:
            # Atom was part of the original unit cell.
            position_i = positions[i]
        # Average.
        positions_in_prim_cell[i % original_number_atoms] += position_i / M

    new_atoms.set_positions(positions_in_prim_cell)

    return new_atoms


from typing import Dict, Iterable, List, Optional, Tuple
from ase import Atoms
from ase.geometry import get_distances

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from scipy.stats import kstest
from sklearn.decomposition import PCA
from kim_tools import KIMTestDriverError


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
    avg_positions_in_prim_cell = np.zeros((original_number_atoms, 3))
    positions_in_prim_cell = np.zeros((number_atoms, 3))

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
        # Average
        avg_positions_in_prim_cell[i % original_number_atoms] += position_i / M
        positions_in_prim_cell[i] = position_i

    new_atoms.set_positions(avg_positions_in_prim_cell)

    # Calculate the distances.
    distances = np.empty((original_number_atoms, M, 3))
    for i in range(number_atoms):
        dr, _ = get_distances(positions_in_prim_cell[i], avg_positions_in_prim_cell[i % original_number_atoms],
                              cell=new_atoms.get_cell(), pbc=True)
        # dr is a distance matrix, here we only have one distance
        assert dr.shape == (1, 1, 3)
        distances[i % original_number_atoms, i // original_number_atoms] = dr[0][0]

    return new_atoms, distances


def test_reduced_distances(reduced_distances: npt.NDArray[float], significance_level: float = 0.05,
                           plot_filename: Optional[str] = None, number_bins: Optional[int] = None ) -> None:
    """Function to test whether the reduced atom positions are normally distributed around their average."""
    assert len(reduced_distances.shape) == 3
    assert reduced_distances.shape[2] == 3

    if plot_filename is not None:
        if number_bins is None:
            raise ValueError("number_bins must be specified if plot_filename is specified")
        if not plot_filename.endswith(".pdf"):
            raise ValueError(f"{plot_filename} is not a PDF file")
        with PdfPages(plot_filename) as pdf:
            for i in range(reduced_distances.shape[0]):
                fig, axs = plt.subplots(1, 3, figsize=(10.0, 4.0))
                for j in range(reduced_distances.shape[2]):
                    axs[j].hist(reduced_distances[i, :, j], bins=number_bins)
                    axs[j].set_xlabel(f"$x_{j}$")
                axs[0].set_ylabel(f"Counts")
                fig.suptitle(f"Atom {i}")
                pdf.savefig()
    else:
        if number_bins is not None:
            raise ValueError("number_bins must not be specified if plot_filename is not specified")

    p_values = np.empty((reduced_distances.shape[0], reduced_distances.shape[2]))
    for i in range(reduced_distances.shape[0]):
        atom_distances = reduced_distances[i]

        # Perform PCA on the xyz distribution.
        pca = PCA(n_components=atom_distances.shape[1])
        pca_components = pca.fit_transform(atom_distances)
        assert pca_components.shape == atom_distances.shape == reduced_distances.shape[1:]

        # Test each component with a KS test.
        for j in range(pca_components.shape[1]):
            component = pca_components[:, j]
            component_mean = np.mean(component)
            assert abs(component_mean) < 1.0e-7
            component_std = np.std(component)
            # Normalize component
            normalized_component = (component - component_mean) / component_std
            assert abs(np.mean(normalized_component)) < 1.0e-7
            assert abs(np.std(normalized_component) - 1.0) < 1.0e-7
            res = kstest(normalized_component, "norm")
            p_values[i, j] = res.pvalue

    if np.any(p_values <= significance_level):
        raise KIMTestDriverError(f"Detected non-normal distribution of reduced atom positions around their average (smallest p value {np.min(p_values)}).")
    else:
        print(f"Detected normal distribution or reduced atom positions around their average (smallest p value {np.min(p_values)}).")








# -*- coding: utf-8 -*-
"""
This module contains the DecayDatabase, DecayMatrix, MatrixBuilder, and DecayDiagram classes.

DecayDatabase:
    - Represents a database of nuclear decay data for a set of nuclides.
    - Provides access to half-lives, decay modes, progeny, branching ratios, decay energies, and related utilities.

DecayMatrix:
    - Stores the decay matrix and related structures for nuclide decay calculations.

MatrixBuilder:
    - Builds the decay matrix and related structures for nuclide decay calculations.

DecayDiagram:
    - Visualizes nuclear decay chains and their reverse (parent) chains for a given nuclide.
    - Provides methods to plot decay diagrams and nuclear charts using networkx and matplotlib.
"""

from typing import Any, Dict, List, Tuple, Union, Optional
from collections import deque

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import sparse

from .nuclide import Nuclide, NuclideStrError
from .utils import HalfLifeColorMap


#---------------------------------
#    class DecayDatabase
#---------------------------------
class DecayDatabase:
    """
    Represents a database of nuclear decay data for a set of nuclides.
    Provides access to half-lives, decay modes, progeny, branching ratios, decay energies, and related utilities.
    """

    def __init__(
            self,
            data_source: str,
            nuclides: np.ndarray,
            half_life_data: np.ndarray,
            decay_constants_data: np.ndarray,
            decay_modes_data: np.ndarray,
            progeny_data: np.ndarray,
            branching_ratios_data: np.ndarray,
            decay_energies_data: np.ndarray
        )-> None:
        """
        Initialize the DecayDatabase with all relevant nuclear decay data arrays.
        Args:
            data_source (str): Name of the data source (e.g., ENDF-B-VIII.1_decay).
            nuclides (np.ndarray): Array of nuclide symbols.
            half_life_data (np.ndarray): Array of half-life data tuples.
            decay_constants_data (np.ndarray): Array of decay constants.
            decay_modes_data (np.ndarray): Array of decay modes for each nuclide.
            progeny_data (np.ndarray): Array of progeny lists for each nuclide.
            branching_ratios_data (np.ndarray): Array of branching ratios for each nuclide.
            decay_energies_data (np.ndarray): Array of decay energies for each nuclide.
        """
        self.data_source = data_source
        self.nuclides = nuclides
        self.half_life_data = half_life_data
        self.decay_constants_data = decay_constants_data
        self.decay_modes_data = decay_modes_data
        self.progeny_data = progeny_data
        self.branching_ratios_data = branching_ratios_data
        self.decay_energies_data = decay_energies_data

        # Map from nuclide symbol to its index in the arrays
        self.nuclide_index_map: Dict[str, int] = {
            nuclide: index for index, nuclide in enumerate(self.nuclides)
        }


    def validate_nuclide_symbol(self, nuclide: str) -> str:
        """
        Validates the nuclide string and returns its standardized symbol.
        Raises:
            NuclideStrError: If the nuclide string is invalid or not found in the database.
        Returns:
            str: Standardized nuclide symbol.
        """
        try:
            symbol = Nuclide(nuclide).nuclide_symbol
        except Exception as e:
            raise NuclideStrError(nuclide, f"Invalid format: {str(e)}")
        if symbol not in self.nuclide_index_map:
            raise NuclideStrError(nuclide, f"Not found in {self.data_source} decay database.")
        return symbol
    

    def get_nuclide_half_life(self, nuclide: str, units: str = 'readable') -> Union[float, str]:
        """
        Retrieves the half-life of a nuclide in the specified units.
        Args:
            nuclide (str): Nuclide symbol or string.
            units (str): 'readable', 's', or other (returns value in original units).
        Returns:
            float or str: Half-life in the requested units or as a readable string.
        """
        nuclide = self.validate_nuclide_symbol(nuclide)
        (
            half_life_second,
            half_life,
            unit,
            half_life_readable
        ) = self.half_life_data[self.nuclide_index_map[nuclide]]

        if units == 'readable':
            return half_life_readable
        if units == 's':
            return half_life_second

        return half_life


    def get_nuclide_progeny(self, nuclide) -> List[str]:
        """
        Returns the list of progeny (daughter nuclides) for a given nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            list[str]: List of progeny nuclide symbols.
        """
        nuclide = self.validate_nuclide_symbol(nuclide)
        return self.progeny_data[self.nuclide_index_map[nuclide]]
    

    def get_nuclide_ancestor(self, nuclide: str) -> List[str]:
        """
        Returns the list of ancestor nuclides for a given nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            list[str]: List of parent nuclide symbols.
        """
        nuclide = self.validate_nuclide_symbol(nuclide)
        ancestors = []
        for parent in self.nuclides:
            progeny = self.get_nuclide_progeny(parent)
            if nuclide in progeny:
                ancestors.append(parent)

        return ancestors
    
    def get_decay_channel(self, parent: str, daughter: str) -> Optional[Tuple[str, float]]:
        """
        Get decay channel (mode and branching ratio) between parent and daughter nuclides.
        Args:
            parent (str): Parent nuclide symbol.
            daughter (str): Daughter nuclide symbol.
        Returns:
            Optional[Tuple[str, float]]: (decay_mode, branching_ratio) if found, None otherwise.
        """
        parent = self.validate_nuclide_symbol(parent)
        daughter = self.validate_nuclide_symbol(daughter)
        progeny = self.get_nuclide_progeny(parent)

        if parent not in self.nuclide_index_map or daughter not in self.nuclide_index_map:
            raise NuclideStrError(f'Parent {parent} or Daughter {daughter}', f"Not found in {self.data_source} decay database.")
        if daughter not in progeny:
            raise ValueError(f"Daughter {daughter} not in progeny of {parent}.")

        daughter_index = progeny.index(daughter)
        decay_mode = self.get_nuclide_decay_modes(parent)[daughter_index]
        branching_ratio = self.get_nuclide_branching_ratios(parent)[daughter_index]

        return decay_mode, branching_ratio
        
        
    def get_nuclide_branching_ratios(self, nuclide: str) -> List[float]:
        """
        Returns the branching ratios for all decay modes of a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            list[float]: Branching ratios for each decay mode.
        """
        nuclide = self.validate_nuclide_symbol(nuclide)
        return self.branching_ratios_data[self.nuclide_index_map[nuclide]]
    
    
    def get_nuclide_decay_modes(self, nuclide: str) -> List[str]:
        """
        Returns the decay modes for a given nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            list[str]: Decay modes for the nuclide.
        """
        nuclide = self.validate_nuclide_symbol(nuclide)
        return self.decay_modes_data[self.nuclide_index_map[nuclide]]
    

    def get_nuclide_decay_constants(self, nuclide: str) -> np.ndarray:
        """
        Returns the decay constant (lambda) for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Decay constant in 1/s.
        """
        nuclide = self.validate_nuclide_symbol(nuclide)
        return self.decay_constants_data[self.nuclide_index_map[nuclide]]
    

    def get_proton_numbers(self) -> np.ndarray:
        """
        Returns an array of proton numbers (Z) for all nuclides in the database.
        Returns:
            np.ndarray: Array of proton numbers.
        """
        return np.array([Nuclide(nuclide).Z for nuclide in self.nuclides])


    def get_neutron_numbers(self) -> np.ndarray:
        """
        Returns an array of neutron numbers (N) for all nuclides in the database.
        Returns:
            np.ndarray: Array of neutron numbers.
        """
        return np.array([Nuclide(nuclide).N for nuclide in self.nuclides])
    

    def get_nuclide_decay_energy(self, nuclide: str, energy_type: str) -> float:
        """
        Returns a specific type of decay energy for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
            energy_type (str): Type of decay energy (e.g., 'EM', 'LP', 'HP', 'Neutrino', etc.).
        Returns:
            float: Decay energy value.
        Raises:
            ValueError: If the energy_type is invalid.
        """
        energy_map = {
            'EM':          0,
            'LP':          1,
            'HP':          2,
            'Neutrino':    3,
            'Gamma':       4,
            'Beta_Minus':  5,
            'Beta_Plus':   6,
            'Alpha':       7,
            'Neutron':     8,
            'Proton':      9,
            'Effective_Q': 10
        }
        if energy_type not in energy_map:
            raise ValueError(f"Invalid energy type: {energy_type}. Use one of {list(energy_map.keys())}.")
        nuclide = self.validate_nuclide_symbol(nuclide)

        return self.decay_energies_data[self.nuclide_index_map[nuclide]][energy_map[energy_type]]


    def get_nuclide_electromagnetic_energy(self, nuclide: str) -> float:
        """
        Returns the electromagnetic decay energy for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Electromagnetic decay energy.
        """
        return self.get_nuclide_decay_energy(nuclide, 'EM')


    def get_nuclide_light_particle_energy(self, nuclide: str) -> float:
        """
        Returns the light particle decay energy for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Light particle decay energy.
        """
        return self.get_nuclide_decay_energy(nuclide, 'LP')


    def get_nuclide_heavy_particle_energy(self, nuclide: str) -> float:
        """
        Returns the heavy particle decay energy for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Heavy particle decay energy.
        """
        return self.get_nuclide_decay_energy(nuclide, 'HP')
    

    def get_nuclide_neutrino_energy(self, nuclide: str) -> float:
        """
        Returns the neutrino decay energy for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Neutrino decay energy.
        """
        return self.get_nuclide_decay_energy(nuclide, 'Neutrino')
    
    
    def get_nuclide_total_decay_energy(self, nuclide: str) -> float:
        """
        Returns the total decay energy released (sum of EM, LP, HP, and Neutrino energies).
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Total decay energy released.
        """
        energy = np.array([
            self.get_nuclide_decay_energy(nuclide, 'EM'),
            self.get_nuclide_decay_energy(nuclide, 'LP'),
            self.get_nuclide_decay_energy(nuclide, 'HP'),
            self.get_nuclide_decay_energy(nuclide, 'Neutrino')
        ])
        return np.nansum(energy)


    def __eq__(self, other: object) -> bool:
        """
        Checks equality for DecayDatabase objects.
        Returns:
            bool: True if all relevant data arrays and source match, False otherwise.
        """
        if not isinstance(other, DecayDatabase):
            return NotImplemented
        return (
            self.data_source == other.data_source
            and np.array_equal(self.nuclides, other.nuclides)
            and np.array_equal(self.half_life_data, other.half_life_data)
            and np.array_equal(self.decay_modes_data, other.decay_modes_data)
            and np.array_equal(self.progeny_data, other.progeny_data)
            and np.array_equal(self.branching_ratios_data, other.branching_ratios_data)
            and self.nuclide_index_map == other.nuclide_index_map
        )


    def __ne__(self, other: object) -> bool:
        """
        Checks inequality for DecayDatabase objects.
        Returns:
            bool: True if not equal, False otherwise.
        """
        if not isinstance(other, DecayDatabase):
            return NotImplemented
        return not self.__eq__(other)


#---------------------------------
#    class DecayMatrix
#---------------------------------
class DecayMatrix:
    """
    Stores the decay matrix and related structures for nuclide decay calculations.
    Contains decay constants, transformation matrices (P and P_inv), initial abundance, and diagonal decay matrix.
    Used for solving Bateman equations in nuclear decay chains.
    """
    
    def __init__(
            self,
            decay_constants: np.ndarray,
            matrix_P: sparse.csc_matrix,
            matrix_P_inv: sparse.csc_matrix
        ) -> None:
        """
        Initializes the DecayMatrix with decay constants and transformation matrices.
        Args:
            decay_constants (np.ndarray): Array of decay constants for each nuclide.
            matrix_P (sparse.csc_matrix): Transformation matrix P (for Bateman solution).
            matrix_P_inv (sparse.csc_matrix): Inverse of transformation matrix P.
        """
        self.decay_constants = decay_constants
        self.matrix_P = matrix_P
        self.matrix_P_inv = matrix_P_inv
        # Initial abundance vector (all zeros by default)
        self.initial_abundance = self._setup_initial_abundance(matrix_P.shape[0])
        # Diagonal matrix for decay constants (all zeros by default, can be set later)
        self.matrix_Lambda = self._setup_matrix_Lambda(matrix_P.shape[0])
        

    @staticmethod
    def _setup_matrix_Lambda(size: int) -> sparse.csc_matrix:
        """
        Creates a diagonal matrix of zeros (placeholder for decay constants).
        Args:
            size (int): Size of the square matrix.
        Returns:
            sparse.csc_matrix: Diagonal sparse matrix of zeros.
        """
        # The diagonal matrix is used for holding decay constants in analytical solutions
        return sparse.csc_matrix(
            (np.zeros(size), (np.arange(size), np.arange(size))),
            shape=(size, size)
        )


    @staticmethod
    def _setup_initial_abundance(size: int) -> np.ndarray:
        """
        Creates an initial abundance vector of zeros.
        Args:
            size (int): Length of the abundance vector.
        Returns:
            np.ndarray: Zero-initialized abundance vector.
        """
        return np.zeros(size, dtype=np.float64)


    def __repr__(self) -> str:
        """
        Returns a string representation indicating the storage type and precision.
        """
        return f"DecayMatrixScipy: data stored in SciPy objects for double precision calculations."


    def __eq__(self, other: object) -> bool:
        """
        Compares two DecayMatrix objects for equality based on their matrices and initial abundance.
        Args:
            other (object): Another DecayMatrix instance.
        Returns:
            bool: True if all matrices and vectors are equal, False otherwise.
        """
        if not isinstance(other, self.__class__):
            return NotImplemented
        # Compare all relevant matrices and vectors for equality
        return (
            _csc_matrix_equal(self.matrix_P, other.matrix_P)
            and _csc_matrix_equal(self.matrix_P_inv, other.matrix_P_inv)
            and _csc_matrix_equal(self.matrix_Exp, other.matrix_Exp)
            and np.array_equal(self.vector_N0, other.vector_N0)
        )


    def __ne__(self, other: object) -> bool:
        """
        Checks inequality for DecayMatrix objects.
        Args:
            other (object): Another DecayMatrix instance.
        Returns:
            bool: True if not equal, False otherwise.
        """
        if not isinstance(other, self.__class__):
            return NotImplemented
        return not self.__eq__(other)


#---------------------------------
#    class MatrixBuilder
#---------------------------------
class MatrixBuilder:
    """
    Builds matrices representing the decay of nuclides based on decay constants, 
    branching ratios, and decay chains. Used for analytical and numerical solutions
    of nuclear decay networks (Bateman equations).

    Attributes:
        decay_database (DecayDatabase): Database containing nuclide decay data.
    """

    def __init__(self, decay_database: DecayDatabase) -> None:
        """
        Initializes the MatrixBuilder with a DecayDatabase instance.
        Args:
            decay_database (DecayDatabase): Contains data about nuclides and decay information.
        """
        self.decay_database = decay_database

    @property
    def matrix_size(self) -> int:
        """
        Returns the size of the decay matrix, which is the number of nuclides in the database.
        Returns:
            int: Number of nuclides.
        """
        return len(self.decay_database.nuclides)

    def build_decay_matrix(self) -> Tuple[sparse.csc_matrix, sparse.csc_matrix]:
        """
        Constructs both the decay matrix A and the transformation matrix P and its inverse.
        Returns:
            tuple: (matrix_P, matrix_P_inv)
                matrix_P (sparse.csc_matrix): Transformation matrix P.
                matrix_P_inv (sparse.csc_matrix): Inverse of matrix P.
        """
        matrix_A = self.build_matrix_A()
        matrix_P = self.build_matrix_P(matrix_A)
        matrix_P_inv = self.build_matrix_P_inv(matrix_A, matrix_P)
        return matrix_P, matrix_P_inv
    
    def build_matrix_A(self) -> sparse.csc_matrix:
        """
        Constructs the decay matrix A, which contains the decay constants for each nuclide.
        The matrix is built based on the decay constants derived from the half-lives of the nuclides,
        decay modes, progeny nuclides, and branching ratios. Diagonal elements are negative total decay rates;
        off-diagonal elements are positive and represent transfer to progeny.
        Returns:
            sparse.csc_matrix: Decay matrix A, where A[i, j] is the decay rate from nuclide j to i.
        """
        (
            half_life_data,
            decay_modes_data,
            progeny_data,
            branching_ratios_data
        ) = self._get_decay_data()

        matrix_A = sparse.lil_matrix((self.matrix_size, self.matrix_size), dtype=np.float64)

        for parent, parent_index in self.decay_database.nuclide_index_map.items():
            half_life = half_life_data[parent_index]
            # If nuclide is stable, decay constant is zero
            if half_life == np.inf:
                decay_constant = 0.0
            else:
                decay_constant = np.log(2) / half_life
            j = self.decay_database.nuclide_index_map[parent]
            matrix_A[j, j] = -decay_constant  # Diagonal: negative total decay rate
            # Fill off-diagonal elements for each decay branch
            for progeny, decay_mode, branching_ratio in zip(progeny_data[parent_index],
                                                            decay_modes_data[parent_index],
                                                            branching_ratios_data[parent_index]):
                # Only consider physical decay modes (skip stable, SF, etc.)
                if decay_mode in ['Stable', 'SF', 'Other', 'γ']:
                    continue
                if progeny not in self.decay_database.nuclides:
                    continue
                i = self.decay_database.nuclide_index_map[progeny]
                matrix_A[i, j] = decay_constant * float(branching_ratio)  # Off-diagonal: decay to progeny
        return matrix_A.tocsc()

    def build_matrix_P(self, matrix_A: sparse.csc_matrix) -> sparse.csc_matrix:
        """
        Constructs the transformation matrix P from the decay matrix A.
        This matrix is used in the analytical solution of the Bateman equations.
        Args:
            matrix_A (sparse.csc_matrix): Decay matrix A.
        Returns:
            sparse.csc_matrix: The transformation matrix P.
        Warnings:
            Prints a warning if the denominator in the calculation of P[i, j] is close to zero, which could lead to numerical instability.
        """
        # Start with identity matrix
        matrix_P = sparse.eye(self.matrix_size, self.matrix_size,
                              dtype=np.float64, format='lil')
        # Build the decay chain indices for rows and columns, which will be used in the recursive calculation
        row_indices_dict, rows_P, cols_P = self._build_decay_chain_indices(matrix_A)
        # Loop over the rows and columns in the decay chain to calculate the values of P
        for index in range(0, rows_P.size):
            i = rows_P[index]
            j = cols_P[index]
            if i == j:
                continue
            Sigma_sum = 0.0
            # Recursively sum over all possible decay paths
            for k in row_indices_dict[j]:
                if k == i:
                    break
                Sigma_sum += matrix_A[i, k] * matrix_P[k, j]
            denominator = matrix_A[j, j] - matrix_A[i, i]
            try:
                if abs(denominator) == 0.0:
                    raise ZeroDivisionError(
                        f'Denominator too small for: '
                        f'i={i} ({self.decay_database.nuclides[i]}), '
                        f'j={j} ({self.decay_database.nuclides[j]}), '
                        f'({matrix_A[j, j]}) - ({matrix_A[i, i]}) = {denominator}.\n'
                    )
                matrix_P[i, j] = Sigma_sum / denominator
            except ZeroDivisionError as e:
                print(f'Warning: {e} Setting matrix_P[{i}, {j}] to 1.0e-10 due to division by zero.')
                matrix_P[i, j] = 1.0e-10
        return matrix_P.tocsc()

    def build_matrix_P_inv(
            self,
            matrix_A: sparse.csc_matrix,
            matrix_P: sparse.csc_matrix
        ) -> sparse.csc_matrix:
        """
        Constructs the inverse of the transformation matrix P.
        This matrix is computed recursively using the elements of matrix P and the decay rates in matrix A.
        Args:
            matrix_A (sparse.csc_matrix): Decay matrix A.
            matrix_P (sparse.csc_matrix): Transformation matrix P.
        Returns:
            sparse.csc_matrix: The inverse of the transformation matrix P.
        """
        matrix_P_inv = sparse.eye(self.matrix_size, self.matrix_size, dtype=np.float64, format="lil")
        row_indices_dict, rows_P, cols_P = self._build_decay_chain_indices(matrix_A)
        for index in range(0, rows_P.size):
            i = rows_P[index]
            j = cols_P[index]
            if i == j:
                continue
            Sigma_sum = 0.0
            # Recursively sum for the inverse
            for k in row_indices_dict[j]:
                if k == i:
                    break
                Sigma_sum -= matrix_P[i, k] * matrix_P_inv[k, j]
            matrix_P_inv[i, j] = Sigma_sum
        return matrix_P_inv.tocsc()

    def _get_decay_data(self) -> Tuple:
        """
        Retrieves the necessary decay data from the decay database.
        Returns:
            tuple: (half_life_data, decay_modes_data, progeny_data, branching_ratios_data)
        """
        return (self.decay_database.half_life_data[:, 0],
                self.decay_database.decay_modes_data,
                self.decay_database.progeny_data,
                self.decay_database.branching_ratios_data
            )

    def _build_decay_chain_indices(self, matrix_A: sparse.csc_matrix) -> Tuple[Dict, np.ndarray, np.ndarray]:
        """
        Builds decay chain indices for the recursive calculation of matrix P.
        Creates a dictionary mapping each nuclide to the set of nuclides that decay into it.
        Args:
            matrix_A (sparse.csc_matrix): Decay matrix A.
        Returns:
            tuple:
                - dict: Mapping from nuclide index to set of indices that decay into it.
                - np.ndarray: Row indices.
                - np.ndarray: Column indices.
        """
        # A dictionary to store decay chain members by index
        row_indices_dict = {}
        for i in range(self.matrix_size-1, -1, -1):
            nonzero_row_indices, _ = matrix_A[:, i].nonzero()
            row_indices_dict[i] = set(nonzero_row_indices)
            # Include decay chain of progeny nuclides
            for j in nonzero_row_indices:
                if j > i:
                    row_indices_dict[i].update(row_indices_dict.get(j, set()))
            # Sort the set of decay chain indices
            row_indices_dict[i] = np.array(sorted(row_indices_dict[i]))
        rows_P = []
        cols_P = []
        # Build row and column indices for the matrix P
        for i in range(self.matrix_size):
            rows_P.extend(row_indices_dict[i])
            cols_P.extend([i] * len(row_indices_dict[i]))
        rows_P = np.array(rows_P, dtype=np.int32)
        cols_P = np.array(cols_P, dtype=np.int32)
        return row_indices_dict, rows_P, cols_P

    def _check_identity(self):
        """
        Checks if the matrix P and its inverse P_inv form an identity matrix.
        Returns:
            bool: True if the product of P and P_inv is an identity matrix, False otherwise.
        """
        matrix_P, matrix_P_inv = self.build_decay_matrix()
        csc = matrix_P @ matrix_P_inv
        if csc.shape[0] != csc.shape[1]:
            return False
        # Check if the product is exactly the identity matrix
        if all(csc.diagonal() == 1) and csc.nnz == csc.shape[0]:
            return True
        else:
            return False


def _csc_matrix_equal(matrix_a: sparse.csc_matrix,
                      matrix_b: sparse.csc_matrix) -> bool:
    """
    Helper function to compare two scipy csc matrices for exact equality.
    Args:
        matrix_a (sparse.csc_matrix): First matrix.
        matrix_b (sparse.csc_matrix): Second matrix.
    Returns:
        bool: True if matrices are exactly equal, False otherwise.
    """
    return (
        np.array_equal(matrix_a.indptr, matrix_b.indptr)
        and np.array_equal(matrix_a.indices, matrix_b.indices)
        and np.array_equal(matrix_a.data, matrix_b.data)
    )
    
    

#---------------------------------
#    class DecayDiagram
#---------------------------------
class DecayDiagram:
    """
    Visualizes nuclear decay chains and their reverse (parent) chains for a given nuclide.
    Provides methods to plot decay diagrams and nuclear charts using networkx and matplotlib.
    """

    def __init__(self, decay_database):
        """
        Initialize the DecayDiagram with a decay database.
        Args:
            decay_database: An instance of DecayDatabase containing decay data.
        """
        self.decay_database = decay_database
            
    
    def find_daughters(self, nuclide: str) -> List[str]:
        """
        Find all possible daughter nuclides that can be produced from the given nuclide.
        Args:
            nuclide (str): Parent nuclide symbol.
        Returns:
            list[str]: List of daughter nuclide symbols.
        """
        return self.decay_database.get_nuclide_progeny(nuclide)
    

    def find_parents(self, nuclide: str) -> List[str]:
        """
        Find all possible parent nuclides that can decay to the given nuclide.
        Args:
            nuclide (str): Daughter nuclide symbol.
        Returns:
            list[str]: List of parent nuclide symbols.
        """
        return self.decay_database.get_nuclide_ancestor(nuclide)

    

    def plot_nuclear_chart(
        self,
        figure: Optional[plt.figure] = None,
        nuclei_linewidths: float = 0.3,
        colorbar: bool = False,
        figsize: Tuple[float, float] = (9, 6),
        dpi: int = 300,
        magic_numbers: List[int] = [2, 8, 20, 28, 50, 82, 126],
        **kwargs
    ) -> plt.figure:
        """
        Plot the nuclear chart (N vs Z) colored by half-life, with magic numbers highlighted.
        Args:
            figure (plt.figure, optional): Existing matplotlib figure to plot on.
            nuclei_linewidths (float): Line width for nuclide boxes.
            colorbar (bool): Whether to show colorbar.
            figsize (Tuple[float, float]): Figure size.
            dpi (int): Dots per inch for the figure.
            magic_numbers (List[int]): List of magic numbers to highlight.
            **kwargs: Additional keyword arguments for plotting.
        Returns:
            plt.figure: The matplotlib figure object.
        Raises:
            ValueError: If no valid nuclide or half-life data is available.
        """

        colormap = HalfLifeColorMap(self.decay_database.half_life_data)
        
        # Group nuclides by half-life ranges
        range_groups = {key: [] for key in colormap.ranges.keys()}
        stable_nuclides = []
        
        # Process nuclides and group them
        for nuclide in self.decay_database.nuclides:
            try:
                parsed = Nuclide(nuclide)
                Z = parsed.Z
                N = parsed.N
                state = parsed.state
                if state:
                    continue  # Skip metastable states
                
                half_life = self.decay_database.get_nuclide_half_life(nuclide, units="s")
                
                if half_life == float('inf'):
                    stable_nuclides.append((N, Z))
                else:
                    range_key = colormap.get_range_for_half_life(half_life)
                    if range_key:
                        range_groups[range_key].append((N, Z, half_life))
                    
            except Exception:
                continue

        if not any(range_groups.values()) and not stable_nuclides:
            raise ValueError("No valid nuclide data for plotting.")

        # Create figure
        fig, ax = colormap._create_figure(figure, figsize, dpi)
        
        # Plot stable nuclei
        colormap._plot_stable_nuclei(ax, stable_nuclides, nuclei_linewidths)
        
        # Plot unstable nuclei by range
        colormap._plot_unstable_nuclei(ax, range_groups, colormap, nuclei_linewidths)
        
        # Get all N and Z values for axis limits
        all_N_values, all_Z_values = colormap._get_axis_limits(range_groups, stable_nuclides)
        
        # Mark magic numbers
        colormap._mark_magic_numbers(ax, magic_numbers, range_groups, stable_nuclides)
        
        # Set axis limits and labels
        colormap._set_axis_properties(ax, all_N_values, all_Z_values)

        colormap._plot_legend(ax, range_groups, stable_nuclides)

        return fig
    

    def plot_decay_chains(
        self,
        nuclide: str,
        label_pos: float = 0.3,
        fig: Optional[mpl.figure.Figure] = None,
        axes: Optional[mpl.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:
        """
        Plot the decay chains starting from the specified nuclide.
        Args:
            nuclide (str): Starting nuclide symbol.
            label_pos (float): Position of edge labels along the edge.
            fig (Optional[Figure]): Existing matplotlib figure to plot on.
            axes (Optional[Axes]): Existing matplotlib axes to plot on.
            kwargs_draw (Optional[Dict]): Additional keyword arguments for nx.draw.
            kwargs_edge_labels (Optional[Dict]): Additional keyword arguments for nx.draw_networkx_edge_labels.
        Returns:
            Tuple[Figure, Axes]: The matplotlib figure and axes objects.
        """
        nuclide = Nuclide(nuclide).nuclide_symbol
        digraph, max_generation, max_xpos = self._build_decay_digraph(nuclide, nx.DiGraph())
        return self.plot(
            digraph,
            max_generation,
            max_xpos,
            label_pos,
            fig,
            axes,
            kwargs_draw,
            kwargs_edge_labels,
        )
    

    def plot_reverse_decay_chains(
        self,
        nuclide: str,
        label_pos: float = 0.3,
        fig: Optional[mpl.figure.Figure] = None,
        axes: Optional[mpl.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:
        """
        Plot the reverse decay chains (parents) for the specified nuclide.
        Args:
            nuclide (str): Target nuclide symbol.
            label_pos (float): Position of edge labels along the edge.
            fig (Optional[Figure]): Existing matplotlib figure to plot on.
            axes (Optional[Axes]): Existing matplotlib axes to plot on.
            kwargs_draw (Optional[Dict]): Additional keyword arguments for nx.draw.
            kwargs_edge_labels (Optional[Dict]): Additional keyword arguments for nx.draw_networkx_edge_labels.
        Returns:
            Tuple[Figure, Axes]: The matplotlib figure and axes objects.
        """
        nuclide = Nuclide(nuclide).nuclide_symbol
        digraph, max_generation, max_xpos = self._build_reverse_decay_digraph(nuclide, nx.DiGraph())
        return self.plot(
            digraph,
            max_generation,
            max_xpos,
            label_pos,
            fig,
            axes,
            kwargs_draw,
            kwargs_edge_labels,
        )
        
    def plot_rProcess_chains(
        self,
        nuclide: str,
        label_pos: float = 0.3,
        fig: Optional[mpl.figure.Figure] = None,
        axes: Optional[mpl.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:
        """
        Plot the r-process chains starting from the specified nuclide.
        Args:
            nuclide (str): Starting nuclide symbol.
            label_pos (float): Position of edge labels along the edge.
            fig (Optional[Figure]): Existing matplotlib figure to plot on.
            axes (Optional[Axes]): Existing matplotlib axes to plot on.
            kwargs_draw (Optional[Dict]): Additional keyword arguments for nx.draw.
            kwargs_edge_labels (Optional[Dict]): Additional keyword arguments for nx.draw_networkx_edge_labels.
        Returns:
            Tuple[Figure, Axes]: The matplotlib figure and axes objects.
        """
        nuclide = Nuclide(nuclide).nuclide_symbol
        digraph, max_generation, max_xpos = self._build_rProcess_reverse_decay_digraph(nuclide, nx.DiGraph())
        return self.plot(
            digraph,
            max_generation,
            max_xpos,
            label_pos,
            fig,
            axes,
            kwargs_draw,
            kwargs_edge_labels,
        )


    @staticmethod
    def plot(
        digraph: nx.DiGraph,
        max_generation: int,
        max_xpos: int,
        label_pos: float = 0.3,
        fig: Optional[mpl.figure.Figure] = None,
        axes: Optional[mpl.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:
        """
        Plot the decay diagram using networkx and matplotlib.
        Args:
            digraph (DiGraph): The networkx DiGraph to plot.
            max_generation (int): Maximum generation in the decay chain.
            max_xpos (int): Maximum x position in the decay chain.
            label_pos (float): Position of edge labels along the edge.
            fig (Optional[Figure]): Existing matplotlib figure to plot on.
            axes (Optional[Axes]): Existing matplotlib axes to plot on.
            kwargs_draw (Optional[Dict]): Additional keyword arguments for nx.draw.
            kwargs_edge_labels (Optional[Dict]): Additional keyword arguments for nx.draw_networkx_edge_labels.
        Returns:
            Tuple[Figure, Axes]: The matplotlib figure and axes objects.
        """
        positions = nx.get_node_attributes(digraph, "pos")
        node_labels = nx.get_node_attributes(digraph, "label")
        edge_labels = nx.get_edge_attributes(digraph, "label")

        # Prepare figure and axes
        fig, axes = _check_fig_axes(
            fig, axes, figsize=(3 * max_xpos + 1.5, 3 * max_generation + 1.5)
        )

        # Set default drawing parameters if not provided
        kwargs_draw = kwargs_draw or {}
        kwargs_draw.setdefault("node_size", 6000)
        kwargs_draw.setdefault("node_color", "#1f77b4")
        kwargs_draw.setdefault("edgecolors", "#000000")
        kwargs_draw.setdefault("edge_color", "#000000")
        kwargs_draw.setdefault("node_shape", "o")
        kwargs_draw.setdefault("alpha", 0.8)
        nx.draw(
            G = digraph,
            pos = positions,
            ax = axes,
            labels = node_labels,
            **kwargs_draw,
        )

        # Set default edge label parameters if not provided
        kwargs_edge_labels = kwargs_edge_labels or {}
        kwargs_edge_labels.setdefault("font_size", 14)
        # kwargs_edge_labels.setdefault("font_color", "red")
        # kwargs_edge_labels.setdefault("font_weight", "bold")
        kwargs_edge_labels.setdefault("bbox", {"boxstyle": None, "ec": (1.0, 1.0, 1.0), "fc": (1.0, 1.0, 1.0)})
        kwargs_edge_labels.setdefault("rotate", False)

        nx.draw_networkx_edge_labels(
            G = digraph,
            pos = positions,
            edge_labels = edge_labels,
            label_pos = label_pos,
            ax = axes,
            **kwargs_edge_labels,
        )

        axes.set_xlim(-0.3, max_xpos + 0.3)
        axes.set_ylim(-max_generation - 0.3, 0.3)

        return fig, axes


    def _build_decay_digraph(
        self,
        nuclide: str,
        digraph: nx.DiGraph
    ) -> nx.DiGraph:
        """
        Build the decay directed graph starting from the specified nuclide.
        Args:
            nuclide (str): Starting nuclide symbol.
            digraph (DiGraph): The networkx DiGraph to build.
        Returns:
            Tuple[DiGraph, int, int]: The graph, maximum generation, and maximum x position.
        """
        nuclide = Nuclide(nuclide).nuclide_symbol
        generation_max_xpos = {0: 0}  # Track max x position for each generation
        dequeue = deque([nuclide])  # Queue for BFS
        generations = deque([0])         # Track generation (depth)
        xpositions = deque([0])          # Track x position for layout
        node_label = (
            _parse_nuclide_label(nuclide)
            + '\n' 
            + str(self.decay_database.get_nuclide_half_life(nuclide, 'readable'))
        )
        
        digraph.add_node(nuclide, generation=0, xpos=0, label=node_label)
        seen = {nuclide}  # Track visited nodes

        while dequeue:
            parent_nuclide = dequeue.popleft()
            generation = generations.popleft() + 1
            xpos = xpositions.popleft()
            if generation not in generation_max_xpos:
                generation_max_xpos[generation] = -1

            daughter_nuclides = self.find_daughters(parent_nuclide)
            xpos = max(xpos, generation_max_xpos[generation] + 1)
            xcounter = 0
            for daughter in daughter_nuclides:
                decay_mode, branching_ratio = self.decay_database.get_decay_channel(parent_nuclide, daughter)
                # Skip spontaneous fission (SF) nodes
                if daughter == "SF" or decay_mode == "SF":
                    continue
                # If daughter not seen, add to graph and queue
                if daughter not in seen:
                    node_label = _parse_nuclide_label(daughter)
                    if daughter in self.decay_database.nuclide_index_map:
                        node_label += f'\n{self.decay_database.get_nuclide_half_life(daughter, "readable")}'
                        if np.isfinite(self.decay_database.get_nuclide_half_life(daughter, 's')):
                            dequeue.append(daughter)
                            generations.append(generation)
                            xpositions.append(xpos + xcounter)
                    digraph.add_node(
                        daughter,
                        generation=generation,
                        xpos=xpos + xcounter,
                        label=node_label,
                    )
                    seen.add(daughter)
                    if xpos + xcounter > generation_max_xpos[generation]:
                        generation_max_xpos[generation] = xpos + xcounter
                    xcounter += 1
                # Add edge with decay mode and branching ratio label
                edge_label = (
                    _parse_decay_mode_label(decay_mode)
                    + '\n'
                    + f'{branching_ratio*100:.2f}%'
                )
                digraph.add_edge(Nuclide(parent_nuclide).nuclide_symbol, daughter, label=edge_label)

        # Assign positions for plotting
        for node in digraph:
            digraph.nodes[node]["pos"] = (
                digraph.nodes[node]["xpos"],
                digraph.nodes[node]["generation"] * -1,
            )

        return digraph, max(generation_max_xpos), max(generation_max_xpos.values())
    

    def _build_reverse_decay_digraph(
        self,
        nuclide: str,
        digraph: nx.DiGraph
    ) -> Tuple[nx.DiGraph, int, int]:
        """
        Build the reverse decay directed graph (parents to the specified nuclide).
        Args:
            nuclide (str): Target nuclide symbol.
            digraph (DiGraph): The networkx DiGraph to build.
        Returns:
            Tuple[DiGraph, int, int]: The graph, maximum generation, and maximum x position.
        """
        nuclide = Nuclide(nuclide).nuclide_symbol
        generation_max_xpos = {0: 0}
        dequeue = deque([nuclide])
        generations = deque([0])
        xpositions = deque([0])
        node_label = (
            _parse_nuclide_label(nuclide)
            + '\n' 
            + str(self.decay_database.get_nuclide_half_life(nuclide, 'readable'))
        )
        digraph.add_node(nuclide, generation=0, xpos=0, label=node_label)
        seen = {nuclide}

        while dequeue:
            daughter_nuclide = dequeue.popleft()
            generation = generations.popleft() + 1
            xpos = xpositions.popleft()
            if generation not in generation_max_xpos:
                generation_max_xpos[generation] = -1

            parent_nuclides = self.find_parents(daughter_nuclide)
            xpos = max(xpos, generation_max_xpos[generation] + 1)
            xcounter = 0
            for parent in parent_nuclides:
                if parent not in seen:
                    node_label = _parse_nuclide_label(parent)
                    if parent in self.decay_database.nuclide_index_map:
                        node_label += f'\n{self.decay_database.get_nuclide_half_life(parent, "readable")}'
                        if np.isfinite(self.decay_database.get_nuclide_half_life(parent, 's')):
                            dequeue.append(parent)
                            generations.append(generation)
                            xpositions.append(xpos + xcounter)
                    digraph.add_node(
                        parent,
                        generation=generation,
                        xpos=xpos + xcounter,
                        label=node_label,
                    )
                    seen.add(parent)
                    if xpos + xcounter > generation_max_xpos[generation]:
                        generation_max_xpos[generation] = xpos + xcounter
                    xcounter += 1
                decay_mode, branching_ratio = self.decay_database.get_decay_channel(parent, daughter_nuclide)
                if decay_mode:
                    edge_label = (
                        _parse_decay_mode_label(decay_mode)
                        + '\n'
                        + f'{branching_ratio*100:.2f}%'
                    )
                    digraph.add_edge(parent, daughter_nuclide, label=edge_label)

        # Assign positions for plotting
        for node in digraph:
            digraph.nodes[node]["pos"] = (
                digraph.nodes[node]["xpos"],
                digraph.nodes[node]["generation"] * -1,
            )

        return digraph, max(generation_max_xpos), max(generation_max_xpos.values()) 



    def _build_rProcess_reverse_decay_digraph(
        self,
        nuclide: str,
        digraph: nx.DiGraph
    ) -> Tuple[nx.DiGraph, int, int]:
        """
        Build the r-process reverse decay directed graph (parents) for the specified nuclide.
        Args:
            nuclide (str): Target nuclide symbol.
            digraph (DiGraph): The networkx DiGraph to build.
        Returns:
            Tuple[DiGraph, int, int]: The graph, maximum generation, and maximum x position.
        """
        nuclide = Nuclide(nuclide).nuclide_symbol
        generation_max_xpos = {0: 0}
        dequeue = deque([nuclide])
        generations = deque([0])
        xpositions = deque([0])
        node_label = (
            _parse_nuclide_label(nuclide)
            + '\n' 
            + str(self.decay_database.get_nuclide_half_life(nuclide, 'readable'))
        )
        digraph.add_node(nuclide, generation=0, xpos=0, label=node_label)
        seen = {nuclide}

        while dequeue:
            daughter_nuclide = dequeue.popleft()
            generation = generations.popleft() + 1
            xpos = xpositions.popleft()
            if generation not in generation_max_xpos:
                generation_max_xpos[generation] = -1

            parent_nuclides = self.find_parents(daughter_nuclide)
            xpos = max(xpos, generation_max_xpos[generation] + 1)
            xcounter = 0
            for parent in parent_nuclides:
                decay_mode, branching_ratio = self.decay_database.get_decay_channel(parent, daughter_nuclide)
                if decay_mode in ['β+&EC', 'SF', 'IT']:
                    continue
                # if self.decay_database.get_nuclide_half_life(parent, 's') < 1.0:
                #     continue
                if parent not in seen:
                    node_label = _parse_nuclide_label(parent)
                    if parent in self.decay_database.nuclide_index_map:
                        node_label += f'\n{self.decay_database.get_nuclide_half_life(parent, "readable")}'
                        if np.isfinite(self.decay_database.get_nuclide_half_life(parent, 's')):
                            dequeue.append(parent)
                            generations.append(generation)
                            xpositions.append(xpos + xcounter)
                    digraph.add_node(
                        parent,
                        generation=generation,
                        xpos=xpos + xcounter,
                        label=node_label,
                    )
                    seen.add(parent)
                    if xpos + xcounter > generation_max_xpos[generation]:
                        generation_max_xpos[generation] = xpos + xcounter
                    xcounter += 1
                if decay_mode:
                    edge_label = (
                        _parse_decay_mode_label(decay_mode)
                        + '\n'
                        + f'{branching_ratio*100:.2f}%'
                    )
                    digraph.add_edge(parent, daughter_nuclide, label=edge_label)

        # Assign positions for plotting
        for node in digraph:
            digraph.nodes[node]["pos"] = (
                digraph.nodes[node]["xpos"],
                digraph.nodes[node]["generation"] * -1,
            )

        return digraph, max(generation_max_xpos), max(generation_max_xpos.values()) 
        

    def plot_local_nuclear_chart_for_rprocess_chain(
        self,
        nuclide: str,
        figure: Optional[plt.figure] = None,
        nuclei_linewidths: float = 0.3,
        figsize: Tuple[float, float] = (9, 6),
        dpi: int = 300,
        magic_numbers: List[int] = [2, 8, 20, 28, 50, 82, 126],
    ) -> plt.figure:
        """
        Plot the local nuclear chart for the r-process chain involving the specified nuclide.
        Args:
            nuclide (str): The starting nuclide symbol.
            figure (plt.figure, optional): Existing matplotlib figure.
            nuclei_linewidths (float): Line width for nuclei.
            figsize (Tuple[float, float]): Figure size.
            dpi (int): Resolution.
            magic_numbers (List[int]): Magic numbers.
        Returns:
            plt.figure: matplotlib figure object.
        """
        digraph, _, _ = self._build_rProcess_reverse_decay_digraph(nuclide, nx.DiGraph())
        nuclides = list(digraph.nodes)
        unstable = []
        stable = []
        for n in nuclides:
            try:
                parsed = Nuclide(n)
                Z = parsed.Z
                N = parsed.N
                state = parsed.state
                if state:
                    continue 
                half_life = self.decay_database.get_nuclide_half_life(n, units="s")
                if half_life == float('inf'):
                    stable.append((N, Z))
                else:
                    unstable.append((N, Z, half_life))
            except Exception:
                continue

        colormap = HalfLifeColorMap(self.decay_database.half_life_data)
        range_groups = {key: [] for key in colormap.ranges.keys()}
        for N, Z, half_life in unstable:
            range_key = colormap.get_range_for_half_life(half_life)
            if range_key:
                range_groups[range_key].append((N, Z, half_life))

        fig, ax = colormap._create_figure(figure, figsize, dpi)
        colormap._plot_stable_nuclei(ax, stable, nuclei_linewidths)
        colormap._plot_unstable_nuclei(ax, range_groups, colormap, nuclei_linewidths)
        pos = {}
        for n in nuclides:
            try:
                parsed = Nuclide(n)
                pos[n] = (parsed.N, parsed.Z)
            except Exception:
                continue
        nx.draw_networkx_edges(
            digraph, pos, ax=ax, arrows=True, arrowstyle='-|>', arrowsize=10, edge_color='black', width=1.0
        )
        edge_labels = nx.get_edge_attributes(digraph, "label")
        nx.draw_networkx_edge_labels(
            digraph, pos, edge_labels=edge_labels, ax=ax,
            font_size=3, font_color='black', label_pos=0.5, rotate=False,
            bbox={"boxstyle": "round,pad=0.1", "fc": "white", "ec": "none", "alpha": 0.7}
        )
        all_N_values, all_Z_values = colormap._get_axis_limits(range_groups, stable)
        colormap._mark_magic_numbers(ax, magic_numbers, range_groups, stable)
        colormap._set_axis_properties(ax, all_N_values, all_Z_values)
        colormap._plot_legend(ax, range_groups, stable)
        return fig


def _check_fig_axes(
    fig_in: Optional[mpl.figure.Figure],
    axes_in: Optional[mpl.axes.Axes],
    **kwargs,
) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:
    """
    Check and create matplotlib figure and axes if needed for plotting.
    Args:
        fig_in (Optional[Figure]): Input figure.
        axes_in (Optional[Axes]): Input axes.
        **kwargs: Additional arguments for plt.subplots.
    Returns:
        Tuple[Figure, Axes]: The figure and axes to use.
    """
    if fig_in is None and axes_in is None:
        fig, axes = plt.subplots(**kwargs)
    elif fig_in is None:
        axes = axes_in
        fig = axes.get_figure()
    elif axes_in is None:
        fig = fig_in
        axes = fig.gca()
    else:
        fig = fig_in
        axes = axes_in
    return fig, axes


def _parse_nuclide_label(nuclide: str) -> str:
    """
    Format a nuclide symbol into a label with superscript mass number for plotting.
    Args:
        nuclide (str): Nuclide symbol in format 'Element-MassNumber'.
    Returns:
        str: Formatted nuclide label with superscript mass number.
    Raises:
        ValueError: If the nuclide format is invalid.
    """
    if nuclide == "SF":
        return "SF"
    if nuclide == "X":
        return "SF"
    if "-" not in nuclide:
        raise ValueError(f"Invalid nuclide format: {nuclide}. Expected format 'Element-MassNumber'.")
    nuclide_unicode_map = {
        "0": "\N{SUPERSCRIPT ZERO}",
        "1": "\N{SUPERSCRIPT ONE}",
        "2": "\N{SUPERSCRIPT TWO}",
        "3": "\N{SUPERSCRIPT THREE}",
        "4": "\N{SUPERSCRIPT FOUR}",
        "5": "\N{SUPERSCRIPT FIVE}",
        "6": "\N{SUPERSCRIPT SIX}",
        "7": "\N{SUPERSCRIPT SEVEN}",
        "8": "\N{SUPERSCRIPT EIGHT}",
        "9": "\N{SUPERSCRIPT NINE}",
        "M": "\N{MODIFIER LETTER SMALL M}",
        "N": "\N{SUPERSCRIPT LATIN SMALL LETTER N}",
        "O": "\N{MODIFIER LETTER SMALL O}",
        "P": "\N{MODIFIER LETTER SMALL P}",
        "Q": "\N{LATIN SMALL LETTER Q}",
        "R": "\N{MODIFIER LETTER SMALL R}",
        "X": "\N{MODIFIER LETTER SMALL X}",
    }
    element_symbol, mass_number = nuclide.split("-")
    mass_number_unicode = "".join([nuclide_unicode_map.get(char, char) for char in mass_number])
    return mass_number_unicode + element_symbol


def _parse_decay_mode_label(decay_mode: str) -> str:
    """
    Format a decay mode string into a label with unicode characters for plotting.
    Args:
        decay_mode (str): Decay mode symbol.
    Returns:
        str: Formatted decay mode label with unicode characters.
    """
    decay_mode_unicode_map = {
        "α": "\N{GREEK SMALL LETTER ALPHA}",
        "β": "\N{GREEK SMALL LETTER BETA}",
        "ε": "\N{GREEK SMALL LETTER EPSILON}",
        "+": "\N{SUPERSCRIPT PLUS SIGN}",
        "-": "\N{SUPERSCRIPT MINUS}",
    }
    for unformatted, formatted in decay_mode_unicode_map.items():
        decay_mode = decay_mode.replace(unformatted, formatted)
    return decay_mode
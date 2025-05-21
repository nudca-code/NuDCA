# -*- coding: utf-8 -*-

from typing import Dict, Tuple, List
from pathlib import Path
import numpy as np
from scipy import sparse

from .decay_database import DecayDatabase, load_decay_database


class DecayMatrix:
    """
    Represents the decay matrix and related structures for nuclide decay calculations.
    Stores decay constants, transformation matrices (P and P_inv), and provides initial abundance and diagonal decay matrix.
    Used for solving the Bateman equations in nuclear decay chains.
    """
    
    def __init__(
            self,
            decay_constants: np.ndarray,
            matrix_P: sparse.csc_matrix,
            matrix_P_inv: sparse.csc_matrix
        ) -> None:
        """
        Initialize the DecayMatrix with decay constants and transformation matrices.
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
        Create a diagonal matrix of zeros (placeholder for decay constants).
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
        Create an initial abundance vector of zeros.
        Args:
            size (int): Length of the abundance vector.
        Returns:
            np.ndarray: Zero-initialized abundance vector.
        """
        return np.zeros(size, dtype=np.float64)


    def __repr__(self) -> str:
        """
        Return a string representation indicating the storage type and precision.
        """
        return f"DecayMatrixScipy: data stored in SciPy objects for double precision calculations."


    def __eq__(self, other: object) -> bool:
        """
        Compare two DecayMatrix objects for equality based on their matrices and initial abundance.
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
        Inequality comparison for DecayMatrix objects.
        Args:
            other (object): Another DecayMatrix instance.
        Returns:
            bool: True if not equal, False otherwise.
        """
        if not isinstance(other, self.__class__):
            return NotImplemented
        return not self.__eq__(other)


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
    def nuclides(self) -> List:
        """
        Returns the list of nuclides from the decay database.
        """
        return self.decay_database.nuclides
    
    @property
    def nuclide_index_map(self) -> Dict:
        """
        Returns the mapping of nuclide names to their respective indices.
        """
        return self.decay_database.nuclide_index_map
    
    @property
    def matrix_size(self) -> int:
        """
        Returns the size of the decay matrix, which is the number of nuclides in the database.
        """
        return len(self.nuclides)

    def build_decay_matrix(self) -> Tuple[sparse.csc_matrix, sparse.csc_matrix]:
        """
        Constructs both the decay matrix A and the matrix P and its inverse.
        Returns:
            Tuple[sparse.csc_matrix, sparse.csc_matrix]:
                A tuple containing the matrix P and its inverse P_inv.
        """
        matrix_A = self.build_matrix_A()
        matrix_P = self.build_matrix_P(matrix_A)
        matrix_P_inv = self.build_matrix_P_inv(matrix_A, matrix_P)
        return matrix_P, matrix_P_inv
    
    def build_matrix_A(self) -> sparse.csc_matrix:
        """
        Constructs the decay matrix A, which contains the decay constants for each nuclide.
        The matrix is built based on the decay constants derived from the half-lives of the nuclides,
        decay modes, progeny nuclides, and branching ratios. It represents the decay rates from one nuclide 
        to another. Diagonal elements are negative total decay rates; off-diagonal elements are positive
        and represent transfer to progeny.
        Returns:
            sparse.csc_matrix: A sparse matrix representing the decay matrix A, where each element A[i, j] 
            represents the decay rate from nuclide j to nuclide i.
        """
        (
            half_life_data,
            decay_modes_data,
            progeny_data,
            branching_ratios_data
        ) = self._get_decay_data()

        matrix_A = sparse.lil_matrix((self.matrix_size, self.matrix_size), dtype=np.float64)

        for parent, parent_index in self.nuclide_index_map.items():
            half_life = half_life_data[parent_index]
            # If nuclide is stable, decay constant is zero
            if half_life == np.inf:
                decay_constant = 0.0
            else:
                decay_constant = np.log(2) / half_life
            j = self.nuclide_index_map[parent]
            matrix_A[j, j] = -decay_constant  # Diagonal: negative total decay rate
            # Fill off-diagonal elements for each decay branch
            for progeny, decay_mode, branching_ratio in zip(progeny_data[parent_index],
                                                            decay_modes_data[parent_index],
                                                            branching_ratios_data[parent_index]):
                # Only consider physical decay modes (skip stable, SF, etc.)
                if decay_mode in ['Stable', 'SF', 'Other', 'γ']:
                    continue
                if progeny not in self.nuclides:
                    continue
                i = self.nuclide_index_map[progeny]
                matrix_A[i, j] = decay_constant * float(branching_ratio)  # Off-diagonal: decay to progeny
        return matrix_A.tocsc()

    def build_matrix_P(self, matrix_A: sparse.csc_matrix) -> sparse.csc_matrix:
        """
        Constructs the matrix P from the decay matrix A.
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
                    raise ZeroDivisionError(f'Denominator too small for i={i} ({self.nuclides[i]}), j={j} ({self.nuclides[j]}): '
                                            f'({matrix_A[j, j]}) - ({matrix_A[i, i]}) = {denominator}.\n')
                matrix_P[i, j] = Sigma_sum / denominator
            except ZeroDivisionError as e:
                print(f'Warning: {e} Setting matrix_P[{i}, {j}] to  due to division by zero.')
                # matrix_P[i, j] = 0.0
            matrix_P[i, j] = Sigma_sum / denominator
        return matrix_P.tocsc()

    def build_matrix_P_inv(self,
                            matrix_A: sparse.csc_matrix,
                            matrix_P: sparse.csc_matrix) -> sparse.csc_matrix:
        """
        Constructs the inverse of the matrix P.
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
            Tuple: (half_life_data, decay_modes_data, progeny_data, branching_ratios_data)
        """
        return (self.decay_database.half_life_data[:, 0],
                self.decay_database.decay_modes_data,
                self.decay_database.progeny_data,
                self.decay_database.branching_ratios_data
            )

    def _build_decay_chain_indices(self, matrix_A: sparse.csc_matrix) -> Tuple[Dict, np.ndarray, np.ndarray]:
        """
        Builds decay chain indices for the recursive calculation of matrix P.
        This function creates a dictionary mapping each nuclide to the set of nuclides that decay into it.
        Args:
            matrix_A (sparse.csc_matrix): Decay matrix A.
        Returns:
            Tuple[Dict, np.ndarray, np.ndarray]:
                - Dictionary mapping each nuclide to a set of nuclides that decay into it.
                - Numpy array of row indices.
                - Numpy array of column indices.
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


def save_decay_matrix(data_source: str='ENDF-B-VIII.1_decay'):
    """
    Save the decay matrix P and its inverse to disk as .npz files for a given data source.
    Args:
        data_source (str): Name of the decay data source.
    """
    decay_database = load_decay_database(data_source)
    matrix_P, matrix_P_inv = MatrixBuilder(decay_database).build_decay_matrix()
    # Save the matrices in compressed sparse format for efficient storage
    sparse.save_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P.npz'), matrix_P.tocsr())
    sparse.save_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P_inv.npz'), matrix_P_inv.tocsr())


def load_decay_matrix(
        data_source: str='ENDF-B-VIII.1_decay'
    ) -> DecayMatrix:
    """
    Load the decay matrix P and its inverse from disk and return a DecayMatrix object.
    Args:
        data_source (str): Name of the decay data source.
    Returns:
        DecayMatrix: The loaded decay matrix object.
    """
    # Optionally regenerate and save the matrices if needed
    # save_decay_matrix(data_source)
    decay_database = load_decay_database(data_source)
    decay_constants = decay_database.decay_constants_data
    matrix_P = sparse.load_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P.npz'))
    matrix_P_inv = sparse.load_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P_inv.npz'))
    return DecayMatrix(
        decay_constants,
        matrix_P,
        matrix_P_inv
    )


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


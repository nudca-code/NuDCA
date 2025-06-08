# -*- coding: utf-8 -*-

"""
Radioactive decay simulation module for nuclear astrophysics calculations.
Provides tools for modeling radioactive decay processes, calculating heating rates,
and managing nuclide abundances in nuclear decay networks.
"""

from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from numba import njit
import matplotlib.pyplot as plt

from .nuclide import Nuclide, NuclideStrError
from .decay_network import DecayDatabase, DecayMatrix
from .constants import NA_CGS, EV_CGS


class RadioactiveDecay:
    """
    Models radioactive decay processes and heating rates in nuclear decay networks.
    
    Provides methods for:
    - Computing time evolution of nuclide abundances
    - Calculating heating rates from radioactive decay
    - Managing and importing/exporting nuclide abundance data
    - Visualizing decay heating rates
    
    Uses a matrix-based approach for efficient computation of large decay networks.
    
    Attributes:
        initial_abundance (Dict[str, float]): Initial abundances of nuclides.
        decay_database (DecayDatabase): Database containing decay properties.
        decay_matrix (DecayMatrix): Matrix representation of decay network.
    """
    def __init__(
        self,
        initial_abundance: Dict[str, float],
        decay_database: DecayDatabase,
        decay_matrix: DecayMatrix
    ) -> None:
        self._initial_abundance = initial_abundance
        self.decay_database = decay_database
        self.decay_matrix = decay_matrix
        self._nuclides = None

    @property
    def initial_abundance(self) -> Dict[str, float]:
        """
        Get the initial nuclide abundances.
        Returns:
            Dict[str, float]: Mapping of nuclide symbol to abundance.
        """
        return self._initial_abundance

    @initial_abundance.setter
    def initial_abundance(self, value: Dict[str, float]):
        """
        Set the initial nuclide abundances.
        Args:
            value (Dict[str, float]): Mapping of nuclide symbol to abundance.
        """
        self._initial_abundance = value
        self._nuclides = None  # Invalidate cache when abundance changes

    @property
    def decay_database(self) -> DecayDatabase:
        """
        Get the decay database instance.
        Returns:
            DecayDatabase: The decay database.
        """
        return self._decay_database
    
    @decay_database.setter
    def decay_database(self, value: DecayDatabase):
        """
        Set the decay database instance.
        Args:
            value (DecayDatabase): The decay database.
        """
        self._decay_database = value

    @property
    def decay_matrix(self) -> DecayMatrix:
        """
        Get the decay matrix instance.
        Returns:
            DecayMatrix: The decay matrix.
        """
        return self._decay_matrix
    
    @decay_matrix.setter
    def decay_matrix(self, value: DecayMatrix):
        """
        Set the decay matrix instance.
        Args:
            value (DecayMatrix): The decay matrix.
        """
        self._decay_matrix = value


    @property
    def nuclides(self) -> List[str]:
        """
        Get the list of nuclide symbols in the initial abundance.
        Returns:
            List[str]: List of nuclide symbols.
        """
        if self._nuclides is None:
            self._nuclides = [
                Nuclide(nuclide).nuclide_symbol for nuclide in self.initial_abundance
            ]
        return self._nuclides
    

    def decay_nuclide_abundances(
        self,
        decay_times: Union[float, np.ndarray]
    ) -> Tuple[List[str], np.ndarray]:
        """
        Compute the time evolution of nuclide abundances in the decay network.
        
        Args:
            decay_times (float or np.ndarray): Single time point or array of time points to compute abundances.
        Returns:
            Tuple[List[str], np.ndarray]:
                - List of nuclide names in the network
                - 2D array of abundances (time x nuclides)
        Note:
            Efficiently computes abundances for all nuclides at all time points using matrix operations.
        """
        nuclide_names, abundances = self.decay_process(decay_times)
        return nuclide_names, abundances

    
    def decay_heating_rates(
        self,
        decay_times: Union[float, np.ndarray],
        energy_type: str = None
    ) -> np.ndarray:
        """
        Calculate heating rates from radioactive decay at specified time points.
        
        Args:
            decay_times (float or np.ndarray): Single time point or array of time points.
            energy_type (str, optional): Type of decay energy to use (e.g., 'beta', 'alpha', 'gamma').
                If None, uses total decay energy.
        Returns:
            np.ndarray: Array of heating rates at each time point.
        Note:
            Heating rates are calculated using the Bateman equations and include
            contributions from all relevant decay channels.
        """
        nuclide_names, abundances = self.decay_process(decay_times)
        indices = np.array([self.decay_database.nuclide_index_map[n] for n in nuclide_names])
        return self._heating_rates(abundances, indices, energy_type)
    

    def decay_process(
        self,
        decay_times: Union[float, np.ndarray]
    ) -> Tuple[List[str], np.ndarray]:
        """
        Compute nuclide abundance evolution using matrix operations (Bateman equations).
        
        Handles:
        - Initial condition setup
        - Matrix exponentiation for time evolution
        - Proper handling of decay chains
        Args:
            decay_times (float or np.ndarray): Time points for abundance calculation.
        Returns:
            Tuple[List[str], np.ndarray]:
                - List of nuclide names
                - 2D array of abundances (time x nuclides)
        Raises:
            ValueError: If decay times are negative or initial conditions are invalid.
            RuntimeError: If matrix operations fail during calculation.
        """
        decay_times = np.atleast_1d(decay_times)
        if np.any(decay_times < 0):
            raise ValueError("Decay time must be non-negative.")
        if not self.initial_abundance:
            raise ValueError("nuclide_abundance is empty.")
        if self.decay_database is None or self.decay_matrix is None:
            raise ValueError("decay_database and decay_matrix must not be None.")
        
        initial_abundance = self.decay_matrix.initial_abundance.copy()
        relevant_indices_set = set()
        for nuclide, quantity in self.initial_abundance.items():
            nuclide = self._validate_nuclide(nuclide)
            index = self.decay_database.nuclide_index_map[nuclide]
            initial_abundance[index] = quantity
            relevant_indices_set.add(index)
            relevant_indices_set.update(self.decay_matrix.matrix_P[:, index].nonzero()[0])
        relevant_indices = np.array(sorted(list(relevant_indices_set)), dtype=int)
        matrix_P = self.decay_matrix.matrix_P
        matrix_P_inv = self.decay_matrix.matrix_P_inv
        decay_constants = self.decay_matrix.decay_constants
        exp_factors = np.exp(-np.outer(decay_times, decay_constants[relevant_indices]))  # (n_times, n_indices)
        matrix_Lambda_base = self.decay_matrix.matrix_Lambda.copy()
        abundance_array = np.zeros((len(decay_times), len(relevant_indices)))
        for i, exp_row in enumerate(exp_factors):
            matrix_Lambda = matrix_Lambda_base.copy()
            matrix_Lambda.data[relevant_indices] = exp_row
            try:
                abundance_array[i] = (
                    matrix_P
                    .dot(matrix_Lambda)
                    .dot(matrix_P_inv)
                    .dot(initial_abundance)
                )[relevant_indices]
            except Exception as e:
                raise RuntimeError(
                    f"Error in abundance calculation at time index {i}: {e}\n"
                    f"Input: decay_times={decay_times}, relevant_indices={relevant_indices}"
                )
            
        return (
            [self.decay_database.nuclides[index] for index in relevant_indices],
            abundance_array
        )


    def plot_heating_rates(
        self,
        decay_times: Union[float, np.ndarray],
        legend: bool = True,
        save_path: str = None,
        xlim: Tuple[float, float] = None,
        ylim: Tuple[float, float] = None,
        energy_type: str = None
    ) -> None:
        """
        Generate a log-log plot of heating rates as a function of decay time.
        
        Args:
            decay_times (float or np.ndarray): Time points for plotting.
            legend (bool): Whether to display the plot legend.
            save_path (str, optional): Path to save the figure (if None, only displays).
            xlim (Tuple[float, float], optional): (xmin, xmax) for x-axis limits.
            ylim (Tuple[float, float], optional): (ymin, ymax) for y-axis limits.
            energy_type (str, optional): Type of decay energy to plot.
        Note:
            The plot uses log-log scaling for better visualization of the
            power-law behavior of radioactive decay heating.
        """
        heating_rates = self.decay_heating_rates(decay_times, energy_type=energy_type)
        x_min, x_max = min(decay_times), max(decay_times)
        y_min, y_max = min(heating_rates), max(heating_rates)
        fig, ax = plt.subplots()
        ax.plot(decay_times, heating_rates, label="Total Heating Rate" if energy_type is None else f"Heating Rate ({energy_type})")
        if xlim:
            ax.set_xlim(*xlim)
        else:
            ax.set_xlim(0.1 * x_min, 2 * x_max)
        if ylim:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(0.1 * y_min, 2 * y_max)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("Decay Time")
        ax.set_ylabel("Heating Rate")
        ax.set_title("Radioactive Decay Heating Rates")
        if legend:
            ax.legend()
        if save_path:
            plt.savefig(save_path)
        plt.show()

    def _validate_nuclide(self, nuclide: str) -> str:
        """
        Validate and standardize a nuclide string.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            str: Standardized nuclide symbol.
        Raises:
            NuclideStrError: If the nuclide string is invalid or not found in the database.
        """
        try:
            symbol = Nuclide(nuclide).nuclide_symbol
        except Exception as e:
            raise NuclideStrError(nuclide, f"Invalid format: {str(e)} | Input: {nuclide}")
        if symbol not in self.decay_database.nuclides:
            raise NuclideStrError(nuclide, f"Not found in decay database. Available: {self.decay_database.nuclides}. Input: {nuclide}")
        return symbol
    

    def _heating_rates(
        self,
        abundances: np.ndarray,
        indices: np.ndarray,
        energy_type: str = None
    ) -> np.ndarray:
        """
        Calculate heating rates from nuclide abundances.
        
        Uses Numba-accelerated computation for efficient heating rate calculation across multiple time points.
        Args:
            abundances (np.ndarray): Array of nuclide abundances (time x nuclides).
            indices (np.ndarray): Indices of relevant nuclides in the decay network.
            energy_type (str, optional): Type of decay energy to use.
        Returns:
            np.ndarray: Array of heating rates at each time point.
        Note:
            The calculation includes Avogadro's number and unit conversions
            to ensure proper physical units in the output.
        """
        decay_constants = self.decay_matrix.decay_constants[indices]
        if energy_type is None:
            decay_energies = np.array([
                self.decay_database.get_nuclide_total_decay_energy(self.decay_database.nuclides[index]) for index in indices
            ])
        else:
            decay_energies = np.array([
                self.decay_database.get_nuclide_decay_energy(self.decay_database.nuclides[index], energy_type) for index in indices
            ])
        decay_constants = np.asarray(decay_constants, dtype=np.float64)
        abundances = np.asarray(abundances, dtype=np.float64)
        decay_energies = np.asarray(decay_energies, dtype=np.float64)
        return calculate_radioactive_heating_rates(decay_constants, abundances, decay_energies)

 

@njit(cache=True)
def calculate_radioactive_heating_rates(
    decay_constants: np.ndarray,
    abundances: np.ndarray,
    decay_energies: np.ndarray
) -> np.ndarray:
    """
    Calculate heating rates from radioactive decay using Numba-accelerated computation.
    
    Implements the core heating rate calculation for radioactive decay processes.
    Computes the total heating rate by summing contributions from all nuclides at each time point,
    taking into account decay constants, abundances, and decay energies.
    
    Formula:
        Heating Rate = NA * λ * Y * E
    where:
        NA = Avogadro's number
        λ = decay constant
        Y = nuclide abundance
        E = decay energy
    Args:
        decay_constants (np.ndarray): Array of decay constants for each nuclide [s^-1].
        abundances (np.ndarray): 2D array of nuclide abundances [mol], shape (n_times, n_nuclides).
        decay_energies (np.ndarray): Array of decay energies for each nuclide [eV].
    Returns:
        np.ndarray: Array of heating rates at each time point [erg/s].
    Notes:
        - Uses Numba JIT compilation for performance optimization.
        - Invalid values (NaN, inf) are automatically skipped in the calculation.
        - Result includes unit conversions from eV to erg using EV_CGS constant.
        - Assumes all input arrays are properly aligned.
    """
    n_times = abundances.shape[0]
    n_nuclides = abundances.shape[1]
    heating_rates = np.zeros(n_times)
    for t in range(n_times):
        for i in range(n_nuclides):
            if (
                not np.isfinite(decay_constants[i]) or
                np.isnan(abundances[t, i]) or
                np.isnan(decay_energies[i])
            ):
                continue
            heating_rates[t] += NA_CGS * decay_constants[i] * abundances[t, i] * decay_energies[i] * EV_CGS
    return heating_rates





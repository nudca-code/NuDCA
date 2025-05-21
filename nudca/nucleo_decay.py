# -*- coding: utf-8 -*-

from typing import Dict, Set, List, Union

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

from .nuclide import Nuclide, NuclideStrError
from .decay_database import DecayDatabase
from .decay_matrix import DecayMatrix
from .constants import NA_CGS, EV_CGS


class RadioactiveDecayBase:
    """
    Base class for modeling radioactive decay processes and heating rates for a collection of nuclides.
    Provides methods for abundance evolution, heating rate calculation, and validation utilities.
    """

    def __init__(
        self,
        nuclide_abundance: Dict[str, int],
        decay_database: DecayDatabase,
        decay_matrix: DecayMatrix
    ) -> None:
        """
        Initialize the base decay model with nuclide abundances, decay database, and decay matrix.
        Args:
            nuclide_abundance (Dict[str, int]): Initial abundances for each nuclide.
            decay_database (DecayDatabase): Database with nuclide decay data.
            decay_matrix (DecayMatrix): Precomputed decay matrix for the nuclide set.
        """
        self.nuclide_abundance = nuclide_abundance
        self.decay_database = decay_database
        self.decay_matrix = decay_matrix
        self._nuclides = None
        

    @property
    def nuclides(self) -> List[str]:
        """
        Retrieve and cache the list of nuclide symbols present in the abundance dictionary.
        Returns:
            List[str]: List of nuclide symbols.
        """
        if self._nuclides is None:
            self._nuclides = [
                Nuclide(nuclide).nuclide_symbol for nuclide in self.nuclide_abundance
            ]
        return self._nuclides


    @property
    def abundance(self, nuclide) -> float:
        """
        Return the abundance of a specific nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Abundance of the nuclide.
        """
        nuclide = self._validate_nuclide(nuclide)
        return self.nuclide_abundance[nuclide]


    def decay_process(self, decay_time: float):
        """
        Calculate the decay of nuclides over a specified time period using the Bateman solution.
        Args:
            decay_time (float): Time period for decay calculation in seconds.
        Returns:
            A new instance of the same class with the decayed abundances.
        Raises:
            ValueError: If decay_time is negative.
        """
        if decay_time < 0:
            raise ValueError("Decay time must be non-negative.")
        # Prepare initial abundance vector and diagonal decay matrix
        N0 = self.decay_matrix.initial_abundance.copy()
        matrix_Lambda = self.decay_matrix.matrix_Lambda.copy()
        # Collect indices of relevant nuclides
        indices_set = set()
        for nuclide, quantity in self.nuclide_abundance.items():
            nuclide = self._validate_nuclide(nuclide)
            index = self.decay_database.nuclide_index_map[nuclide]
            N0[index] = quantity
            # Get all indices that might be affected by this nuclide's decay
            indices_set.update(self.decay_matrix.matrix_P[:, index].nonzero()[0])
        indices = np.array(list(indices_set), dtype=int)
        # Exponential decay for each relevant nuclide
        matrix_Lambda.data[indices] = np.exp(
            -decay_time * self.decay_matrix.decay_constants[indices]
        )
        # Bateman solution: P * Lambda * P_inv * N0
        decayed_abundance_array = (
            self.decay_matrix.matrix_P
            .dot(matrix_Lambda)
            .dot(self.decay_matrix.matrix_P_inv)
            .dot(N0)
        )
        decayed_abundance = dict(
            zip(
                self.decay_database.nuclides[indices],
                decayed_abundance_array[indices]
            )
        )
        return self.__class__(decayed_abundance, self.decay_database, self.decay_matrix)


    @property
    def heating_rate(self) -> float:
        """
        Calculate the total radioactive heating rate (all nuclides, all decay types).
        Returns:
            float: Total heating rate in erg/s.
        """
        return np.sum([
            self._heating_rate_by_nuclide(nuclide)
            for nuclide in self.nuclide_abundance
        ])
    
    
    def _heating_rate_by_nuclide_and_type(self, nuclide: str, energy_type: str) -> float:
        """
        Calculate the heating rate for a specific nuclide and decay energy type.
        Args:
            nuclide (str): Nuclide symbol.
            energy_type (str): Type of decay energy (e.g. 'EM', 'LP', ...).
        Returns:
            float: Heating rate in erg/s for the given nuclide and energy type.
        """
        nuclide = self._validate_nuclide(nuclide)
        index = self.decay_database.nuclide_index_map[nuclide]
        decay_constant = self.decay_matrix.decay_constants[index]
        decay_abundance = self.nuclide_abundance[nuclide]
        decay_energy = self.decay_database.decay_energy(nuclide, energy_type)
        if (
            not np.isfinite(decay_constant) or
            np.isnan(decay_constant) or
            np.isnan(decay_abundance) or 
            np.isnan(decay_energy) 
        ):
            return 0.0
        # Heating rate: N_A * lambda * N * Q * eV->erg
        heating_rate = NA_CGS * decay_constant * decay_abundance * decay_energy * EV_CGS
        return heating_rate
        

    def _heating_rate_by_nuclide(
        self,
        nuclide: str,
    ) -> float:
        """
        Calculate the total heating rate for a specific nuclide (all decay types summed).
        Args:
            nuclide (str): Nuclide symbol.
        Returns:
            float: Heating rate in erg/s for the given nuclide.
        """
        nuclide = self._validate_nuclide(nuclide)
        index = self.decay_database.nuclide_index_map[nuclide]
        decay_constant = self.decay_matrix.decay_constants[index]
        decay_abundance = self.nuclide_abundance[nuclide]
        decay_energy = self.decay_database.decay_energy_released(nuclide)
        if (
            not np.isfinite(decay_constant) or
            np.isnan(decay_constant) or
            np.isnan(decay_abundance) or 
            np.isnan(decay_energy) 
        ):
            return 0.0
        heating_rate = NA_CGS * decay_constant * decay_abundance * decay_energy * EV_CGS
        return heating_rate
        
    def _heating_rate_by_type(self, energy_type: str) -> float:
        """
        Calculate the total heating rate for a specific decay energy type (summed over all nuclides).
        Args:
            energy_type (str): Type of decay energy (e.g. 'EM', 'LP', ...).
        Returns:
            float: Total heating rate in erg/s for the given energy type.
        """
        return np.sum([
            self._heating_rate_by_nuclide_and_type(nuclide, energy_type)
            for nuclide in self.nuclide_abundance
        ])

    def _validate_nuclide(self, nuclide: str) -> str:
        """
        Validate nuclide and return its standardized symbol.
        Args:
            nuclide (str): Nuclide string or symbol.
        Returns:
            str: Standardized nuclide symbol.
        Raises:
            NuclideStrError: If nuclide is invalid or not found in the database.
        """
        try:
            symbol = Nuclide(nuclide).nuclide_symbol
        except Exception as e:
            raise NuclideStrError(nuclide, f"Invalid format: {str(e)}")
        if symbol not in self.decay_database.nuclides:
            raise NuclideStrError(nuclide, f"Not found in decay database.")
        return symbol


class RadioactiveDecay(RadioactiveDecayBase):
    """
    Main class for time-dependent radioactive decay and heating rate calculations.
    Provides time series, parallelized batch calculations, and plotting utilities.
    """
    def __init__(self, 
                 nuclide_abundance: Dict[str, int],
                 decay_database: DecayDatabase,
                 decay_matrix: DecayMatrix
                ) -> None:
        """
        Initialize the RadioactiveDecay model.
        Args:
            nuclide_abundance (Dict[str, int]): Initial abundances for each nuclide.
            decay_database (DecayDatabase): Database with nuclide decay data.
            decay_matrix (DecayMatrix): Precomputed decay matrix for the nuclide set.
        """
        super().__init__(nuclide_abundance, decay_database, decay_matrix)
        
    def decay_time_series_pandas(
        self,
        time_period: Union[float, np.ndarray],
        time_units: str = "s",
        time_scale: str = "linear",
        npoints: int = 501,
    ) -> pd.DataFrame:
        """
        Calculate decay time series and return as a pandas DataFrame.
        Args:
            time_period (float or np.ndarray): Either a single float value or an array of time points.
            time_units (str): Units of time (default: "s").
            time_scale (str): Scale of time points ("linear" or "log").
            npoints (int): Number of points for automatic time series generation.
        Returns:
            pd.DataFrame: DataFrame with time series data (abundances for each nuclide).
        """
        if isinstance(time_period, (int, float)):
            tmin = 0.0 if time_scale == "linear" else 0.1
            time_points = (
                np.linspace(tmin, time_period, num=npoints)
                if time_scale == "linear"
                else np.logspace(np.log10(tmin), np.log10(time_period), num=npoints)
            )
        else:
            time_points = np.asarray(time_period)
            if np.any(time_points < 0):
                raise ValueError("All time points must be non-negative")
        # Use parallel processing for decay calculations
        decayed_states = Parallel(n_jobs=-1)(
            delayed(self.decay_process)(t) for t in time_points
        )
        # Collect data from all decayed states
        time_column = f"Time ({time_units})"
        data = {time_column: time_points}
        # Get all unique nuclides from all states
        all_nuclides = set()
        for state in decayed_states:
            all_nuclides.update(state.nuclide_abundance.keys())
        # Collect abundances for each nuclide
        for nuclide in sorted(all_nuclides):
            data[nuclide] = [
                state.nuclide_abundance.get(nuclide, 0.0) 
                for state in decayed_states
            ]
        return pd.DataFrame(data).set_index(time_column)

    def decay_time_series(
        self,
        time_period: Union[float, np.ndarray],
        time_units: str = "s",
        time_scale: str = "linear",
        npoints: int = 501,
    ) -> tuple[list[float], dict[str, list[float]]]:
        """
        Calculate decay time series and return as tuple of time points and abundances.
        Args:
            time_period (float or np.ndarray): Either a single float value or an array of time points.
            time_units (str): Units of time (default: "s").
            time_scale (str): Scale of time points ("linear" or "log").
            npoints (int): Number of points for automatic time series generation.
        Returns:
            tuple: (list of time points, dict mapping nuclide symbols to lists of abundances)
        """
        df = self.decay_time_series_pandas(
            time_period=time_period,
            time_units=time_units,
            time_scale=time_scale,
            npoints=npoints,
        )
        return (list(df.index), df.to_dict(orient="list"))
    
    def _calc_heating_rate_by_type(self, decay_time, energy_type):
        """
        Helper for parallelized heating rate calculation by energy type at a given time.
        """
        return self.decay_process(decay_time)._heating_rate_by_type(energy_type)
    
    def _calc_heating_rate(self, decay_time: float) -> float:
        """
        Helper for parallelized total heating rate calculation at a given time.
        """
        return self.decay_process(decay_time).heating_rate  
    
    def heating_rates_by_type(self, decay_times, energy_type):
        """
        Calculate heating rates for a sequence of times and a specific energy type (parallelized).
        Args:
            decay_times (array-like): Sequence of decay times.
            energy_type (str): Decay energy type.
        Returns:
            list: Heating rates at each time.
        """
        heating_rates = Parallel(
            n_jobs=-1,
            # backend='threading'
        )(
            delayed(self._calc_heating_rate_by_type)(time, energy_type) for time in decay_times
        )
        return heating_rates

    def decay_heating_rates(self, decay_times) -> List[float]:
        """
        Calculate total heating rates for a sequence of times (parallelized).
        Args:
            decay_times (array-like): Sequence of decay times.
        Returns:
            list: Total heating rates at each time.
        Raises:
            ValueError: If any decay time is negative.
        """
        if np.any(decay_times < 0):
            raise ValueError("All decay times must be non-negative.")
        total_heating_rates = Parallel(
            n_jobs=-1,
            # backend='threading'
        )(
            delayed(self._calc_heating_rate)(time) for time in decay_times
        )
        return total_heating_rates

    def plot_heating_rates(self, decay_times):
        """
        Plot the total heating rate as a function of decay time (log-log scale).
        Args:
            decay_times (array-like): Sequence of decay times.
        """
        heating_rates = self.decay_heating_rates(decay_times)
        x_min, x_max = min(decay_times), max(decay_times)
        y_min, y_max = min(heating_rates), max(heating_rates)
        fig, ax = plt.subplots()
        ax.plot(decay_times, heating_rates)
        ax.set_xlim(0.1*x_min, 2*x_max)
        ax.set_ylim(0.1*y_min, 2*y_max)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.show()
        # return fig



    
    
    
       

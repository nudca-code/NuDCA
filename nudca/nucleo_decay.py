# -*- coding: utf-8 -*-

from typing import Dict, Set, List, Union

import numpy as np
import pandas as pd
from scipy import sparse
from joblib import Parallel, delayed

import matplotlib.pyplot as plt

from .nuclide import Nuclide, NuclideStrError
from .decay_database import DecayDatabase
from .decay_matrix import DecayMatrix
from .constants import NA_CGS, EV_CGS


class RadioactiveDecayBase:

    def __init__(
        self,
        nuclide_abundance: Dict[str, int],
        decay_database: DecayDatabase,
        decay_matrix: DecayMatrix
    ) -> None:
        self.nuclide_abundance = nuclide_abundance
        self.decay_database = decay_database
        self.decay_matrix = decay_matrix

        self._nuclides = None
        

    @property
    def nuclides(self) -> List[str]:
        """Retrieve and cache the list of nuclide symbols."""
        if self._nuclides is None:
            self._nuclides = [
                Nuclide(nuclide).nuclide_symbol for nuclide in self.nuclide_abundance
            ]
        return self._nuclides


    @property
    def abundance(self, nuclide) -> float:
        """Return the abundances of nuclides."""
        nuclide = self._validate_nuclide(nuclide)
        return self.nuclide_abundance[nuclide]


    def decay_process(self, decay_time: float):
        """
        Calculate the decay of nuclides over a specified time period.
        
        Args:
            decay_time (float): Time period for decay calculation in seconds.
            
        Returns:
            A new instance with the decayed abundances.
            
        Raises:
            ValueError: If decay_time is negative.
        """
        if decay_time < 0:
            raise ValueError("Decay time must be non-negative.")
            
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
        
        matrix_Lambda.data[indices] = np.exp(
            -decay_time * self.decay_matrix.decay_constants[indices]
        )
        
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
        return np.sum([
            self._heating_rate_by_nuclide(nuclide)
            for nuclide in self.nuclide_abundance
        ])
    
    
    def _heating_rate_by_nuclide_and_type(self, nuclide: str, energy_type: str) -> float:
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
        
        heating_rate = NA_CGS * decay_constant * decay_abundance * decay_energy * EV_CGS
        
        return heating_rate
        


    def _heating_rate_by_nuclide(
        self,
        nuclide: str,
    ) -> float:
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
        return np.sum([
            self._heating_rate_by_nuclide_and_type(nuclide, energy_type)
            for nuclide in self.nuclide_abundance
        ])



    def _validate_nuclide(self, nuclide: str) -> str:
        """Validate nuclide and return its standardized symbol."""
        try:
            symbol = Nuclide(nuclide).nuclide_symbol
        except Exception as e:
            raise NuclideStrError(nuclide, f"Invalid format: {str(e)}")
        if symbol not in self.decay_database.nuclides:
            raise NuclideStrError(nuclide, f"Not found in decay database.")
        return symbol


    
    




class RadioactiveDecay(RadioactiveDecayBase):
    def __init__(self, 
                 nuclide_abundance: Dict[str, int],
                 decay_database: DecayDatabase,
                 decay_matrix: DecayMatrix
                ) -> None:
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
            time_period: Either a single float value or an array of time points
            time_units: Units of time (default: "s")
            time_scale: Scale of time points ("linear" or "log")
            npoints: Number of points for automatic time series generation
            
        Returns:
            DataFrame with time series data
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
            time_period: Either a single float value or an array of time points
            time_units: Units of time (default: "s")
            time_scale: Scale of time points ("linear" or "log")
            npoints: Number of points for automatic time series generation
            
        Returns:
            Tuple containing:
            - List of time points
            - Dictionary mapping nuclide symbols to lists of abundances
        """
        df = self.decay_time_series_pandas(
            time_period=time_period,
            time_units=time_units,
            time_scale=time_scale,
            npoints=npoints,
        )

        return (list(df.index), df.to_dict(orient="list"))
    
    

    def _calc_heating_rate_by_type(self, decay_time, energy_type):
        return self.decay_process(decay_time)._heating_rate_by_type(energy_type)
    
    
    def _calc_heating_rate(self, decay_time: float) -> float:
        return self.decay_process(decay_time).heating_rate  
    

    def heating_rates_by_type(self, decay_times, energy_type):
        heating_rates = Parallel(
            n_jobs=-1,
            # backend='threading'
        )(
            delayed(self._calc_heating_rate_by_type)(time, energy_type) for time in decay_times
        )
        return heating_rates



    def decay_heating_rates(self, decay_times) -> List[float]:
        if np.any(decay_times < 0):
            raise ValueError("All decay times must be non-negative.")
        
        total_heating_rates = Parallel(
            n_jobs=-1,
            # backend='threading'
        )(
            delayed(self._calc_heating_rate)(time) for time in decay_times

        )

        return total_heating_rates




    # def heating_rates_by_nuclide_and_type(self, decay_times, nuclide, energy_type):
    #     def calc_heating_rate_by_nuclide_and_type(decay_time):
    #         return self.decay_process(decay_time)._heating_rate_by_nuclide_and_type(nuclide, energy_type)
        
    #     heating_rates = Parallel(n_jobs=-1)(
    #         delayed(calc_heating_rate_by_nuclide_and_type)(time) for time in decay_times
    #     )

    #     return heating_rates



    # def heating_rates_by_nuclide(self, decay_times, nuclide):
        
    #     def calc_heating_rate_by_nuclide(decay_time):
    #         return self.decay_process(decay_time)._heating_rate_by_nuclide(nuclide)
        
    #     heating_rates = Parallel(n_jobs=-1)(
    #         delayed(calc_heating_rate_by_nuclide)(time) for time in decay_times
    #     )

    #     return heating_rates
    
    
    
    # def heating_rates_by_type(self, decay_times, energy_type):
        
    #     def calc_heating_rate_by_type(decay_time):
    #         return self.decay_process(decay_time)._heating_rate_by_type(energy_type)
        
    #     heating_rates = Parallel(n_jobs=-1)(
    #         delayed(calc_heating_rate_by_type)(time) for time in decay_times
    #     )
        
    #     return heating_rates
        
    
    def plot_heating_rates(self, decay_times):
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
        
        
    



    
    
    
       

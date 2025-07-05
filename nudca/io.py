# -*- coding: utf-8 -*-

"""
This module contains the Inputer, Outputer, and AbundanceEstimator classes for reading and writing nuclide abundances and decay data.

load_decay_database:
    - Loads a DecayDatabase object from a compressed NumPy file for the specified data source.
load_decay_matrix:
    - Loads a DecayMatrix object from a compressed NumPy file for the specified data source.
serialize_decay_database:
    - Saves a DecayDatabase object to a compressed NumPy file for the specified data source.
serialize_decay_matrix:
    - Saves a DecayMatrix object to a compressed NumPy file for the specified data source.

Inputer:
    - Reads nuclide abundances from file.
    - Supports CSV, TSV, Excel, and JSON formats.

Outputer:
    - Saves nuclide abundances to file.

AbundanceEstimator:
    - Estimates the initial abundances of parent nuclides based on decay chains and final abundances.
"""

import re
from pathlib import Path
from collections import deque
from typing import Any, Optional, Union, Dict, List, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sparse

from .decay_network import (
    DecayDatabase,
    DecayMatrix,
    MatrixBuilder
)
from .utils import DecayDatabaseManager
from .nuclide import Nuclide, Z_to_element


def load_decay_database(data_source: str = 'ENDF-B-VIII.1_decay') -> DecayDatabase:
    """
    Load a DecayDatabase object from a compressed NumPy file for the specified data source.
    Args:
        data_source (str): Name of the decay data source.
    Returns:
        DecayDatabase: Loaded decay database object.
    Raises:
        FileNotFoundError: If the data file does not exist.
        KeyError: If required fields are missing in the data file.
    """
    try:
        database = np.load(
            DecayDatabaseManager(data_source).data_path.joinpath(f'{data_source}.npz'),
            allow_pickle=True,
        )
    except FileNotFoundError:
        raise FileNotFoundError(f"Decay database file for {data_source} not found.")
    except KeyError as e:
        raise KeyError(f"Required field {e} not found in decay database file for {data_source}.")
    nuclides = database['nuclides']
    half_life_data = database['half_life_data']
    decay_constants_data = database['decay_constants_data']
    decay_modes_data = database['decay_modes_data']
    progeny_data = database['progeny_data']
    branching_ratios_data = database['branching_ratios_data']
    decay_energies_data = database['decay_energies_data']

    return DecayDatabase(
        data_source,
        nuclides,
        half_life_data,
        decay_constants_data,
        decay_modes_data,
        progeny_data,
        branching_ratios_data,
        decay_energies_data
    )


def load_decay_matrix(
        data_source: str='ENDF-B-VIII.1_decay'
    ) -> DecayMatrix:
    """
    Load the decay matrix P and its inverse from disk and return a DecayMatrix object.
    Args:
        data_source (str): Name of the decay data source.
    Returns:
        DecayMatrix: The loaded decay matrix object.
    Raises:
        FileNotFoundError: If the matrix files do not exist.
    """
    decay_database = load_decay_database(data_source)
    decay_constants = decay_database.decay_constants_data
    try:
        matrix_P = sparse.load_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P.npz'))
        matrix_P_inv = sparse.load_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P_inv.npz'))
    except FileNotFoundError:
        raise FileNotFoundError("Decay matrix files not found.")

    return DecayMatrix(
        decay_constants,
        matrix_P,
        matrix_P_inv
    )


def serialize_decay_database(data_source: str='ENDF-B-VIII.1_decay') -> None:
    """
    Load, sort, and save the decay database in both CSV and compressed NumPy formats.
    Prepares the data for fast loading and use in DecayDatabase.
    """
    try:
        df = pd.read_json(Path(__file__).resolve().parent.joinpath('data', f'{data_source}.json'), orient='records')
        df = DecayDatabaseManager(data_source).sort_nuclides_data(df)
        df.to_csv(Path(__file__).resolve().parent.joinpath('data', f'{data_source}_sorted.csv'), index=False)
        df.to_json(Path(__file__).resolve().parent.joinpath('data', f'{data_source}_sorted.json'), orient='records', indent=4)
    except FileNotFoundError:
        raise FileNotFoundError(f"Decay database file for {data_source} not found.")
    except KeyError as e:
        raise KeyError(f"Required field {e} not found in decay database file for {data_source}.")
    (
        nuclides,
        half_life_data,
        decay_constants_data,
        decay_modes_data,
        progeny_data,
        branching_ratios_data,
        decay_energies_data
    ) = DecayDatabaseManager(data_source).transform_radionuclide_data(df)
    np.savez_compressed(
        Path(__file__).resolve().parent.joinpath('data', f'{data_source}.npz'),
        nuclides=np.asarray(nuclides, dtype=object),
        half_life_data=np.asarray(half_life_data, dtype=object),
        decay_constants_data=np.asarray(decay_constants_data, dtype=np.float64),
        decay_modes_data=np.asarray(decay_modes_data, dtype=object),
        progeny_data=np.asarray(progeny_data, dtype=object),
        branching_ratios_data=np.asarray(branching_ratios_data, dtype=object),
            decay_energies_data=np.asarray(decay_energies_data, dtype=np.float64) 
        )
    


def serialize_decay_matrix(data_source: str='ENDF-B-VIII.1_decay'):
    """
    Save the decay matrix P and its inverse to disk as .npz files for a given data source.
    Args:
        data_source (str): Name of the decay data source.
    Returns:
        None
    Raises:
        FileNotFoundError: If the decay database file does not exist.
    """
    try:
        decay_database = load_decay_database(data_source)
        matrix_P, matrix_P_inv = MatrixBuilder(decay_database).build_decay_matrix()
        # Save the matrices in compressed sparse format for efficient storage
        sparse.save_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P.npz'), matrix_P.tocsr())
        sparse.save_npz(Path(__file__).resolve().parent.joinpath('data', 'decay_matrix_P_inv.npz'), matrix_P_inv.tocsr())
    except FileNotFoundError:
        raise FileNotFoundError(f"Decay database file for {data_source} not found.")
    except KeyError as e:
        raise KeyError(f"Required field {e} not found in decay database file for {data_source}.")



class AbundanceEstimator:
    """A class for estimating initial abundances based on decay chains and final abundances.

    This class helps estimate the initial abundances of parent nuclides based on
    known decay chains and the final abundance of a stable daughter nuclide.

    Attributes:
        decay_database: The database containing decay information
    """

    def __init__(self, decay_database: DecayDatabase):
        """Initialize the AbundanceEstimator.

        Args:
            decay_database: The database containing decay information
        """
        self.decay_database = decay_database
        
        
    def initial_abundances_rProcess(self, final_abundances: Dict[str, float]) -> Dict[str, float]:
        """Estimate the initial abundances of parent nuclides based on decay chains and final abundances.
        
        Args:
            nuclides: List of nuclides
            final_abundances: Dictionary mapping nuclides to their final abundances
        Returns:
            Dictionary mapping initial nuclides to their combined abundances
        """
        combined_initial_abundances = {}
        for nuclide, final_abundance in final_abundances.items():
            decay_chains = self.build_rProcess_chains(nuclide)
            chain_contributions = self._calculate_chain_contributions(decay_chains)
            for chain, contribution in chain_contributions.items():
                initial_nuclide = chain[0]
                abundance = final_abundance * contribution
                if initial_nuclide in combined_initial_abundances:
                    combined_initial_abundances[initial_nuclide] += abundance
                else:
                    combined_initial_abundances[initial_nuclide] = abundance
        return combined_initial_abundances
    
    
    def build_rProcess_chains(self, nuclide: str) -> List[List[str]]:
        nuclide = Nuclide(nuclide).nuclide_symbol
        decay_chains = []
        stack = deque([(nuclide, [nuclide])])

        while stack:
            current, path = stack.pop()  # pop() for DFS, popleft() for BFS
            parents = self.decay_database.get_nuclide_ancestor(current)
            valid_parents = []
            for parent in parents:
                decay_mode, _ = self.decay_database.get_decay_channel(parent, current)
                if decay_mode in ['β+&EC', 'SF']:
                    continue
                # if self.decay_database.get_nuclide_half_life(parent, 's') < 1.0:
                #     continue
                if (parent in self.decay_database.nuclide_index_map and 
                    np.isfinite(self.decay_database.get_nuclide_half_life(parent, 's'))):
                    valid_parents.append(parent)
            if not valid_parents:
                decay_chains.append(path[::-1])
            else:
                for parent in valid_parents:
                    stack.append((parent, path + [parent]))
        return decay_chains
    
    
    def _calculate_chain_contributions(self, decay_chains: List[List[str]]) -> Dict[Tuple[str, ...], float]:
        """Calculate the contribution ratios of each decay chain to the final nuclide.

        Args:
            decay_chains: List of decay chains

        Returns:
            Dict[Tuple[str, ...], float]: Dictionary mapping chains to their contribution ratios
        """
        chain_contributions = {}
        total_contribution = 0.0

        for chain in decay_chains:
            contribution = 1.0
            for i in range(len(chain) - 1):
                parent = chain[i]
                daughter = chain[i + 1]
                
                daughters = self.decay_database.get_nuclide_progeny(parent)
                daughter_index = daughters.index(daughter)
                branching_ratio = self.decay_database.get_nuclide_branching_ratios(parent)[daughter_index]
                
                contribution *= branching_ratio

            chain_contributions[tuple(chain)] = contribution
            total_contribution += contribution

        for chain in chain_contributions:
            chain_contributions[chain] /= total_contribution

        return chain_contributions


    def calculate_initial_abundances(
        self,
        decay_chains: List[List[str]],
        final_nuclide: str,
        final_abundance: float
    ) -> Dict[str, float]:
        """Calculate initial abundances based on decay chains and final abundance.

        Args:
            decay_chains: List of decay chains, each chain is a list of nuclides
            final_nuclide: The final stable nuclide
            final_abundance: The abundance of the final nuclide

        Returns:
            Dict[str, float]: Dictionary mapping initial nuclides to their abundances
        """

        chain_contributions = self._calculate_chain_contributions(decay_chains)

        initial_abundances = {}
        for chain, contribution in chain_contributions.items():
            initial_nuclide = chain[0]
            initial_abundances[initial_nuclide] = final_abundance * contribution

        return initial_abundances




class Inputer:
    """
    Input utility for reading initial nuclide abundances from file.
    Supports CSV, TSV, Excel, and JSON formats.
    """
    def __init__(self, decay_database: DecayDatabase):
        self.decay_database = decay_database

    @staticmethod
    def extract_mass_number(nuclide):
        match = re.search(r'\d+', nuclide)
        return int(match.group()) if match else float('inf')

    @staticmethod
    def to_dict(Z, A, Y):
        nuclide_Y = {}
        for z, a, y in zip(Z, A, Y):
            element = Z_to_element(z)
            nuclide = Nuclide(f'{element}{int(a)}').nuclide_symbol
            nuclide_Y[nuclide] = y

        return nuclide_Y

    def get_sorted_Y(self, nuclide_Y):
        sorted_Y = dict(sorted(nuclide_Y.items(), key=lambda x: self.extract_mass_number(x[0])))    
        A, Y = [], []
        for key, value in sorted_Y.items():
            A_decay = Nuclide(key).A
            Y_decay = value
            A.append(A_decay)
            Y.append(Y_decay)
        return A, Y

    def filter_Y(self, Z, A, Y):
        nuclide_dict = self.to_dict(Z, A, Y)
        nuclides = self.decay_database.nuclides
        common_nuclides = set(nuclide_dict.keys()) & set(nuclides)
        filtered_dict = {key: nuclide_dict[key] for key in common_nuclides}
        
        return filtered_dict
    
    
    def filter_nuclide_abundance(self, nuclide_abundance: Dict[str, float]) -> Dict[str, float]:
        nuclides = self.decay_database.nuclides
        common_nuclides = set(nuclide_abundance.keys()) & set(nuclides)
        filtered_dict = {key: nuclide_abundance[key] for key in common_nuclides}
        return filtered_dict




    @staticmethod
    def read_abundance(
        filename: Union[str, Path],
        filetype: Optional[str] = None,
        atomic_number_col: str = "Z",
        atomic_mass_col: str = "A",
        abundance_col: str = "Y"
    ) -> Dict[str, float]:
        """
        Read initial nuclide abundances from file.
        Args:
            filename: Input file path.
            filetype: 'csv', 'xlsx', or 'json'. If None, inferred from filename.
            atomic_number_col: Column name for atomic number (default 'Z').
            atomic_mass_col: Column name for atomic mass (default 'A').
            abundance_col: Column name for abundance (default 'Y').
        Returns:
            Dict[str, float]: Mapping from nuclide symbol to abundance.
        """
        if filetype is None:
            filetype = str(filename).split('.')[-1].lower()
        if filetype == 'csv':
            df = pd.read_csv(filename)
        elif filetype in ('xls', 'xlsx'):
            df = pd.read_excel(filename)
        elif filetype == 'json':
            df = pd.read_json(filename)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")
        # Try to find the correct columns
        # cols = [col.lower() for col in df.columns]
        Z_col = atomic_number_col if atomic_number_col in df.columns else None
        A_col = atomic_mass_col if atomic_mass_col in df.columns else None
        Y_col = abundance_col if abundance_col in df.columns else None
        if Z_col is None:
            for c in df.columns:
                if c.lower() in ["atomic number", "atomic_number", "z"]:
                    Z_col = c
                    break
        if Y_col is None:
            for c in df.columns:
                if c.lower() in ["abundance", "quantity", "amount", "y"]:
                    Y_col = c
                    break
        if Z_col is None or Y_col is None:
            raise ValueError(f"Could not find nuclide or abundance columns in {filename}")
        # Standardize nuclide symbols
        nuclide_abundance = {}
        for _, row in df.iterrows():
            symbol = Nuclide(row[Z_col]).nuclide_symbol
            nuclide_abundance[symbol] = float(row[Y_col])
        return nuclide_abundance
    

    @staticmethod
    def read_abundance_from_file(filename: Union[str, Path], filetype: Optional[str] = None) -> Dict[str, float]:
        """
        Read nuclide abundances from file.
        Args:
            filename: Input file path.
            filetype: 'csv', 'xlsx', or 'json'. If None, inferred from filename.
        Returns:
            Dict[str, float]: Mapping from nuclide symbol to abundance.
        """
        pass
            
            

    def import_nuclide_abundance(self, filepath: str, fmt: str = 'csv') -> None:
        """
        Import nuclide abundances from a file.
        
        Supports multiple file formats (CSV, JSON, Excel).
        The file must contain 'nuclide' and 'abundance' columns.
        Args:
            filepath (str): Path to the input file.
            fmt (str): File format ('csv', 'json', 'excel').
        Raises:
            FileNotFoundError: If the specified file doesn't exist.
            ValueError: If the file format is invalid or required columns are missing.
        """
        fmt = fmt.lower()
        if fmt == 'csv':
            df = pd.read_csv(filepath)
        elif fmt == 'json':
            df = pd.read_json(filepath)
        elif fmt in ('xlsx', 'excel'):
            df = pd.read_excel(filepath)
        else:
            raise ValueError(f"Unsupported format: {fmt}")
        if not {'nuclide', 'abundance'}.issubset(df.columns):
            raise ValueError(f"File must contain columns: 'nuclide', 'abundance'. Found: {df.columns}")
        abundance_dict = dict(zip(df['nuclide'], df['abundance']))
        abundance_dict = {str(k): float(v) for k, v in abundance_dict.items()}
        self.nuclide_abundance = abundance_dict



    def export_nuclide_abundance(self, filepath: str, fmt: str = 'csv') -> None:
        """
        Export current nuclide abundances to a file.
        
        Supports multiple file formats (CSV, JSON, Excel).
        The output file contains 'nuclide' and 'abundance' columns.
        Args:
            filepath (str): Path to save the output file.
            fmt (str): File format ('csv', 'json', 'excel').
        Raises:
            ValueError: If the specified format is not supported.
        """
        fmt = fmt.lower()
        df = pd.DataFrame({
            'nuclide': list(self.nuclide_abundance.keys()),
            'abundance': list(self.nuclide_abundance.values())
        })
        if fmt == 'csv':
            df.to_csv(filepath, index=False)
        elif fmt == 'json':
            df.to_json(filepath, orient='records', force_ascii=False)
        elif fmt in ('xlsx', 'excel'):
            df.to_excel(filepath, index=False)
        else:
            raise ValueError(f"Unsupported format: {fmt}")


    def export_decay_abundance(
        self,
        decay_times: Union[float, np.ndarray],
        filepath: str,
        fmt: str = 'csv',
        index_as_time: bool = True
    ) -> None:
        """
        Export time evolution of nuclide abundances to a file.
        
        Provides flexible output formats for decay calculations,
        supporting both wide (time as index) and long (melted) formats.
        Args:
            decay_times (float or np.ndarray): Time points for abundance calculation.
            filepath (str): Path to save the output file.
            fmt (str): File format ('csv', 'json', 'excel').
            index_as_time (bool): If True, uses time points as index; if False, uses long format.
        Raises:
            ValueError: If the specified format is not supported.
        """
        nuclide_names, abundances = self.decay_nuclide_abundances(decay_times)
        df = pd.DataFrame(abundances, columns=nuclide_names, index=np.atleast_1d(decay_times))
        df.index.name = 'time'
        if not index_as_time:
            df = df.reset_index().melt(id_vars='time', var_name='nuclide', value_name='abundance')
        fmt = fmt.lower()
        if fmt == 'csv':
            df.to_csv(filepath, index=True)
        elif fmt == 'json':
            df.to_json(filepath, orient='records', force_ascii=False)
        elif fmt in ('xlsx', 'excel'):
            df.to_excel(filepath, index=True)
        else:
            raise ValueError(f'Unsupported format: {fmt}')














class Outputer:
    """
    Output utility for saving calculation results (abundances, heating rates) to file.
    Supports CSV, Excel, and JSON formats.
    """
    def __init__(self):
        pass


    @staticmethod
    def save_abundance(
        abundance: Union[Dict[str, float], pd.DataFrame],
        filename: Union[str, Path],
        filetype: Optional[str] = None,
        index: bool = True
    ):
        """
        Save nuclide abundance (dict or DataFrame) to file.
        Args:
            abundance: Dict or DataFrame of nuclide abundances (can be time series DataFrame).
            filename: Output file path.
            filetype: 'csv', 'xlsx', or 'json'. If None, inferred from filename.
            index: Whether to write DataFrame index (default True).
        """
        if filetype is None:
            filetype = str(filename).split('.')[-1].lower()
        if isinstance(abundance, dict):
            df = pd.DataFrame(list(abundance.items()), columns=["nuclide", "abundance"])
        else:
            df = abundance
        if filetype == 'csv':
            df.to_csv(filename, index=index)
        elif filetype in ('xls', 'xlsx'):
            df.to_excel(filename, index=index)
        elif filetype == 'json':
            df.to_json(filename, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")

    @staticmethod
    def save_heating_rate(
        heating_rate: Union[List[float], np.ndarray, pd.Series, pd.DataFrame],
        times: Optional[Union[List[float], np.ndarray]] = None,
        filename: Union[str, Path] = "heating_rate.csv",
        filetype: Optional[str] = None
    ):
        """
        Save heating rate (with optional time axis) to file.
        Args:
            heating_rate: List/array/Series/DataFrame of heating rates.
            times: Optional time points (same length as heating_rate).
            filename: Output file path.
            filetype: 'csv', 'xlsx', or 'json'. If None, inferred from filename.
        """
        if filetype is None:
            filetype = str(filename).split('.')[-1].lower()
        if isinstance(heating_rate, (list, np.ndarray, pd.Series)):
            if times is not None:
                df = pd.DataFrame({"time": times, "heating_rate": heating_rate})
            else:
                df = pd.DataFrame({"heating_rate": heating_rate})
        else:
            df = heating_rate
        if filetype == 'csv':
            df.to_csv(filename, index=False)
        elif filetype in ('xls', 'xlsx'):
            df.to_excel(filename, index=False)
        elif filetype == 'json':
            df.to_json(filename, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")













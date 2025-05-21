# -*- coding: utf-8 -*-

import ast
from pathlib import Path
from typing import Dict, List, Tuple, Union, Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .nuclide import Nuclide, NuclideStrError


class DecayDatabase:
    """
    A class representing a database of nuclear decay data for a set of nuclides.
    Provides access to half-lives, decay modes, progeny, branching ratios, decay energies, and plotting utilities.
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


    def _validate_nuclide(self, nuclide: str) -> str:
        """
        Validate the nuclide string and return its standardized symbol.
        Raises NuclideStrError if invalid or not found in the database.
        """
        try:
            symbol = Nuclide(nuclide).nuclide_symbol
        except Exception as e:
            raise NuclideStrError(nuclide, f"Invalid format: {str(e)}")
        if symbol not in self.nuclide_index_map:
            raise NuclideStrError(nuclide, f"Not found in {self.data_source} decay database.")
        return symbol
    

    def half_life(self, nuclide: str, units: str = 'readable') -> Union[float, str]:
        """
        Get the half-life of a nuclide in the specified units.
        Args:
            nuclide (str): Nuclide symbol or string.
            units (str): 'readable', 's', or other (returns value in original units).
        Returns:
            float or str: Half-life in the requested units or as a readable string.
        """
        nuclide = self._validate_nuclide(nuclide)
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


    def progeny(self, nuclide) -> List[str]:
        """
        Get the list of progeny (daughter nuclides) for a given nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            list[str]: List of progeny nuclide symbols.
        """
        nuclide = self._validate_nuclide(nuclide)
        return self.progeny_data[self.nuclide_index_map[nuclide]]
    
    
    def branching_ratios(self, nuclide: str) -> List[float]:
        """
        Get the branching ratios for all decay modes of a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            List: Branching ratios for each decay mode.
        """
        nuclide = self._validate_nuclide(nuclide)
        return self.branching_ratios_data[self.nuclide_index_map[nuclide]]
    
    
    def decay_modes(self, nuclide: str) -> List[str]:
        """
        Get the decay modes for a given nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            List: Decay modes for the nuclide.
        """
        nuclide = self._validate_nuclide(nuclide)
        return self.decay_modes_data[self.nuclide_index_map[nuclide]]
    

    def decay_constants(self, nuclide: str) -> np.ndarray:
        """
        Get the decay constant (lambda) for a nuclide.
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Decay constant in 1/s.
        """
        nuclide = self._validate_nuclide(nuclide)
        return self.decay_constants_data[self.nuclide_index_map[nuclide]]
    

    def proton_numbers(self) -> np.ndarray:
        """
        Get the array of proton numbers (Z) for all nuclides in the database.
        Returns:
            np.ndarray: Array of proton numbers.
        """
        return np.array([Nuclide(nuclide).Z for nuclide in self.nuclides])


    def neutron_numbers(self) -> np.ndarray:
        """
        Get the array of neutron numbers (N) for all nuclides in the database.
        Returns:
            np.ndarray: Array of neutron numbers.
        """
        return np.array([Nuclide(nuclide).N for nuclide in self.nuclides])
    

    def decay_energy(self, nuclide: str, energy_type: str) -> float:
        """
        Get a specific type of decay energy for a nuclide.
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
            raise ValueError(f"Invalid energy type: {energy_type}.\
                             Use 'EM', 'LP', 'HP' or 'Neutrino'.")
        nuclide = self._validate_nuclide(nuclide)

        return self.decay_energies_data[self.nuclide_index_map[nuclide]][energy_map[energy_type]]


    def decay_energy_EM(self, nuclide: str) -> float:
        """
        Get the electromagnetic decay energy for a nuclide.
        """
        return self.decay_energy(nuclide, 'EM')


    def decay_energy_LP(self, nuclide: str) -> float:
        """
        Get the light particle decay energy for a nuclide.
        """
        return self.decay_energy(nuclide, 'LP')


    def decay_energy_HP(self, nuclide: str) -> float:
        """
        Get the heavy particle decay energy for a nuclide.
        """
        return self.decay_energy(nuclide, 'HP')
    

    def decay_energy_neutrino(self, nuclide: str) -> float:
        """
        Get the neutrino decay energy for a nuclide.
        """
        return self.decay_energy(nuclide, 'Neutrino')
    
    
    def decay_energy_released(self, nuclide: str) -> float:
        """
        Get the total decay energy released (sum of EM, LP, HP, and Neutrino energies).
        Args:
            nuclide (str): Nuclide symbol or string.
        Returns:
            float: Total decay energy released.
        """
        energy = np.array([
            self.decay_energy(nuclide, 'EM'),
            self.decay_energy(nuclide, 'LP'),
            self.decay_energy(nuclide, 'HP'),
            self.decay_energy(nuclide, 'Neutrino')
        ])
        return np.nansum(energy)
    


    def plot_nuclear_chart(
        self,
        figure: Optional[plt.figure] = None,
        cmap: mpl.colors.Colormap = cm.viridis,
        nuclei_linewidths: float = 0.8,
        colorbar: bool = True,
        figsize: Tuple[float, float] = (12, 8),
        dpi: int = 250,
        min_half_life: float = 1.0,
        max_half_life: float = 4.36e17,
        magic_numbers: List[int] = [2, 8, 20, 28, 50, 82, 126],
        **kwargs
    ) -> plt.figure:
        """
        Plot the nuclear chart (N vs Z) colored by half-life, with magic numbers highlighted.
        Args:
            figure (plt.figure, optional): Existing matplotlib figure to plot on.
            cmap (mpl.colors.Colormap): Colormap for half-life.
            nuclei_linewidths (float): Line width for nuclide boxes.
            colorbar (bool): Whether to show colorbar.
            figsize (Tuple[float, float]): Figure size.
            dpi (int): Dots per inch for the figure.
            min_half_life (float): Minimum half-life for color normalization.
            max_half_life (float): Maximum half-life for color normalization.
            magic_numbers (List): List of magic numbers to highlight.
            **kwargs: Additional keyword arguments for plotting.
        Returns:
            plt.figure: The matplotlib figure object.
        Raises:
            ValueError: If no valid nuclide or half-life data is available.
        """
        nuclide_map = {}
        for nuclide in self.nuclides:
            try:
                parsed = Nuclide(nuclide)
                Z = parsed.Z
                N = parsed.N
                state = parsed.state
                if state:
                    continue  # Skip metastable states
                half_life = self.half_life(nuclide, units="s")
                nuclide_map[(N, Z)] = half_life
            except Exception:
                continue

        if not nuclide_map:
            raise ValueError("No valid nuclide data for plotting.")

        neutron_numbers = np.array([N for (N, Z) in nuclide_map.keys()])
        proton_numbers = np.array([Z for (N, Z) in nuclide_map.keys()])
        half_lives = np.array([nuclide_map[(N, Z)] for (N, Z) in nuclide_map.keys()])
        half_life_values = half_lives[np.isfinite(half_lives) & (half_lives > 0) & (half_lives < max_half_life)]
        if not len(half_life_values):
            raise ValueError("No valid half-life data.")
        max_half_life = half_life_values.max()
        norm = mpl.colors.LogNorm(vmin=min_half_life, vmax=max_half_life)

        if figure is None:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        else:
            fig = figure
            ax = fig.gca()
        ax.set_aspect('equal')

        rectangles = []
        colors = []
        for (N, Z) in nuclide_map.keys():
            half_life = nuclide_map[(N, Z)]
            rectangles.append(Rectangle((N-0.5, Z-0.5), 1, 1, fill=True))
            # Color coding for special cases
            if half_life == np.inf:
                colors.append((0.0, 0.0, 0.0, 1.0))  # black for stable
            elif half_life >= 0 and half_life < min_half_life:
                colors.append((0.5, 0.5, 0.5, 1.0))  # gray for very short-lived
            elif half_life >= max_half_life:
                colors.append((1.0, 1.0, 0.0, 1.0))  # yellow for very long-lived
            else:
                colors.append(cmap(norm(half_life)))
        pc = PatchCollection(rectangles, facecolor=colors, edgecolor='black', linewidths=nuclei_linewidths, zorder=1)
        ax.add_collection(pc)

        # Highlight magic numbers
        for magic_number in magic_numbers:
            mask_N = neutron_numbers == magic_number
            if np.any(mask_N):
                Z_min, Z_max = proton_numbers[mask_N].min(), proton_numbers[mask_N].max()
                ax.add_patch(Rectangle((magic_number-0.5, Z_min-0.5), 1, Z_max - Z_min + 1,
                                       fill=False, edgecolor='black', linewidth=2, linestyle='--', zorder=5))
                
            mask_Z = proton_numbers == magic_number
            if np.any(mask_Z):
                N_min, N_max = neutron_numbers[mask_Z].min(), neutron_numbers[mask_Z].max()
                ax.add_patch(Rectangle((N_min - 0.5, magic_number - 0.5), N_max - N_min + 1, 1,
                                       fill=False, edgecolor='black', linewidth=2, linestyle='--', zorder=5))

        ax.set_xlim(neutron_numbers.min() - 1, neutron_numbers.max() + 5)
        ax.set_ylim(proton_numbers.min() - 1, proton_numbers.max() + 5)
        ax.set_xlabel('Neutron Number (N)', fontsize=16)
        ax.set_ylabel('Proton Number (Z)', fontsize=16)

        if colorbar:
            mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes(position='right', size="5%", pad=0.01, aspect=1)
            cbar = fig.colorbar(mappable, ax=ax, cax=cax)
            cbar.set_label('Half-life (s)', fontsize=12)

        return fig
    


    def __eq__(self, other: object) -> bool:
        """
        Equality comparison for DecayDatabase objects.
        Returns True if all relevant data arrays and source match.
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
        Inequality comparison for DecayDatabase objects.
        """
        if not isinstance(other, DecayDatabase):
            return NotImplemented
        return not self.__eq__(other)


    # def __repr__(self) -> str:
    #     return (
    #         f"Decay dataset: {self.dataset_name}, contains SymPy data: "
    #         f"{self._sympy_data is not None}"
    #     )


class DecayDatabaseManager:
    """
    Manager class for loading, saving, and processing nuclear decay data files.
    Handles sorting, conversion, and serialization of decay data for use in DecayDatabase.
    """
    
    RADIONUCLIDE_LABEL            = 'Radionuclide'
    MASS_NUMBER_LABEL             = 'A'
    ATOMIC_NUMBER_LABEL           = 'Z'
    IS_STABLE_LABEL               = 'Is_Stable'
    HALF_LIFE_SECOND_LABEL        = 'Half_Life_Second'
    HALF_LIFE_READABLE_LABEL      = 'Half_Life_Readable'
    HALF_LIFE_VALUE_LABEL         = 'Half_Life_Value'
    HALF_LIFE_UNIT_LABEL          = 'Half_Life_Unit'
    DECAY_MODES_LABEL             = 'Decay_Modes'
    NUM_DECAY_MODES_LABEL         = 'Num_Decay_Modes'
    BRANCHINF_RATIOS_LABEL        = 'Branching_Ratios'
    PROGENY_LABEL                 = 'Progeny'
    DECAY_ENERGY_EM_LABEL         = 'Decay_Energy_EM'
    DECAY_ENERGY_LP_LABEL         = 'Decay_Energy_LP'
    DECAY_ENERGY_HP_LABEL         = 'Decay_Energy_HP'
    DECAY_ENERGY_NEUTRINO_LABEL   = 'Decay_Energy_Neutrino'
    DECAY_ENERGY_GAMMA_LABEL      = 'Decay_Energy_Gamma'
    DECAY_ENERGY_BETA_MINUS_LABEL = 'Decay_Energy_Beta_Minus'
    DECAY_ENERGY_BETA_PLUS_LABEL  = 'Decay_Energy_Beta_Plus'
    DECAY_ENERGY_ALPHA_LABEL      = 'Decay_Energy_Alpha'
    DECAY_ENERGY_NEUTRON_LABEL    = 'Decay_Energy_Neutron'
    DECAY_ENERGY_PROTON_LABEL     = 'Decay_Energy_Proton'
    DECAY_EFFECTIVE_Q_LABEL       = 'Decay_Effective_Q'
    
    
    
    def __init__(self, data_source: str = 'ENDF-B-VIII.1_decay'):
        """
        Initialize the DecayDatabaseManager with a data source name.
        Args:
            data_source (str): Name of the decay data source.
        """
        self.data_source = data_source
        self.data_path = Path(__file__).resolve().parent.joinpath('data')


    def save_decay_database(self) -> None:
        """
        Load, sort, and save the decay database in both CSV and compressed NumPy formats.
        This prepares the data for fast loading and use in DecayDatabase.
        """
        df = pd.read_json(self.data_path.joinpath(f'{self.data_source}.json'), orient='records')
        df = self._sort_algorithm(df)
        df.to_csv(self.data_path.joinpath(f'{self.data_source}_sorted.csv'), index=False)
        df.to_json(self.data_path.joinpath(f'{self.data_source}_sorted.json'), orient='records', indent=4)
        
        (
            nuclides,
            half_life_data,
            decay_constants_data,
            decay_modes_data,
            progeny_data,
            branching_ratios_data,
            decay_energies_data
        ) = self._process_radionuclide_data(df)

        np.savez_compressed(
            self.data_path.joinpath(f'{self.data_source}.npz'),
            nuclides=np.asarray(nuclides),
            half_life_data=half_life_data,
            decay_constants_data=decay_constants_data,
            decay_modes_data=np.asarray(decay_modes_data, dtype=object),
            progeny_data=np.asarray(progeny_data, dtype=object),
            branching_ratios_data=np.asarray(branching_ratios_data, dtype=object),
            decay_energies_data=np.asarray(decay_energies_data, dtype=object) 
        )


    def _process_radionuclide_data(
            self,
            df: pd.DataFrame
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Process the DataFrame of radionuclide data into arrays for DecayDatabase.
        Args:
            df (pd.DataFrame): DataFrame containing radionuclide data.
        Returns:
            Tuple: Arrays for nuclides, half-life, decay constants, decay modes, progeny, branching ratios, and decay energies.
        """
        nuclides = [Nuclide(nuclide).nuclide_symbol for nuclide in df[self.RADIONUCLIDE_LABEL].values]

        half_life_data = np.array([
            (
                np.float64(half_life_second),
                np.float64(half_life_value),
                half_life_unit,
                half_life_readable
            )
            for half_life_second, half_life_value, half_life_unit, half_life_readable in zip(
                df[self.HALF_LIFE_SECOND_LABEL],
                df[self.HALF_LIFE_VALUE_LABEL],
                df[self.HALF_LIFE_UNIT_LABEL],
                df[self.HALF_LIFE_READABLE_LABEL]
            )
        ], dtype=object)

        # Calculate decay constants (lambda = ln(2) / half-life)
        decay_constants_data: np.ndarray = np.array(
            [
                np.log(2.0) / half_life_values[0] for half_life_values in half_life_data
            ],
            dtype = np.float64
        )

        decay_energies_data = np.array([
            (
                decay_energy_EM, decay_energy_LP, decay_energy_HP, decay_energy_neutrino,
                decay_energy_gamma, decay_energy_beta_minus, decay_energy_beta_plus,
                decay_energy_alpha, decay_energy_neutron, decay_energy_proton, decay_effective_q
            )
            for decay_energy_EM, decay_energy_LP, decay_energy_HP, decay_energy_neutrino,
                decay_energy_gamma, decay_energy_beta_minus, decay_energy_beta_plus,
                decay_energy_alpha, decay_energy_neutron, decay_energy_proton, decay_effective_q
            in zip(
                df[self.DECAY_ENERGY_EM_LABEL], df[self.DECAY_ENERGY_LP_LABEL], df[self.DECAY_ENERGY_HP_LABEL],
                df[self.DECAY_ENERGY_NEUTRINO_LABEL], df[self.DECAY_ENERGY_GAMMA_LABEL],
                df[self.DECAY_ENERGY_BETA_MINUS_LABEL], df[self.DECAY_ENERGY_BETA_PLUS_LABEL],
                df[self.DECAY_ENERGY_ALPHA_LABEL], df[self.DECAY_ENERGY_NEUTRON_LABEL],
                df[self.DECAY_ENERGY_PROTON_LABEL], df[self.DECAY_EFFECTIVE_Q_LABEL]
            )
        ], dtype=np.float64)

        num_nuclides = len(nuclides)
        decay_modes_data = [[] for _ in range(num_nuclides)]
        progeny_data = [[] for _ in range(num_nuclides)]
        branching_ratios_data = [[] for _ in range(num_nuclides)]

        num_decay_modes_list = df[self.NUM_DECAY_MODES_LABEL].values
        decay_modes_list = df[self.DECAY_MODES_LABEL].values
        progeny_list = df[self.PROGENY_LABEL].values
        branching_ratios_list = df[self.BRANCHINF_RATIOS_LABEL].values

        for i in range(num_nuclides):
            progeny, branching_ratios, decay_modes = [], [], []
            num_decay_modes = num_decay_modes_list[i]

            if num_decay_modes > 0:
                progeny = progeny_list[i][:num_decay_modes]
                branching_ratios = branching_ratios_list[i][:num_decay_modes]
                decay_modes = decay_modes_list[i][:num_decay_modes]

                # Sort decay modes by branching ratio (descending), then progeny, then mode
                sorted_triples = sorted(
                    zip(branching_ratios, progeny, decay_modes),
                    key=lambda x: (-x[0], x[1], x[2])
                )

                branching_ratios_data[i], progeny_data[i], decay_modes_data[i] = map(list, zip(*sorted_triples))
                
        database = (
            np.asarray(nuclides),
            half_life_data,
            decay_constants_data,
            decay_modes_data,
            progeny_data,
            branching_ratios_data,
            decay_energies_data,
        )

        return database


    def _sort_algorithm(self, df: pd.DataFrame, decay_mode: str = 'β-') -> pd.DataFrame:
        """
        Sort the DataFrame of nuclides for consistent ordering, prioritizing stability and decay chains.
        Args:
            df (pd.DataFrame): DataFrame to sort.
            decay_mode (str): Decay mode to use for sorting chains (default: 'β-').
        Returns:
            pd.DataFrame: Sorted DataFrame.
        """
        df['Is_Stable_Sort'] = df[self.IS_STABLE_LABEL].apply(lambda x: 1 if x else 0)

        df.sort_values(
            by=["Is_Stable_Sort", self.MASS_NUMBER_LABEL, self.ATOMIC_NUMBER_LABEL, "Metastable_State"],
            ascending=[True, False, False, False],
            inplace=True
        )

        # Set infinite half-life for stable nuclides
        df.loc[df[self.IS_STABLE_LABEL] == True, self.HALF_LIFE_SECOND_LABEL] = np.inf
        df.loc[df[self.IS_STABLE_LABEL] == True, self.HALF_LIFE_VALUE_LABEL] = np.inf

        sorted_df = self._process_sort(df, decay_mode=decay_mode)
        sorted_df.reset_index(drop=True, inplace=True)
        sorted_df = sorted_df[df.columns.drop('Is_Stable_Sort')]

        # Convert string representations of lists to actual lists
        for column in [self.BRANCHINF_RATIOS_LABEL, self.PROGENY_LABEL, self.DECAY_MODES_LABEL]:
            sorted_df = self._convert_str_to_list(sorted_df, column)

        return sorted_df


    def _process_sort(self, df: pd.DataFrame, decay_mode: str = 'β-') -> pd.DataFrame:
        """
        Sort nuclides in the DataFrame so that decay chains (by the specified mode) are in order.
        Args:
            df (pd.DataFrame): DataFrame to process.
            decay_mode (str): Decay mode to use for sorting.
        Returns:
            pd.DataFrame: Sorted DataFrame.
        """
        df_dict = df.to_dict('records')
        swapped = True

        while swapped:
            swapped = False
            nuclide_to_index = {row[self.RADIONUCLIDE_LABEL]: idx for idx, row in enumerate(df_dict)}

            for i in range(len(df_dict)):
                current_row = df_dict[i]
                decay_modes_irow = current_row[self.DECAY_MODES_LABEL]

                if decay_mode not in decay_modes_irow:
                    continue

                mode_indices = [idx for idx, mode in enumerate(decay_modes_irow) if mode == decay_mode]

                for mode_index in mode_indices:
                    daughter_nuclide = current_row[self.PROGENY_LABEL][mode_index]
                    if daughter_nuclide not in nuclide_to_index:
                        continue

                    j = nuclide_to_index[daughter_nuclide]
                    daughter_row = df_dict[j]

                    # Ensure parent precedes daughter in the chain
                    if (
                        current_row[self.MASS_NUMBER_LABEL] == daughter_row[self.MASS_NUMBER_LABEL] and
                        daughter_row[self.ATOMIC_NUMBER_LABEL] == current_row[self.ATOMIC_NUMBER_LABEL] + 1
                    ):
                        if i > j:
                            df_dict[i], df_dict[j] = df_dict[j], df_dict[i]
                            swapped = True
                            break
                if swapped:
                    break

        return pd.DataFrame(df_dict)


    def _convert_str_to_list(self, df: pd.DataFrame, column_name: str) -> pd.DataFrame:
        """
        Convert string representations of lists in a DataFrame column to actual Python lists.
        Args:
            df (pd.DataFrame): DataFrame to process.
            column_name (str): Name of the column to convert.
        Returns:
            pd.DataFrame: DataFrame with converted column.
        """
        for index, value in df[column_name].items():
            if isinstance(value, str):
                try:
                    if value.startswith('[') and value.endswith(']'):
                        df.at[index, column_name] = ast.literal_eval(value)
                    else:
                        print(f"Warning: '{value}' is not a valid list string at index {index}")
                except (ValueError, SyntaxError) as e:
                    print(f"Warning: Failed to convert value at index {index}: '{value}' ({str(e)})")
        return df


def load_decay_database(data_source: str = 'ENDF-B-VIII.1_decay') -> DecayDatabase:
    """
    Load a DecayDatabase object from a compressed NumPy file for the specified data source.
    Args:
        data_source (str): Name of the decay data source.
    Returns:
        DecayDatabase: Loaded decay database object.
    """
    # DecayDatabaseManager(data_source).save_decay_database()

    database = np.load(
        DecayDatabaseManager(data_source).data_path.joinpath(f'{data_source}.npz'),
        allow_pickle=True,
    )
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


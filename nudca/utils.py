from typing import Dict, List, Tuple, Union, Optional
from pathlib import Path

import ast
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle, Patch

from .nuclide import Nuclide

# Label constants for nuclear data columns
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
DECAY_ENERGY_SF_LABEL         = 'Decay_Energy_SF'
DECAY_ENERGY_GAMMA_LABEL      = 'Decay_Energy_Gamma'
DECAY_ENERGY_BETA_MINUS_LABEL = 'Decay_Energy_Beta_Minus'
DECAY_ENERGY_BETA_PLUS_LABEL  = 'Decay_Energy_Beta_Plus'
DECAY_ENERGY_ALPHA_LABEL      = 'Decay_Energy_Alpha'
DECAY_ENERGY_NEUTRON_LABEL    = 'Decay_Energy_Neutron'
DECAY_ENERGY_PROTON_LABEL     = 'Decay_Energy_Proton'
DECAY_EFFECTIVE_Q_LABEL       = 'Decay_Effective_Q'


class DecayDatabaseManager:
    """
    Manages loading, saving, and processing of nuclear decay data files.
    Handles sorting, conversion, and serialization of decay data for use in DecayDatabase.
    """
    
    def __init__(self, data_source: str = 'ENDF-B-VIII.1_decay'):
        """
        Initialize the DecayDatabaseManager with a data source name.
        Args:
            data_source (str): Name of the decay data source.
        """
        self.data_source = data_source
        self.data_path = Path(__file__).resolve().parent.joinpath('data')


    def transform_radionuclide_data(
            self,
            df: pd.DataFrame
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Convert a DataFrame of radionuclide data into arrays for DecayDatabase.
        Args:
            df (pd.DataFrame): DataFrame containing radionuclide data.
        Returns:
            Tuple: Arrays for nuclides, half-life, decay constants, decay modes, progeny, branching ratios, and decay energies.
        """
        nuclides = [Nuclide(nuclide).nuclide_symbol for nuclide in df[RADIONUCLIDE_LABEL].values]
        half_life_data = np.array([
            (
                np.float64(half_life_second),
                np.float64(half_life_value),
                half_life_unit,
                half_life_readable
            )
            for half_life_second, half_life_value, half_life_unit, half_life_readable in zip(
                df[HALF_LIFE_SECOND_LABEL],
                df[HALF_LIFE_VALUE_LABEL],
                df[HALF_LIFE_UNIT_LABEL],
                df[HALF_LIFE_READABLE_LABEL]
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
                decay_energy_EM, decay_energy_LP, decay_energy_HP, decay_energy_neutrino, decay_energy_SF,
                decay_energy_gamma, decay_energy_beta_minus, decay_energy_beta_plus,
                decay_energy_alpha, decay_energy_neutron, decay_energy_proton, decay_effective_q
            )
            for decay_energy_EM, decay_energy_LP, decay_energy_HP, decay_energy_neutrino, decay_energy_SF,
                decay_energy_gamma, decay_energy_beta_minus, decay_energy_beta_plus,
                decay_energy_alpha, decay_energy_neutron, decay_energy_proton, decay_effective_q
            in zip(
                df[DECAY_ENERGY_EM_LABEL], df[DECAY_ENERGY_LP_LABEL], df[DECAY_ENERGY_HP_LABEL],
                df[DECAY_ENERGY_NEUTRINO_LABEL], df[DECAY_ENERGY_SF_LABEL], df[DECAY_ENERGY_GAMMA_LABEL],
                df[DECAY_ENERGY_BETA_MINUS_LABEL], df[DECAY_ENERGY_BETA_PLUS_LABEL],
                df[DECAY_ENERGY_ALPHA_LABEL], df[DECAY_ENERGY_NEUTRON_LABEL],
                df[DECAY_ENERGY_PROTON_LABEL], df[DECAY_EFFECTIVE_Q_LABEL]
            )
        ], dtype=np.float64)
        num_nuclides = len(nuclides)
        decay_modes_data = [[] for _ in range(num_nuclides)]
        progeny_data = [[] for _ in range(num_nuclides)]
        branching_ratios_data = [[] for _ in range(num_nuclides)]
        num_decay_modes_list = df[NUM_DECAY_MODES_LABEL].values
        decay_modes_list = df[DECAY_MODES_LABEL].values
        progeny_list = df[PROGENY_LABEL].values
        branching_ratios_list = df[BRANCHINF_RATIOS_LABEL].values
        for i in range(num_nuclides):
            progeny, branching_ratios, decay_modes = [], [], []
            num_decay_modes = num_decay_modes_list[i]
            if num_decay_modes > 0:
                progeny = progeny_list[i][:num_decay_modes]
                branching_ratios = branching_ratios_list[i][:num_decay_modes]
                decay_modes = decay_modes_list[i][:num_decay_modes]
                # Sort decay modes by branching ratio (descending), then progeny, then mode
                # sorted_triples = sorted(
                #     zip(branching_ratios, progeny, decay_modes),
                #     key=lambda x: (-x[0], x[1], x[2])
                # )
                # branching_ratios_data[i], progeny_data[i], decay_modes_data[i] = map(list, zip(*sorted_triples))
                
                branching_ratios_data[i], progeny_data[i], decay_modes_data[i] = branching_ratios, progeny, decay_modes
        database = (
            nuclides,
            half_life_data,
            decay_constants_data,
            decay_modes_data,
            progeny_data,
            branching_ratios_data,
            decay_energies_data,
        )
        return database

    def sort_nuclides_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Sort the DataFrame of nuclides for consistent ordering, prioritizing stability and decay chains.
        Args:
            df (pd.DataFrame): DataFrame to sort.
        Returns:
            pd.DataFrame: Sorted DataFrame.
        """
        df.sort_values(by=[MASS_NUMBER_LABEL], ascending=[False], inplace=True)
        # Set infinite half-life for stable nuclides
        df.loc[df[IS_STABLE_LABEL] == True, HALF_LIFE_SECOND_LABEL] = np.inf
        df.loc[df[IS_STABLE_LABEL] == True, HALF_LIFE_VALUE_LABEL] = np.inf

        sorted_df = NuclideSorter(df).sort_nuclides_decay_chains()

        # Convert string representations of lists to actual lists
        for column in [BRANCHINF_RATIOS_LABEL, PROGENY_LABEL, DECAY_MODES_LABEL]:
            sorted_df = self.parse_string_to_list(sorted_df, column)
        unstable_df = sorted_df[sorted_df[IS_STABLE_LABEL] == False]
        stable_df = sorted_df[sorted_df[IS_STABLE_LABEL] == True]
        final_df = pd.concat([unstable_df, stable_df], ignore_index=True)
        final_df.reset_index(drop=True, inplace=True)
        return final_df


    @staticmethod
    def parse_string_to_list(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
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




class NuclideSorter:
    # Class-level constants for column labels to avoid hardcoding strings throughout the code


    def __init__(self, df: pd.DataFrame):
        """
        Initialize the NuclideSorter with a DataFrame containing nuclide data.
        The DataFrame is expected to have columns for mass number, nuclide symbol, and progeny list.
        Args:
            df (pd.DataFrame): DataFrame with columns 'A', 'Radionuclide', and 'Progeny'.
        """
        self.df = df


    @staticmethod
    def parent_to_progeny_map(df_subgroup: pd.DataFrame) -> Dict[str, List[str]]:
        """
        Construct a mapping from each parent nuclide to its list of progeny nuclides.
        Args:
            df (pd.DataFrame): DataFrame with columns for nuclide and progeny.
        Returns:
            Dict[str, List[str]]: Dictionary mapping parent nuclide symbols to lists of their progeny symbols.
        """
        return {row[RADIONUCLIDE_LABEL]: row[PROGENY_LABEL] for _, row in df_subgroup.iterrows()}


    def topological_sort(self, df_subgroup: Optional[pd.DataFrame] = None) -> List[str]:
        """
        Perform a topological sort of nuclides with the same mass number (A) based on parent-progeny relationships.
        This ensures that for any parent-progeny relationship, the parent appears before the progeny in the sorted list.
        Args:
            df (pd.DataFrame, optional): DataFrame containing nuclides with the same mass number. If None, use self.df.
        Returns:
            List[str]: Topologically sorted list of nuclide symbols.
        Raises:
            ValueError: If a cyclic dependency is detected in the parent-progeny relationships.
        """
        sorted_nuclides = []
        # Work on a copy to avoid modifying the original DataFrame
        working_df = (df_subgroup if df_subgroup is not None else self.df).copy()
        parent_to_progeny = self.parent_to_progeny_map(working_df)

        while not working_df.empty:
            acyclic = False
            # Iterate through all nuclides in the current working set
            for nuclide in working_df[RADIONUCLIDE_LABEL].tolist():
                nuclide_symbol = Nuclide(nuclide).nuclide_symbol
                # Find all parents of the current nuclide
                parents = [parent for parent, progeny in parent_to_progeny.items() if nuclide_symbol in progeny]
                # If the nuclide has no parents, it can be safely added to the sorted list
                if not parents:
                    sorted_nuclides.append(nuclide)
                    # Remove the nuclide from the working set and mapping
                    working_df = working_df[working_df[RADIONUCLIDE_LABEL] != nuclide]
                    parent_to_progeny.pop(nuclide, None)
                    acyclic = True
                    break
            # If no acyclic node was found, there must be a cycle in the data
            if not acyclic:
                raise ValueError(f"Cyclic dependency or data error detected. "
                                 f"Remaining nuclides: {working_df[RADIONUCLIDE_LABEL].tolist()}")
        return sorted_nuclides



    def sort_nuclides_decay_chains(self) -> pd.DataFrame:
        """
        Sort all nuclides by descending mass number (A), and within each mass number group,
        perform a topological sort based on parent-progeny relationships.
        Args:
            output_path (Optional[str]): If provided, save the sorted DataFrame to this CSV file.
        Returns:
            pd.DataFrame: The fully sorted DataFrame, with nuclides grouped by mass number and sorted topologically within each group.
        """
        sorted_df = []
        for mass_number in sorted(self.df[MASS_NUMBER_LABEL].unique(), reverse=True):
            df_subgroup = self.df[self.df[MASS_NUMBER_LABEL] == mass_number]
            sorted_nuclides = self.topological_sort(df_subgroup)
            df_subgroup_sorted = df_subgroup.set_index(RADIONUCLIDE_LABEL).reindex(sorted_nuclides).reset_index()
            sorted_df.append(df_subgroup_sorted)

        final_sorted_df = pd.concat(sorted_df, ignore_index=True)
        return final_sorted_df




class HalfLifeColorMap:
    """Configuration for half-life color mapping in nuclear chart."""
    
    def __init__(
        self,
        half_life_data: Optional[np.ndarray] = None
    ) -> None:
        self.half_life_data = half_life_data
        self.half_life_seconds = np.array([float(x[0]) for x in self.half_life_data])
        self.ranges = self._create_ranges()
        
    @property
    def min_half_life(self) -> float:
        # return np.where(self.half_life_seconds == 0)
        return np.min(
            self.half_life_seconds[np.isfinite(self.half_life_seconds) & (self.half_life_seconds > 0)]
        )
    
    @property
    def max_half_life(self) -> float:
        return np.max(
            self.half_life_seconds[np.isfinite(self.half_life_seconds) & (self.half_life_seconds > 0)]
        )
    
    
    def _create_figure(
        self,
        figure: Optional[plt.figure],
        figsize: Tuple[float, float],
        dpi: int
    ) -> Tuple[plt.figure, plt.Axes]:
        """Create or get existing figure and axes."""
        if figure is None:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        else:
            fig = figure
            ax = fig.gca()
        ax.set_aspect('equal')
        return fig, ax

    def _plot_stable_nuclei(
        self,
        ax: plt.Axes,
        stable_nuclides: List[Tuple[int, int]],
        nuclei_linewidths: float
    ) -> None:
        """Plot stable nuclei in black."""
        for N, Z in stable_nuclides:
            ax.add_patch(
                Rectangle(
                    (N-0.5, Z-0.5), 1, 1,
                    facecolor='black',
                    edgecolor='black',
                    linewidth=nuclei_linewidths
                )
            )

    def _plot_unstable_nuclei(
        self,
        ax: plt.Axes,
        range_groups: Dict[str, List[Tuple[int, int, float]]],
        color_map: 'HalfLifeColorMap',
        nuclei_linewidths: float
    ) -> None:
        """Plot unstable nuclei colored by their half-life ranges."""
        for range_key, nuclides in range_groups.items():
            if not nuclides:
                continue
                
            cmap, norm = color_map.get_cmap_and_norm(range_key)
            
            for N, Z, value in nuclides:
                color = cmap(norm(value))
                ax.add_patch(Rectangle(
                    (N-0.5, Z-0.5), 1, 1,
                    facecolor=color,
                    edgecolor='black',
                    linewidth=nuclei_linewidths,
                    # alpha=0.6
                ))

    def _get_axis_limits(
        self,
        range_groups: Dict[str, List[Tuple[int, int, float]]],
        stable_nuclides: List[Tuple[int, int]]
    ) -> Tuple[List[int], List[int]]:
        """Get all N and Z values for axis limits."""
        all_N_values = []
        all_Z_values = []
        
        for nuclides in range_groups.values():
            for N, Z, _ in nuclides:
                all_N_values.append(N)
                all_Z_values.append(Z)
                
        for N, Z in stable_nuclides:
            all_N_values.append(N)
            all_Z_values.append(Z)
            
        return all_N_values, all_Z_values


    def _mark_magic_numbers(
        self,
        ax: plt.Axes,
        magic_numbers: List[int],
        range_groups: Dict[str, List[Tuple[int, int, float]]],
        stable_nuclides: List[Tuple[int, int]]
    ) -> None:
        """Mark magic numbers on the nuclear chart."""
        for magic_number in magic_numbers:
            magic_N_nuclides = []
            magic_Z_nuclides = []
            
            # Check stable nuclides
            for N, Z in stable_nuclides:
                if N == magic_number:
                    magic_N_nuclides.append((N, Z))
                if Z == magic_number:
                    magic_Z_nuclides.append((N, Z))
            
            # Check unstable nuclides
            for nuclides in range_groups.values():
                for N, Z, _ in nuclides:
                    if N == magic_number:
                        magic_N_nuclides.append((N, Z))
                    if Z == magic_number:
                        magic_Z_nuclides.append((N, Z))
            
            # Highlight neutron magic number
            if magic_N_nuclides:
                Z_values = [Z for _, Z in magic_N_nuclides]
                Z_min, Z_max = min(Z_values), max(Z_values)
                ax.add_patch(
                    Rectangle(
                        (magic_number-0.5, Z_min-0.5), 1, Z_max - Z_min + 1,
                        fill=False,
                        edgecolor='black',
                        linewidth=0.8,
                        linestyle='--',
                        zorder=5
                    )
                )
            
            # Mark proton magic number
            if magic_Z_nuclides:
                N_values = [N for N, _ in magic_Z_nuclides]
                N_min, N_max = min(N_values), max(N_values)
                ax.add_patch(Rectangle(
                    (N_min-0.5, magic_number-0.5), N_max - N_min + 1, 1,
                    fill=False, edgecolor='black', linewidth=0.8, linestyle='--', zorder=5
                ))

    def _set_axis_properties(
        self,
        ax: plt.Axes,
        all_N_values: List[int],
        all_Z_values: List[int],
        fontsize: int = 14
    ) -> None:
        """Set axis limits and labels."""
        ax.set_xlim(min(all_N_values) - 1, max(all_N_values) + 5)
        ax.set_ylim(min(all_Z_values) - 1, max(all_Z_values) + 5)
        ax.set_xlabel('Neutron Number (N)', fontsize=fontsize)
        ax.set_ylabel('Proton Number (Z)', fontsize=fontsize)
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        
    
    def _create_ranges(self) -> Dict[str, Dict]:
        """Create color mapping ranges for different half-life intervals."""
        # half_life_intervals = [
        #     (self.min_half_life, 1.0e-3, 'X1', 'Greys', 0.30, 0.40),       # < 1 ms
        #     (1.0e-3, 1.0e-2, 'X2', 'Greens', 0.30, 0.32),                  # 1 ms - 10 ms
        #     (1.0e-2, 1.0e-1, 'X4', 'Greens', 0.50, 0.52),                  # 10 ms - 100 ms
        #     (1.0e-1, 1.0e0, 'X5', 'Greens', 0.70, 0.72),                   # 100 ms - 1 s
        #     (1.0e0, 1.0e1, 'X6', 'Greens', 0.90, 0.92),                    # 1 s - 10 s
        #     (1.0e1, 60.0, 'X7', 'Blues', 0.30, 0.32),                      # 10 s - 1 min
        #     (60.0, 3600.0, 'X8', 'Blues', 0.50, 0.52),                     # 1 min - 1 h
        #     (3600.0, 86400.0, 'X9', 'Blues', 0.70, 0.72),                  # 1 h - 1 day
        #     (86400.0, 3.1536e7, 'X10', 'Blues', 0.90, 0.92),               # 1 day - 1 year
        #     (3.1536e7, 3.1536e13, 'X11', 'Purples', 0.70, 0.72),           # 1 year - 1 My
        #     (3.1536e13, 3.1536e16, 'X12', 'Purples', 0.80, 0.82),          # 1 My - 10 Gy
        #     (3.1536e16, self.max_half_life, 'X13', 'Purples', 0.98, 1.0)   # > 10 Gy
        # ]
        half_life_intervals = [
            (self.min_half_life, 1.0e-3, 'X1', 'Greys', 0.4, 0.42),          # < 1 ms
            (1.0e-3, 1.0e-2, 'X2', 'spring_r', 0.4, 0.42),                   # 1 ms - 10 ms
            (1.0e-2, 1.0e-1, 'X4', 'summer_r', 0.30, 0.32),                  # 10 ms - 100 ms
            (1.0e-1, 1.0e0, 'X5', 'summer_r', 0.60, 0.62),                   # 100 ms - 1 s
            (1.0e0, 1.0e1, 'X6', 'summer_r', 0.90, 0.92),                    # 1 s - 10 s
            (1.0e1, 60.0, 'X7', 'winter_r', 0.30, 0.32),                     # 10 s - 1 min
            (60.0, 3600.0, 'X8', 'Blues', 0.50, 0.52),                       # 1 min - 1 h
            (3600.0, 86400.0, 'X9', 'Blues', 0.70, 0.72),                    # 1 h - 1 day
            (86400.0, 3.1536e7, 'X10', 'Purples', 0.70, 0.72),               # 1 day - 1 year
            (3.1536e7, 3.1536e13, 'X11', 'Purples', 0.98, 1.0),              # 1 year - 1 My
            (3.1536e13, 3.1536e16, 'X12', 'plasma_r', 0.98, 1.0),            # 1 My - 10 Gy
            (3.1536e16, self.max_half_life, 'X13', 'inferno_r', 0.86, 0.88)  # > 10 Gy
        ]


        # Calculate number of nuclides in each range
        half_life_seconds = np.array([float(x[0]) for x in self.half_life_data])
        range_counts = {}
        for min_time, max_time, key, _, _, _ in half_life_intervals:
            count = np.sum((half_life_seconds >= min_time) & (half_life_seconds < max_time))
            range_counts[key] = max(64, min(256, count * 2))  # Ensure N is between 64 and 256

        ranges = {}
        for min_time, max_time, key, cmap_name, color_start, color_end in half_life_intervals:
            ranges[key] = {
                'min': min_time,
                'max': max_time,
                'cmap': mpl.colors.LinearSegmentedColormap.from_list(
                    f'viridis_{key.lower()}',
                    [getattr(cm, cmap_name)(color_start), getattr(cm, cmap_name)(color_end)],
                    N=range_counts[key]
                )
            }
        return ranges

    def get_range_for_half_life(self, half_life: float) -> Optional[str]:
        """Get the appropriate range key for a given half-life value."""
        for range_key, range_info in self.ranges.items():
            if range_info['min'] <= half_life < range_info['max']:
                return range_key
        return None

    def get_cmap_and_norm(self, range_key: str) -> Tuple[mpl.colors.Colormap, mpl.colors.Normalize]:
        """Get colormap and normalization for a given range key."""
        range_info = self.ranges[range_key]
        return range_info['cmap'], mpl.colors.LogNorm(
            vmin=range_info['min'],
            vmax=range_info['max']
        )
        
        
        
    def _plot_legend(
        self,
        ax: plt.Axes,
        range_groups: Dict[str, List[Tuple[int, int, float]]],
        stable_nuclides: List[Tuple[int, int]],
        fontsize: int = 11,
        title: str = "ENDF-B-VIII.1 Nuclides",
        loc: str = 'best',
        show_stable: bool = True,
    ) -> None:
        """Plot legend for the nuclear chart."""
        legend_handles = []
        total_count = 0
        intervals = [
            (self.min_half_life, 1.0e-3, 'X1', '< 1 ms'),
            (1.0e-3, 1.0e-2, 'X2', '1 ms - 10 ms'),
            (1.0e-2, 1.0e-1, 'X4', '10 ms - 100 ms'),
            (1.0e-1, 1.0e0, 'X5', '100 ms - 1 s'),
            (1.0e0, 1.0e1, 'X6', '1 s - 10 s'),
            (1.0e1, 60.0, 'X7', '10 s - 1 min'),
            (60.0, 3600.0, 'X8', '1 min - 1 h'),
            (3600.0, 86400.0, 'X9', '1 h - 1 day'),
            (86400.0, 3.1536e7, 'X10', '1 day - 1 year'),
            (3.1536e7, 3.1536e13, 'X11', '1 year - 1 My'),
            (3.1536e13, 3.1536e16, 'X12', '1 My - 10 Gy'),
            (3.1536e16, self.max_half_life, 'X13', '> 10 Gy'),
        ]
        for min_time, max_time, key, label in intervals:
            count = len(range_groups.get(key, []))
            if count == 0:
                continue
            total_count += count
            cmap, norm = self.get_cmap_and_norm(key)
            mid_value = np.exp((np.log(min_time) + np.log(max_time)) / 2)
            color = cmap(norm(mid_value))
            legend_handles.append(
                Patch(color=color, label=f"{label}: {count}")
            )
        # Add stable nuclides (black)
        if show_stable and stable_nuclides:
            legend_handles.append(
                Patch(color='black', label=f'Stable: {len(stable_nuclides)}')
            )
            total_count += len(stable_nuclides)
        ax.legend(
            handles=legend_handles,
            title=f"{title} (Total: {total_count})",
            fontsize=fontsize,
            title_fontsize=fontsize+1,
            loc=loc,
            ncol=2,
            frameon=False
        )


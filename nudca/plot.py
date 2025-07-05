
from typing import Any, Dict, List, Tuple, Union, Optional
from collections import deque

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

from .nuclide import Nuclide
from .utils import HalfLifeColorMap

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
        digraph, max_generation, max_xpos = self._build_generic_decay_digraph(
            nuclide,
            nx.DiGraph(),
            traversal_direction="forward",
            decay_mode_filter=None,
            edge_direction="forward"
        )
        
        return self._plot(
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
        digraph, max_generation, max_xpos = self._build_generic_decay_digraph(
            nuclide,
            nx.DiGraph(),
            traversal_direction="reverse",
            decay_mode_filter=None,
            edge_direction="forward"
        )
        return self._plot(
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
        digraph, max_generation, max_xpos = self._build_generic_decay_digraph(
            nuclide,
            nx.DiGraph(),
            traversal_direction="reverse",
            decay_mode_filter=['β+&EC', 'SF'],
            edge_direction="forward"
        )
        return self._plot(
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
    def _plot(
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


    def _build_generic_decay_digraph(
        self,
        nuclide: str,
        digraph: nx.DiGraph,
        traversal_direction: str = "forward",
        decay_mode_filter: Optional[List[str]] = None,
        edge_direction: str = "forward",
    ) -> Tuple[nx.DiGraph, int, int]:
        """
        Build a generic decay directed graph.
        
        Args:
            nuclide (str): Starting or target nuclide symbol.
            digraph (DiGraph): The networkx DiGraph to build.
            direction (str): "forward" for parent->daughter, "reverse" for daughter->parent.
            decay_mode_filter (list): List of decay modes to skip.
            edge_direction (str): "forward" for parent->daughter, "reverse" for parent<-daughter.
            
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
            current = dequeue.popleft()
            generation = generations.popleft() + 1
            xpos = xpositions.popleft()
            if generation not in generation_max_xpos:
                generation_max_xpos[generation] = -1

            if traversal_direction == "forward":
                next_nuclides = self.find_daughters(current)
            else:
                next_nuclides = self.find_parents(current)

            xpos = max(xpos, generation_max_xpos[generation] + 1)
            xcounter = 0
            for next_nuclide in next_nuclides:
                try:
                    if traversal_direction == "forward":
                        parent, daughter = current, next_nuclide
                    else:
                        parent, daughter = next_nuclide, current
                    decay_mode, branching_ratio = self.decay_database.get_decay_channel(parent, daughter)
                except Exception:
                    continue

                if decay_mode_filter and decay_mode in decay_mode_filter:
                    continue
                if next_nuclide == "SF" or decay_mode == "SF":
                    continue

                if next_nuclide not in seen:
                    node_label = _parse_nuclide_label(next_nuclide)
                    if next_nuclide in self.decay_database.nuclide_index_map:
                        node_label += f'\n{self.decay_database.get_nuclide_half_life(next_nuclide, "readable")}'
                        if np.isfinite(self.decay_database.get_nuclide_half_life(next_nuclide, 's')):
                            dequeue.append(next_nuclide)
                            generations.append(generation)
                            xpositions.append(xpos + xcounter)
                    digraph.add_node(
                        next_nuclide,
                        generation=generation,
                        xpos=xpos + xcounter,
                        label=node_label,
                    )
                    seen.add(next_nuclide)
                    if xpos + xcounter > generation_max_xpos[generation]:
                        generation_max_xpos[generation] = xpos + xcounter
                    xcounter += 1

                edge_label = (
                    _parse_decay_mode_label(decay_mode)
                    + '\n'
                    + f'{branching_ratio*100:.2f}%'
                )
                if edge_direction == "forward":
                    digraph.add_edge(parent, daughter, label=edge_label)
                else:
                    digraph.add_edge(daughter, parent, label=edge_label)

        # Assign positions for plotting
        for node in digraph:
            digraph.nodes[node]["pos"] = (
                digraph.nodes[node]["xpos"],
                digraph.nodes[node]["generation"] * -1,
            )

        return digraph, max(generation_max_xpos), max(generation_max_xpos.values())


    def plot_chain_nz_chart(
        self,
        nuclide: str,
        figsize: tuple = (4, 3),
        dpi: int = 300,
        nuclei_linewidths: float = 0.5,
        **kwargs
    ):
        # digraph, _, _ = self._build_generic_decay_digraph(
        #     Nuclide(nuclide).nuclide_symbol,
        #     nx.DiGraph(),
        #     traversal_direction='forward',
        #     decay_mode_filter=['SF'],
        #     edge_direction="forward"
        # )
        
        digraph, _, _ = self._build_generic_decay_digraph(
            Nuclide(nuclide).nuclide_symbol,
            nx.DiGraph(),
            traversal_direction='reverse',
            decay_mode_filter=['β+&EC', 'SF', 'IT'],
            edge_direction="forward"
        )
        
        
        chain_nuclides = list(digraph.nodes)
        
        colormap = HalfLifeColorMap(self.decay_database.half_life_data)
        range_groups = {key: [] for key in colormap.ranges.keys()}
        stable_nuclides = []
        nz_map = {}
        for n in chain_nuclides:
            try:
                parsed = Nuclide(n)
                Z, N, state = parsed.Z, parsed.N, parsed.state
                if state:
                    continue
                nz_map[n] = (N, Z)
                half_life = self.decay_database.get_nuclide_half_life(n, units="s")
                if half_life == float('inf'):
                    stable_nuclides.append((N, Z))
                else:
                    range_key = colormap.get_range_for_half_life(half_life)
                    if range_key:
                        range_groups[range_key].append((N, Z, half_life))
            except Exception:
                continue

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        for N, Z in stable_nuclides:
            ax.add_patch(
                Rectangle(
                    (N-0.35, Z-0.35), 0.7, 0.7,
                    facecolor='black',
                    edgecolor='black',
                    linewidth=nuclei_linewidths
                )
            )
            for nuclide, (nN, nZ) in nz_map.items():
                if (nN, nZ) == (N, Z):
                    ax.text(
                        N, Z, _parse_nuclide_label(nuclide),
                        ha='center', va='center', color='white', fontsize=6, fontweight='bold', zorder=10
                    )
                    break

        for range_key, nuclides in range_groups.items():
            if not nuclides:
                continue
            cmap, norm = colormap.get_cmap_and_norm(range_key)
            for N, Z, value in nuclides:
                # color = cmap(norm(value))
                ax.add_patch(Rectangle(
                    (N-0.35, Z-0.35), 0.7, 0.7,
                    facecolor='white',
                    edgecolor='black',
                    linewidth=nuclei_linewidths,
                    # alpha=0.6
                ))

                for nuclide, (nN, nZ) in nz_map.items():
                    if (nN, nZ) == (N, Z):
                        ax.text(
                            N, Z, _parse_nuclide_label(nuclide),
                            ha='center', va='center', color='black', fontsize=6, fontweight='bold', zorder=10
                        )
                        break

        for parent, daughter in digraph.edges:
            if parent in nz_map and daughter in nz_map:
                N1, Z1 = nz_map[parent]
                N2, Z2 = nz_map[daughter]
                ax.annotate(
                    '',
                    xy=(N2, Z2), xytext=(N1, Z1),
                    arrowprops=dict(
                        arrowstyle="-|>", color="black", lw=0.5, shrinkA=6, shrinkB=6, mutation_scale=8
                    ),
                    zorder=5
                )

        all_N = [N for group in range_groups.values() for (N, Z, _) in group] + [N for (N, Z) in stable_nuclides]
        all_Z = [Z for group in range_groups.values() for (N, Z, _) in group] + [Z for (N, Z) in stable_nuclides]
        if all_N and all_Z:
            ax.set_xlim(min(all_N)-1, max(all_N)+1)
            ax.set_ylim(min(all_Z)-1, max(all_Z)+1)
            ax.xaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.set_xlabel("Neutron Number (N)")
        ax.set_ylabel("Proton Number (Z)")

        plt.tight_layout()
        return fig, ax


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
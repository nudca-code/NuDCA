# -*- coding: utf-8 -*-

from collections import deque
from typing import Any, Tuple, Dict, Optional, Union, List

import numpy as np
import networkx as nx

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors

from .nuclide import Nuclide


class DecayDiagram:

    def __init__(self, nuclide, decay_database):
        self.nuclide = Nuclide(nuclide).nuclide_symbol
        self.decay_database = decay_database
        self.nuclide_index = self.decay_database.nuclide_index_map[self.nuclide]
        

    @property
    def half_life(self, units: str = "readable") -> Union[float, str]:

        return self.decay_database.half_life(self.nuclide)


    @property
    def progeny(self) -> list[str]:

        return self.decay_database.progeny(self.nuclide)


    @property
    def branching_ratios(self) -> list[float]:
  
        return self.decay_database.branching_ratios_data[self.nuclide_index]


    @property
    def decay_modes(self) -> list[str]:

        return self.decay_database.decay_modes_data[self.nuclide_index]


    def plot(
        self,
        label_pos: float = 0.3,
        fig: Optional[matplotlib.figure.Figure] = None,
        axes: Optional[matplotlib.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ):
        digraph, max_generation, max_xpos = self._build_decay_digraph(nx.DiGraph())

        positions = nx.get_node_attributes(digraph, "pos")
        node_labels = nx.get_node_attributes(digraph, "label")
        edge_labels = nx.get_edge_attributes(digraph, "label")

        fig, axes = _check_fig_axes(
            fig, axes, figsize=(3 * max_xpos + 1.5, 3 * max_generation + 1.5)
        )

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
        digraph: nx.DiGraph
    ) -> nx.DiGraph:

        generation_max_xpos = {0: 0}
        dequeue = deque([self.nuclide])
        generations = deque([0])
        xpositions = deque([0])
        node_label = (
            _parse_nuclide_label(self.nuclide)
            + '\n' 
            + str(self.decay_database.half_life(self.nuclide, 'readable'))
        )
        
        digraph.add_node(self.nuclide, generation=0, xpos=0, label=node_label)
        seen = {self.nuclide}

        while len(dequeue) > 0:
            parent_nuclide = dequeue.popleft()
            generation = generations.popleft() + 1
            xpos = xpositions.popleft()
            if generation not in generation_max_xpos:
                generation_max_xpos[generation] = -1

            progeny = self.decay_database.progeny(parent_nuclide)
            branching_ratios = self.decay_database.branching_ratios(parent_nuclide)
            decay_modes = self.decay_database.decay_modes(parent_nuclide)

            xpos = max(xpos, generation_max_xpos[generation] + 1)
            xcounter = 0
            for idx, daughter in enumerate(progeny):
                
                # ---------
                if daughter == "SF" or decay_modes[idx] == "SF":
                    continue
                
                if daughter not in seen:
                    node_label = _parse_nuclide_label(daughter)
                    if daughter in self.decay_database.nuclide_index_map:
                        node_label += f'\n{self.decay_database.half_life(daughter, "readable")}'
                        if np.isfinite(self.decay_database.half_life(daughter, 's')):
                            dequeue.append(daughter)
                            generations.append(generation)
                            xpositions.append(xpos + xcounter)
                    
                    #----------
                    # if daughter == "SF":
                        # daughter = Nuclide(parent_nuclide).nuclide_symbol + "_SF"

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

                edge_label = (
                    _parse_decay_mode_label(decay_modes[idx])
                    + '\n'
                    + f'{branching_ratios[idx]*100:.2f}%'
                )
                digraph.add_edge(Nuclide(parent_nuclide).nuclide_symbol, daughter, label=edge_label)

        for node in digraph:
            digraph.nodes[node]["pos"] = (
                digraph.nodes[node]["xpos"],
                digraph.nodes[node]["generation"] * -1,
            )

        return digraph, max(generation_max_xpos), max(generation_max_xpos.values())








class ReverseDecayDiagram:
    """A class for visualizing reverse decay chains of nuclides.

    This class provides functionality to create and display reverse decay diagrams
    showing all possible parent nuclides that can decay to a given nuclide,
    including their branching ratios and decay modes.

    Attributes:
        nuclide (str): The nuclide symbol in the format 'Element-MassNumber'
        decay_database: The database containing decay information
        nuclide_index (int): The index of the nuclide in the database
    """

    def __init__(self, nuclide: str, decay_database: Any):
        """Initialize the ReverseDecayDiagram.

        Args:
            nuclide (str): The nuclide symbol in the format 'Element-MassNumber'
            decay_database: The database containing decay information
        """
        self.nuclide = Nuclide(nuclide).nuclide_symbol
        self.decay_database = decay_database
        self.nuclide_index = self.decay_database.nuclide_index_map[self.nuclide]

    @property
    def half_life(self, units: str = "readable") -> Union[float, str]:
        """Get the half-life of the nuclide.

        Args:
            units (str): The units for the half-life ('readable' or 's')

        Returns:
            Union[float, str]: The half-life in the specified units
        """
        return self.decay_database.half_life(self.nuclide, units)

    def find_parents(self, nuclide: str) -> List[str]:
        """Find all possible parent nuclides that can decay to the given nuclide.

        Args:
            nuclide (str): The daughter nuclide symbol

        Returns:
            List[str]: List of parent nuclides
        """
        parents = []
        nuclide = Nuclide(nuclide).nuclide_symbol
    
        for parent in self.decay_database.nuclide_index_map.keys():
            progeny = self.decay_database.progeny(parent)
            if nuclide in progeny:
                parents.append(parent)

        return parents

    def get_decay_info(self, parent: str, daughter: str) -> Optional[Tuple[str, float]]:
        """Get decay information between parent and daughter nuclides.

        Args:
            parent (str): The parent nuclide symbol
            daughter (str): The daughter nuclide symbol

        Returns:
            Optional[Tuple[str, float]]: Tuple of (decay_mode, branching_ratio) if found, None otherwise
        """
        parent = Nuclide(parent).nuclide_symbol
        daughter = Nuclide(daughter).nuclide_symbol

        if parent not in self.decay_database.nuclide_index_map:
            return None

        parent_index = self.decay_database.nuclide_index_map[parent]
        progeny = self.decay_database.progeny(parent)
        
        if daughter not in progeny:
            return None

        daughter_index = progeny.index(daughter)
        decay_mode = self.decay_database.decay_modes_data[parent_index][daughter_index]
        branching_ratio = self.decay_database.branching_ratios_data[parent_index][daughter_index]

        return decay_mode, branching_ratio

    def plot(
        self,
        label_pos: float = 0.3,
        fig: Optional[matplotlib.figure.Figure] = None,
        axes: Optional[matplotlib.axes.Axes] = None,
        kwargs_draw: Optional[Dict[str, Any]] = None,
        kwargs_edge_labels: Optional[Dict[str, Any]] = None,
    ):
        digraph, max_generation, max_xpos = self._build_reverse_decay_digraph(nx.DiGraph())

        positions = nx.get_node_attributes(digraph, "pos")
        node_labels = nx.get_node_attributes(digraph, "label")
        edge_labels = nx.get_edge_attributes(digraph, "label")

        fig, axes = _check_fig_axes(
            fig, axes, figsize=(3 * max_xpos + 1.5, 3 * max_generation + 1.5)
        )

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



    def _build_reverse_decay_digraph(
        self,
        digraph: nx.DiGraph
    ) -> Tuple[nx.DiGraph, int, int]:
        """Build the reverse decay directed graph.

        Args:
            digraph (DiGraph): The networkx DiGraph to build

        Returns:
            Tuple[DiGraph, int, int]: The graph, maximum generation, and maximum x position
        """
        generation_max_xpos = {0: 0}
        dequeue = deque([self.nuclide])
        generations = deque([0])
        xpositions = deque([0])
        node_label = (
            _parse_nuclide_label(self.nuclide)
            + '\n' 
            + str(self.decay_database.half_life(self.nuclide, 'readable'))
        )
        
        digraph.add_node(self.nuclide, generation=0, xpos=0, label=node_label)
        seen = {self.nuclide}

        while len(dequeue) > 0:
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
                        node_label += f'\n{self.decay_database.half_life(parent, "readable")}'
                        if np.isfinite(self.decay_database.half_life(parent, 's')):
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

                decay_info = self.get_decay_info(parent, daughter_nuclide)
                if decay_info:
                    decay_mode, branching_ratio = decay_info
                    edge_label = (
                        _parse_decay_mode_label(decay_mode)
                        + '\n'
                        + f'{branching_ratio*100:.2f}%'
                    )
                    digraph.add_edge(parent, daughter_nuclide, label=edge_label)

        for node in digraph:
            digraph.nodes[node]["pos"] = (
                digraph.nodes[node]["xpos"],
                digraph.nodes[node]["generation"] * -1,
            )

        return digraph, max(generation_max_xpos), max(generation_max_xpos.values()) 







def _check_fig_axes(
    fig_in: Optional[matplotlib.figure.Figure],
    axes_in: Optional[matplotlib.axes.Axes],
    **kwargs,
) -> Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """Check and create figure and axes if needed.

    Args:
        fig_in (Optional[Figure]): Input figure
        axes_in (Optional[Axes]): Input axes
        **kwargs: Additional arguments for plt.subplots

    Returns:
        Tuple[Figure, Axes]: The figure and axes to use
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
    """Parse nuclide symbol into a formatted label.

    Args:
        nuclide (str): Nuclide symbol in format 'Element-MassNumber'

    Returns:
        str: Formatted nuclide label with superscript mass number
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
    """Parse decay mode into a formatted label.

    Args:
        decay_mode (str): Decay mode symbol

    Returns:
        str: Formatted decay mode label with unicode characters
    """
    decay_mode_unicode_map = {
        "α": "\N{GREEK SMALL LETTER ALPHA}",
        "β": "\N{GREEK SMALL LETTER BETA}",
        "ε": "\N{GREEK SMALL LETTER EPSILON}",
        "+": "\N{SUPERSCRIPT PLUS SIGN}",
        "-": "\N{SUPERSCRIPT MINUS}",
        # "12C": "\N{SUPERSCRIPT ONE}\N{SUPERSCRIPT TWO}C",
        # "14C": "\N{SUPERSCRIPT ONE}\N{SUPERSCRIPT FOUR}C",
        # "20O": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT ZERO}O",
        # "23F": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT THREE}F",
        # "22Ne": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT TWO}Ne",
        # "24Ne": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT FOUR}Ne",
        # "25Ne": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT FIVE}Ne",
        # "26Ne": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT SIX}Ne",
        # "28Mg": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT EIGHT}Mg",
        # "29Mg": "\N{SUPERSCRIPT TWO}\N{SUPERSCRIPT NINE}Mg",
        # "30Mg": "\N{SUPERSCRIPT THREE}\N{SUPERSCRIPT ZERO}Mg",
        # "32Si": "\N{SUPERSCRIPT THREE}\N{SUPERSCRIPT TWO}Si",
        # "34Si": "\N{SUPERSCRIPT THREE}\N{SUPERSCRIPT FOUR}Si",
    }

    for unformatted, formatted in decay_mode_unicode_map.items():
        decay_mode = decay_mode.replace(unformatted, formatted)

    return decay_mode


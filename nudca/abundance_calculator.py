# -*- coding: utf-8 -*-

from typing import Dict, List, Tuple, Optional
import numpy as np

from .nuclide import Nuclide
from .decay_database import DecayDatabase


class AbundanceCalculator:
    """A class for calculating initial abundances based on decay chains and final abundances.

    This class helps calculate the initial abundances of parent nuclides based on
    known decay chains and the final abundance of a stable daughter nuclide.

    Attributes:
        decay_database: The database containing decay information
    """

    def __init__(self, decay_database: DecayDatabase):
        """Initialize the AbundanceCalculator.

        Args:
            decay_database: The database containing decay information
        """
        self.decay_database = decay_database

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
        # 验证输入
        for chain in decay_chains:
            if not self._validate_chain(chain, final_nuclide):
                raise ValueError(f"Invalid decay chain: {chain}")

        # 计算每条链的贡献比
        chain_contributions = self._calculate_chain_contributions(decay_chains)

        # 计算初始丰度
        initial_abundances = {}
        for chain, contribution in chain_contributions.items():
            initial_nuclide = chain[0]
            initial_abundances[initial_nuclide] = final_abundance * contribution

        return initial_abundances

    def _validate_chain(self, chain: List[str], final_nuclide: str) -> bool:
        """Validate if a decay chain is valid and leads to the final nuclide.

        Args:
            chain: List of nuclides in the decay chain
            final_nuclide: The expected final nuclide

        Returns:
            bool: True if the chain is valid, False otherwise
        """
        if not chain or chain[-1] != final_nuclide:
            return False

        for i in range(len(chain) - 1):
            parent = chain[i]
            daughter = chain[i + 1]
            
            if parent not in self.decay_database.nuclide_index_map:
                return False
                
            progeny = self.decay_database.progeny(parent)
            if daughter not in progeny:
                return False

        return True

    def _calculate_chain_contributions(self, decay_chains: List[List[str]]) -> Dict[Tuple[str, ...], float]:
        """Calculate the contribution ratio of each decay chain.

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
                
                # 获取衰变分支比
                parent_index = self.decay_database.nuclide_index_map[parent]
                progeny = self.decay_database.progeny(parent)
                daughter_index = progeny.index(daughter)
                branching_ratio = self.decay_database.branching_ratios_data[parent_index][daughter_index]
                
                contribution *= branching_ratio

            chain_contributions[tuple(chain)] = contribution
            total_contribution += contribution

        # 归一化贡献比
        for chain in chain_contributions:
            chain_contributions[chain] /= total_contribution

        return chain_contributions


def calculate_initial_abundances(
    decay_database: DecayDatabase,
    decay_chains: List[List[str]],
    final_nuclide: str,
    final_abundance: float
) -> Dict[str, float]:
    """Helper function to calculate initial abundances.

    Args:
        decay_database: The database containing decay information
        decay_chains: List of decay chains
        final_nuclide: The final stable nuclide
        final_abundance: The abundance of the final nuclide

    Returns:
        Dict[str, float]: Dictionary mapping initial nuclides to their abundances
    """
    calculator = AbundanceCalculator(decay_database)
    return calculator.calculate_initial_abundances(decay_chains, final_nuclide, final_abundance) 
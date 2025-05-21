#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
NuDCA: A Numerical Code for Nuclear Decay Chains in Astrophysics

This package provides tools for analyzing nuclear decay chains and simulating
kilonovae phenomena in astrophysics. It includes comprehensive nuclear decay
databases, efficient matrix calculations for decay chains, and visualization
tools for nuclear charts.

Main Features:
-------------
- Nuclear decay database with comprehensive nuclide information
- Efficient sparse matrix operations for decay chains
- Nuclear chart visualization
- Kilonovae analysis tools
- Support for ENDF-B-VIII.1 decay data

Example:
--------
>>> from nudca import DecayDatabase, DecayMatrix, RadioactiveDecay
>>> db = DecayDatabase('ENDF-B-VIII.1_decay')
>>> half_life = db.half_life('U-238')
>>> matrix = DecayMatrix.from_database(db)
>>> decay = RadioactiveDecay({'U-238': 1.0}, db, matrix)
>>> heating_rates = decay.decay_heating_rates([1e6, 1e7, 1e8])
"""

__version__ = '0.1.0'
__author__ = 'Qiuhong-Chen'
__email__ = 'chohonche@163.com'
__license__ = 'MIT'
__url__ = 'https://github.com/chohonche/nudca'

# Core classes and functions
from .decay_database import (
    DecayDatabase,
    DecayDatabaseManager,
    load_decay_database
)
from .decay_matrix import (
    DecayMatrix,
    MatrixBuilder,
    load_decay_matrix
)
from .nuclide import (
    Nuclide,
    NuclideStrError
)
from .nucleo_decay import (
    RadioactiveDecay,
    RadioactiveDecayBase
)

# Constants and utilities
from .constants import (
    NA_CGS,
    EV_CGS,
    # Add other constants as needed
)

# Import specific functions from io module
# from .io import (
#     load_nuclide_data,
#     save_nuclide_data
# )

# Import visualization tools
from .decay_diagram import (
    DecayDiagram,
    # plot_nuclear_chart
)

# Kilonovae module
from . import kilonovae

# Make these available at package level
__all__ = [
    # Core classes
    'DecayDatabase',
    'DecayDiagram',
    'DecayDatabaseManager',
    'DecayMatrix',
    'MatrixBuilder',
    'Nuclide',
    'NuclideStrError',
    'RadioactiveDecay',
    'RadioactiveDecayBase',
    
    # Functions
    'load_decay_database',
    'load_decay_matrix',
    # 'load_nuclide_data',
    # 'save_nuclide_data',
    # 'plot_decay_chain',
    # 'plot_nuclear_chart',
    
    # Submodules
    'kilonovae',
]

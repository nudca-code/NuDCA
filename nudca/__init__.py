#!/usr/bin/env python
# -*- coding: utf-8 -*-


__version__ = '0.1.0'
__author__ = 'Qiuhong-Chen'
__email__ = 'chohonche@163.com'
__license__ = 'MIT'
__url__ = 'https://github.com/nudca-code/NuDCA'


# Core classes and functions
from .nuclide import (
    Nuclide,
    NuclideStrError
)
from .decay_network import (
    DecayDatabase,
    DecayMatrix,
    MatrixBuilder,
    DecayDiagram
)
from .decay_core import (
    RadioactiveDecay
)
from .utils import (
    DecayDatabaseManager
)
from .io import (
    load_decay_database,
    load_decay_matrix,
    serialize_decay_database,
    serialize_decay_matrix,
    Inputer,
    Outputer
)

# Kilonovae module
from . import kilonovae

# Make these available at package level
__all__ = [
    # Core classes
    'Nuclide',
    'NuclideStrError',
    'DecayDatabase',
    'DecayDatabaseManager',
    'MatrixBuilder',
    'DecayMatrix',
    'DecayDiagram',
    'RadioactiveDecay',
    # Functions
    'load_decay_database',
    'load_decay_matrix',
    'serialize_decay_database',
    'serialize_decay_matrix',
    'Inputer',
    'Outputer',
    # Submodules
    'kilonovae',
]

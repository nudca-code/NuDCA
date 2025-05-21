# NuDCA API Reference

This document provides detailed information about the NuDCA (Nuclear Decay Chains in Astrophysics) API.

## Table of Contents

- [Core Classes](#core-classes)
  - [DecayDatabase](#decaydatabase)
  - [DecayMatrix](#decaymatrix)
  - [MatrixBuilder](#matrixbuilder)
  - [Nuclide](#nuclide)
- [Utility Functions](#utility-functions)
- [Data Types](#data-types)

## Core Classes

### DecayDatabase

The main class for accessing nuclear decay data.

```python
class DecayDatabase:
    def __init__(
        self,
        data_source: str,
        nuclides: np.ndarray,
        half_life_data: np.ndarray,
        decay_modes_data: np.ndarray,
        progeny_data: np.ndarray,
        branching_ratios_data: np.ndarray,
        decay_energies_data: np.ndarray
    ) -> None
```

#### Methods

##### half_life
```python
def half_life(self, nuclide: str, units: str = 'readable') -> Union[float, str]
```

Get the half-life of a nuclide.

- **Parameters:**
  - `nuclide`: Nuclide symbol (e.g., 'U-238')
  - `units`: Output units ('readable', 's', or default)
- **Returns:**
  - Half-life in specified units

##### progeny
```python
def progeny(self, nuclide: str) -> list[str]
```
Get the progeny nuclides of a given nuclide.

- **Parameters:**
  - `nuclide`: Nuclide symbol
- **Returns:**
  - List of progeny nuclide symbols

##### decay_modes
```python
def decay_modes(self, nuclide: str) -> List
```
Get the decay modes of a nuclide.

- **Parameters:**
  - `nuclide`: Nuclide symbol
- **Returns:**
  - List of decay modes

##### decay_energy
```python
def decay_energy(self, nuclide: str, energy_type: str) -> float
```
Get specific decay energy by type.

- **Parameters:**
  - `nuclide`: Nuclide symbol
  - `energy_type`: Type of energy ('EM', 'LP', 'HP', 'Neutrino', etc.)
- **Returns:**
  - Decay energy in MeV

### DecayMatrix

Class for handling decay matrix calculations.

```python
class DecayMatrix:
    def __init__(
        self,
        decay_constants: np.ndarray,
        matrix_P: sparse.csc_matrix,
        matrix_P_inv: sparse.csc_matrix
    ) -> None
```

#### Properties

- `decay_constants`: Array of decay constants
- `matrix_P`: Decay matrix P
- `matrix_P_inv`: Inverse of matrix P
- `initial_abundance`: Initial abundance vector
- `matrix_Lambda`: Lambda matrix

### MatrixBuilder

Class for building decay matrices.

```python
class MatrixBuilder:
    def __init__(self, decay_database: DecayDatabase) -> None
```

#### Methods

##### build_decay_matrix
```python
def build_decay_matrix(self) -> Tuple[sparse.csc_matrix, sparse.csc_matrix]
```
Constructs the decay matrix and its inverse.

- **Returns:**
  - Tuple of (matrix_P, matrix_P_inv)

##### build_matrix_A
```python
def build_matrix_A(self) -> sparse.csc_matrix
```
Builds the decay matrix A.

- **Returns:**
  - Sparse matrix A

### Nuclide

Class for handling nuclide information.

```python
class Nuclide:
    def __init__(self, nuclide: Union[str, int]) -> None
```

#### Properties

- `Z`: Atomic number
- `A`: Mass number
- `N`: Neutron number
- `state`: Metastable state
- `id`: Canonical nuclide ID

## Utility Functions

### load_decay_database
```python
def load_decay_database(data_source: str = 'ENDF-B-VIII.1_decay') -> DecayDatabase
```
Load decay database from specified source.

### load_decay_matrix
```python
def load_decay_matrix(data_source: str = 'ENDF-B-VIII.1_decay') -> DecayMatrix
```
Load decay matrix from specified source.

## Data Types

### Energy Types
- `EM`: Electromagnetic energy
- `LP`: Light particle energy
- `HP`: Heavy particle energy
- `Neutrino`: Neutrino energy
- `Gamma`: Gamma ray energy
- `Beta_Minus`: Beta minus decay energy
- `Beta_Plus`: Beta plus decay energy
- `Alpha`: Alpha decay energy
- `Neutron`: Neutron emission energy
- `Proton`: Proton emission energy
- `Effective_Q`: Effective Q-value

### Decay Modes
- `β-`: Beta minus decay
- `β+`: Beta plus decay
- `α`: Alpha decay
- `EC`: Electron capture
- `IT`: Isomeric transition
- `SF`: Spontaneous fission
- `n`: Neutron emission
- `p`: Proton emission
- `Stable`: Stable nuclide

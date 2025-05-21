# NuDCA: A Numerical Code for Nuclear Decay Chains in Astrophysics

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen)](https://nudca.readthedocs.io/)

NuDCA is a powerful scientific computing library for nuclear decay and kilonovae research. It provides comprehensive tools for analyzing radioactive decay chains, calculating decay matrices, and simulating kilonovae phenomena.

## Features

- **Nuclear Decay Database**
  - Comprehensive database of nuclides and their properties
  - Support for ENDF-B-VIII.1 decay data
  - Half-life, decay modes, branching ratios, and decay energies
  - Nuclear chart visualization

- **Decay Matrix Calculations**
  - Efficient sparse matrix operations for decay chains
  - Support for complex decay networks
  - High-precision numerical calculations
  - Optimized for large-scale simulations

- **Kilonovae Analysis**
  - Simulation of r-process nucleosynthesis
  - Analysis of kilonovae light curves
  - Energy deposition calculations
  - Nuclear heating rate computations

## Installation

```bash
# Using pip
pip install nudca

# Using poetry
poetry add nudca
```

## Quick Start

```python
from nudca import DecayDatabase, DecayMatrix

# Load decay database
db = DecayDatabase('ENDF-B-VIII.1_decay')

# Get nuclide information
half_life = db.half_life('U-238')
progeny = db.progeny('U-238')
decay_energy = db.decay_energy('U-238', 'Alpha')

# Create decay matrix
matrix = DecayMatrix.from_database(db)

# Plot nuclear chart
db.plot_nuclear_chart(
    min_lifetime=1.0,
    max_lifetime=1e10,
    cmap='viridis'
)
```

## Documentation

For detailed documentation, please visit our [documentation site](https://nudca.readthedocs.io/).

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## Citation

If you use NuDCA in your research, please cite:

```bibtex
@software{nudca,
  author = {Qiuhong-Chen},
  title = {NuDCA: A Numerical Code for Nuclear Decay Chains in Astrophysics},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/yourusername/nudca}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- ENDF/B-VIII.1 nuclear data library
- SciPy and NumPy communities
- All contributors and users of this project

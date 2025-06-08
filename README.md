# NuDCA: A Numerical Code for Nuclear Decay Chains in Astrophysics

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen)](https://nudca-code.github.io/nudca.github.io/)

`NuDCA` is an open-source Python program specifically designed for calculating the decay chains of radioactive nuclides. We have developed a versatile computational framework that efficiently handles decay chains with complex branching ratios. This framework employs a matrix-based numerical solving method, ensuring both analytical accuracy and computational efficiency, enabling rapid resolution of decay chain differential equations.

The primary goal of `NuDCA` is to provide researchers with a flexible and reliable tool for simulating the accumulation, decay, and associated physical phenomena of radioactive nuclides. The program supports a wide range of complex scenarios, including the splitting and merging of multi-branch decay chains and the superposition of various decay modes. This flexibility makes it invaluable across fields such as nuclear science, geochronology, and astrophysics. As an open-source Python-based software, `NuDCA` offers excellent extensibility, allowing researchers to easily customize its functionality to meet specific requirements.

In extreme astrophysical environments like core-collapse supernovae (CCSN) and neutron star mergers (NSM), significant energy is released through radioactive decay, powering transient phenomena such as supernovae (SN) and kilonovae (KN). To address the needs of this research area, `NuDCA` features a specialized module for calculating the light curves of kilonovae. This capability provides critical support for studying the synthesis of r-process elements, the physical properties of ejecta, and the radiation characteristics of kilonovae.

Furthermore, the development of `NuDCA` emphasizes user experience and community collaboration. With its intuitive interface and comprehensive documentation, even beginners can quickly start using the program. As an open-source project, `NuDCA` fosters collaboration between academia and industry, aiming to drive technological advancements and scientific research, establishing itself as a benchmark tool in the field of radioactive decay simulation.

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
pip install nudca
```

## Documentation

For detailed documentation, please visit our [documentation site](https://nudca-code.github.io/nudca.github.io/).



## Citation

If you use NuDCA in your research, please cite:

```bibtex
@software{nudca,
  author = {Qiuhong-Chen and Menghua-Chen},
  title = {NuDCA: A Numerical Code for Nuclear Decay Chains in Astrophysics},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/nudca-code/NuDCA.git}
}
```

## Acknowledgments

- ENDF/B-VIII.1 nuclear data library
- SciPy and NumPy communities
- All contributors and users of this project

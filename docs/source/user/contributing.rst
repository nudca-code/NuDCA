Contributor's guide
===================

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin: 2em 0;">
       <h2 style="color: white; margin-top: 0;">🤝 Contributing to NuDCA</h2>
       <p style="font-size: 1.1em; margin-bottom: 0;">Join the community and help advance nuclear decay calculations and kilonova modeling in Python</p>
   </div>

Welcome to the NuDCA contributor's guide! We appreciate your interest in contributing to this open-source project. Whether you're fixing a bug, adding a feature, improving documentation, or suggesting enhancements, your contributions are valuable to the nuclear astrophysics community.

.. raw:: html

   <div style="background: #e6f3ff; padding: 1.5em; border-radius: 8px; margin: 2em 0;">
       <h3 style="color: #0066cc; margin-top: 0;">🎯 Quick Start</h3>
       <ol style="color: #2c5282;">
           <li>Fork the repository on GitHub</li>
           <li>Clone your fork locally</li>
           <li>Create a feature branch</li>
           <li>Make your changes</li>
           <li>Run tests and ensure they pass</li>
           <li>Submit a pull request</li>
       </ol>
   </div>

Types of Contributions
======================

We welcome all types of contributions:

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1.5em; margin: 2em 0;">
       <div style="background: #f0fff4; padding: 1.5em; border-radius: 8px; border-left: 4px solid #48bb78;">
           <h4 style="color: #22543d; margin-top: 0;">🐛 Bug Reports</h4>
           <p style="color: #2f855a;">Help us identify and fix issues</p>
       </div>
       
       <div style="background: #fef5e7; padding: 1.5em; border-radius: 8px; border-left: 4px solid #ed8936;">
           <h4 style="color: #c05621; margin-top: 0;">✨ Feature Requests</h4>
           <p style="color: #dd6b20;">Suggest new functionality</p>
       </div>
       
       <div style="background: #f0f8ff; padding: 1.5em; border-radius: 8px; border-left: 4px solid #4299e1;">
           <h4 style="color: #2b6cb0; margin-top: 0;">📝 Documentation</h4>
           <p style="color: #2c5282;">Improve guides and examples</p>
       </div>
       
       <div style="background: #faf5ff; padding: 1.5em; border-radius: 8px; border-left: 4px solid #9f7aea;">
           <h4 style="color: #553c9a; margin-top: 0;">🔬 Code</h4>
           <p style="color: #6b46c1;">Bug fixes and new features</p>
       </div>
   </div>

Getting Started
===============

Development Environment Setup
-----------------------------

1. **Fork and Clone**

   .. code-block:: bash

      # Fork the repository on GitHub
      # Then clone your fork
      git clone https://github.com/YOUR_USERNAME/NuDCA.git
      cd NuDCA
      
      # Add upstream remote
      git remote add upstream https://github.com/nudca-code/NuDCA.git

2. **Create Virtual Environment**

   .. code-block:: bash

      python -m venv nudca-dev
      source nudca-dev/bin/activate  # Linux/macOS
      # or
      nudca-dev\Scripts\activate     # Windows

3. **Install Development Dependencies**

   .. code-block:: bash

      pip install -e ".[dev]"

4. **Verify Setup**

   .. code-block:: bash

      python -c "import nudca; print('Setup successful!')"
      pytest tests/ -v

Branch Management
-----------------

**Create Feature Branch:**

.. code-block:: bash

   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/issue-description

**Keep Your Fork Updated:**

.. code-block:: bash

   git fetch upstream
   git checkout main
   git merge upstream/main
   git push origin main

Contributing Workflow
=====================

Code Contributions
------------------

**1. Choose an Issue**

* Browse `GitHub Issues <https://github.com/nudca-code/NuDCA/issues>`_
* Look for "good first issue" labels for beginners
* Comment on issues you'd like to work on

**2. Development Process**

.. code-block:: bash

   # Create and switch to feature branch
   git checkout -b feature/my-feature
   
   # Make your changes
   # Add/modify code, tests, documentation
   
   # Run tests frequently
   pytest tests/
   
   # Check code style
   black nudca/ tests/
   flake8 nudca/ tests/

**3. Commit Guidelines**

Follow conventional commit format:

.. code-block:: bash

   git commit -m "feat: add support for custom nuclear data import"
   git commit -m "fix: resolve matrix computation issue for large networks"
   git commit -m "docs: update installation instructions"
   git commit -m "test: add unit tests for decay calculations"

**Commit Types:**
* ``feat``: New features
* ``fix``: Bug fixes  
* ``docs``: Documentation changes
* ``test``: Test additions/modifications
* ``refactor``: Code refactoring
* ``perf``: Performance improvements

**4. Push and Create Pull Request**

.. code-block:: bash

   git push origin feature/my-feature

Then create a pull request on GitHub with:
* Clear title and description
* Reference to related issues
* Screenshots for UI changes
* Test results

Code Standards
==============

Style Guidelines
----------------

NuDCA follows Python community standards:

**PEP 8 Compliance:**

.. code-block:: bash

   # Auto-format with black
   black nudca/ tests/
   
   # Check style with flake8
   flake8 nudca/ tests/

**Type Hints:**

.. code-block:: python

   from typing import List, Dict, Optional, Union
   
   def calculate_decay_rates(
       nuclides: List[str],
       times: np.ndarray,
       database: DecayDatabase
   ) -> Dict[str, np.ndarray]:
       """Calculate decay rates for given nuclides."""
       ...

**Docstrings (NumPy Style):**

.. code-block:: python

   def decay_process(
       self, 
       times: np.ndarray, 
       method: str = 'matrix_exponential'
   ) -> Dict[str, np.ndarray]:
       """
       Calculate nuclide abundances over time.
       
       Parameters
       ----------
       times : np.ndarray
           Time points at which to calculate abundances.
       method : str, optional
           Calculation method to use, by default 'matrix_exponential'.
           
       Returns
       -------
       Dict[str, np.ndarray]
           Dictionary mapping nuclide names to abundance arrays.
           
       Examples
       --------
       >>> calc = RadioactiveDecay({'U238': 1.0}, db, matrix)
       >>> abundances = calc.decay_process([0, 1e6, 1e9])
       """

Testing Standards
-----------------

**Test Structure:**

.. code-block:: python

   import pytest
   import numpy as np
   from nudca import RadioactiveDecay, load_decay_database
   
   class TestRadioactiveDecay:
       """Test suite for RadioactiveDecay class."""
       
       @pytest.fixture
       def setup_database(self):
           """Fixture to load test database."""
           return load_decay_database()
       
       def test_single_isotope_decay(self, setup_database):
           """Test decay of single isotope."""
           calc = RadioactiveDecay({'Co60': 1.0}, setup_database)
           result = calc.decay_process([0, 1e6])
           
           assert len(result) > 0
           assert 'Co60' in result
           assert result['Co60'][0] == 1.0  # Initial condition
           assert result['Co60'][1] < 1.0   # Decay occurred
       
       @pytest.mark.parametrize("isotope,expected_products", [
           ('U238', ['Th234', 'Pa234', 'U234']),
           ('Co60', ['Ni60']),
       ])
       def test_decay_products(self, setup_database, isotope, expected_products):
           """Test that correct decay products are generated."""
           calc = RadioactiveDecay({isotope: 1.0}, setup_database)
           result = calc.decay_process([0, 1e12])
           
           for product in expected_products:
               assert product in result

**Test Coverage:**

.. code-block:: bash

   # Run tests with coverage
   pytest --cov=nudca --cov-report=html tests/
   
   # View coverage report
   open htmlcov/index.html

**Performance Tests:**

.. code-block:: python

   def test_large_network_performance(benchmark):
       """Benchmark performance with large decay networks."""
       db = load_decay_database()
       matrix = load_decay_matrix()
       
       # Create calculation with many isotopes
       initial = {f'U{A}': 1.0 for A in range(230, 240)}
       calc = RadioactiveDecay(initial, db, matrix)
       
       def run_calculation():
           return calc.decay_process([0, 1e6, 1e9])
       
       result = benchmark(run_calculation)
       assert len(result) > 0

Documentation Standards
=======================

Documentation Types
-------------------

**API Documentation:**
* Comprehensive docstrings for all public functions
* Parameter and return type documentation
* Usage examples

**User Guide:**
* Conceptual explanations
* Tutorial-style content
* Best practices

**Examples:**
* Complete working examples
* Jupyter notebooks
* Real-world applications

**Building Documentation:**

.. code-block:: bash

   cd docs/
   make html
   
   # Live rebuild during development
   pip install sphinx-autobuild
   sphinx-autobuild source build/html

Writing Guidelines
------------------

**Clear and Concise:**
* Use simple, direct language
* Explain scientific concepts clearly
* Provide context for astrophysical applications

**Examples and Code:**
* Include runnable code examples
* Show expected outputs
* Cover common use cases

**Cross-References:**
* Link to related functions and concepts
* Reference external papers and data sources
* Connect to broader scientific context

Review Process
==============

Pull Request Guidelines
-----------------------

**Before Submitting:**

.. code-block:: bash

   # Ensure tests pass
   pytest tests/
   
   # Check code style
   black --check nudca/ tests/
   flake8 nudca/ tests/
   
   # Update documentation if needed
   cd docs/ && make html

**PR Description Template:**

.. code-block:: markdown

   ## Description
   Brief description of changes
   
   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Documentation update
   - [ ] Performance improvement
   
   ## Testing
   - [ ] Tests added/updated
   - [ ] All tests pass
   - [ ] Documentation updated
   
   ## Related Issues
   Closes #123

**Review Process:**

1. **Automated Checks**: CI/CD runs tests and style checks
2. **Code Review**: Maintainers review for quality and design
3. **Testing**: Manual testing for complex features
4. **Documentation**: Ensure documentation is complete
5. **Approval**: At least one maintainer approval required

Community Guidelines
====================

Code of Conduct
---------------

We are committed to providing a welcoming and inclusive environment for all contributors:

* **Be respectful** in all interactions
* **Be collaborative** and help newcomers
* **Be constructive** in feedback and criticism
* **Focus on the science** and technical merit
* **Welcome diverse perspectives** and backgrounds

Communication Channels
----------------------

* **GitHub Issues**: Bug reports and feature requests
* **GitHub Discussions**: Questions and general discussion
* **Email**: chohonche@163.com for private matters
* **Documentation**: This guide and the user manual

Recognition
-----------

Contributors are recognized through:

* **GitHub Contributor Graph**: Automatic tracking of contributions
* **Acknowledgments**: Major contributors listed in documentation
* **Release Notes**: Contributors mentioned in version releases
* **CITATION.cff**: Academic citation tracking

Specific Contribution Areas
===========================

Nuclear Data
------------

**Data Sources:**
* ENDF-B-VIII.1 evaluations
* Experimental measurements
* Theoretical calculations

**Validation:**
* Cross-checks between data sources
* Comparison with experimental results
* Energy conservation verification

**Format Standards:**
* JSON for structured data
* CSV for tabular data
* Documentation of data provenance

.. code-block:: python

   # Example: Adding new nuclear data
   def add_custom_decay_data(database, nuclide_data):
       """Add custom decay data to database."""
       # Validate data format
       validate_decay_data(nuclide_data)
       
       # Check energy conservation
       verify_energy_conservation(nuclide_data)
       
       # Add to database
       database.add_nuclide(nuclide_data)

Algorithms and Methods
----------------------

**Numerical Methods:**
* Matrix exponential algorithms
* Time stepping schemes
* Convergence criteria

**Performance Optimization:**
* Sparse matrix operations
* Vectorization strategies
* Memory management

**New Physics:**
* Modified decay rates in extreme environments
* Fission fragment distributions
* Custom thermalization models

Visualization and Analysis
--------------------------

**Plotting Functions:**
* Nuclear chart visualizations
* Decay chain diagrams
* Abundance evolution plots

**Analysis Tools:**
* Statistical analysis functions
* Uncertainty quantification
* Parameter sensitivity studies

**Interactive Features:**
* Jupyter notebook widgets
* Real-time parameter adjustment
* 3D visualizations

Release Process
===============

Version Planning
----------------

**Semantic Versioning:**
* **Major (X.0.0)**: Breaking changes
* **Minor (X.Y.0)**: New features (backward compatible)
* **Patch (X.Y.Z)**: Bug fixes

**Release Schedule:**
* Major releases: Annual
* Minor releases: Quarterly
* Patch releases: As needed

**Feature Planning:**
* Community input through GitHub Issues
* Roadmap discussions in GitHub Discussions
* Scientific priorities from user feedback

Release Workflow
-----------------

**Pre-release:**
1. Feature freeze
2. Beta testing period
3. Documentation review
4. Performance benchmarking

**Release:**
1. Version number update
2. Changelog generation
3. Tag creation
4. Package building and upload
5. GitHub release notes

**Post-release:**
1. Announcement to community
2. Update citations and references
3. Begin next release planning

Getting Help
============

For Contributors
----------------

* **Technical Questions**: GitHub Discussions
* **Bug Reports**: GitHub Issues
* **Design Discussions**: GitHub Issues with "enhancement" label
* **Direct Contact**: chohonche@163.com

For Users
---------

* **Usage Questions**: GitHub Discussions
* **Documentation**: :doc:`user_guide` and :doc:`api`
* **Examples**: :doc:`examples` and Jupyter notebooks
* **Installation Help**: :doc:`installation`

Learning Resources
==================

Nuclear Physics Background
--------------------------

* **Textbooks**: 
  * "Introduction to Nuclear and Particle Physics" by Das & Ferbel
  * "Nuclear Physics: Principles and Applications" by Lilley

* **Review Papers**:
  * r-process nucleosynthesis reviews
  * Kilonova modeling papers
  * Nuclear data evaluation methods

Python and Scientific Computing
-------------------------------

* **NumPy Documentation**: https://numpy.org/doc/
* **SciPy Lecture Notes**: https://scipy-lectures.org/
* **Scientific Python Ecosystem**: https://scientific-python.org/

Contributing to Open Source
---------------------------

* **GitHub Guides**: https://guides.github.com/
* **Open Source Guides**: https://opensource.guide/
* **Python Developer's Guide**: https://devguide.python.org/

.. raw:: html

   <div style="background: #d4edda; padding: 2em; border-radius: 8px; margin: 3em 0; text-align: center;">
       <h3 style="color: #155724; margin-top: 0;">🎉 Thank You!</h3>
       <p style="color: #155724; font-size: 1.1em; margin-bottom: 0;">
           Every contribution, no matter how small, helps advance nuclear astrophysics research and education. 
           Thank you for being part of the NuDCA community!
       </p>
   </div> 
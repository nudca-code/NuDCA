Building from source
====================

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin: 2em 0;">
       <h2 style="color: white; margin-top: 0;">🔨 Building NuDCA from Source</h2>
       <p style="font-size: 1.1em; margin-bottom: 0;">Complete guide for developers to build and install NuDCA from source code</p>
   </div>

This document provides detailed instructions for building NuDCA from source. This is primarily intended for developers who want to contribute to NuDCA or need the latest development version.

Prerequisites
=============

System Requirements
-------------------

**Operating Systems:**
* Linux (Ubuntu 18.04+, CentOS 7+, etc.)
* macOS (10.14+)
* Windows (10+ with WSL2 recommended)

**Python Requirements:**
* Python 3.8 or higher
* pip (latest version recommended)

**Build Tools:**

For Linux/macOS:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get update
   sudo apt-get install git build-essential python3-dev
   
   # CentOS/RHEL
   sudo yum install git gcc gcc-c++ python3-devel
   
   # macOS (with Homebrew)
   brew install git python

For Windows:

.. code-block:: powershell

   # Install Visual Studio Build Tools
   # Download from: https://visualstudio.microsoft.com/downloads/
   
   # Install Git for Windows
   # Download from: https://git-scm.com/download/win

Python Dependencies
-------------------

**Required packages:**

.. code-block:: bash

   pip install numpy scipy matplotlib astropy pandas

**Optional packages for development:**

.. code-block:: bash

   pip install pytest pytest-cov black flake8 mypy sphinx

Getting the Source Code
=======================

Clone from GitHub
------------------

.. code-block:: bash

   git clone https://github.com/nudca-code/NuDCA.git
   cd NuDCA

For development, fork the repository first:

.. code-block:: bash

   # Fork on GitHub, then clone your fork
   git clone https://github.com/YOUR_USERNAME/NuDCA.git
   cd NuDCA
   
   # Add upstream remote
   git remote add upstream https://github.com/nudca-code/NuDCA.git

Development Installation
========================

Virtual Environment Setup
--------------------------

**Using venv (recommended):**

.. code-block:: bash

   # Create virtual environment
   python -m venv nudca-dev
   
   # Activate (Linux/macOS)
   source nudca-dev/bin/activate
   
   # Activate (Windows)
   nudca-dev\Scripts\activate

**Using conda:**

.. code-block:: bash

   conda create -n nudca-dev python=3.9
   conda activate nudca-dev

Editable Installation
---------------------

Install NuDCA in development mode:

.. code-block:: bash

   # Install in editable mode with development dependencies
   pip install -e ".[dev]"
   
   # Alternative: install dependencies separately
   pip install -e .
   pip install -r requirements-dev.txt

This creates an "editable" installation where changes to the source code are immediately reflected without reinstalling.

Verify Installation
-------------------

.. code-block:: python

   import nudca
   print(nudca.__version__)
   
   # Run basic functionality test
   db = nudca.load_decay_database()
   print(f"Loaded {len(db.nuclides)} nuclides")

Building Documentation
======================

Prerequisites
-------------

.. code-block:: bash

   pip install sphinx sphinx-design sphinx-copybutton pydata-sphinx-theme myst-parser

Building HTML Documentation
---------------------------

.. code-block:: bash

   cd docs
   make html
   
   # Open in browser (Linux/macOS)
   open build/html/index.html
   
   # Windows
   start build/html/index.html

Building Other Formats
-----------------------

.. code-block:: bash

   # PDF (requires LaTeX)
   make latexpdf
   
   # EPUB
   make epub
   
   # Clean build files
   make clean

Live Documentation Server
-------------------------

For continuous development:

.. code-block:: bash

   # Install sphinx-autobuild
   pip install sphinx-autobuild
   
   # Start live server
   sphinx-autobuild docs/source docs/build/html

Testing
=======

Running Tests
-------------

**Basic test suite:**

.. code-block:: bash

   # Run all tests
   pytest
   
   # Run with coverage
   pytest --cov=nudca
   
   # Run specific test file
   pytest tests/test_decay_network.py
   
   # Run specific test
   pytest tests/test_decay_network.py::test_decay_database_creation

**Parallel execution:**

.. code-block:: bash

   # Install pytest-xdist
   pip install pytest-xdist
   
   # Run tests in parallel
   pytest -n auto

Test Categories
---------------

**Unit tests:**

.. code-block:: bash

   pytest tests/unit/

**Integration tests:**

.. code-block:: bash

   pytest tests/integration/

**Performance tests:**

.. code-block:: bash

   pytest tests/performance/ --benchmark-only

Writing Tests
-------------

**Test file structure:**

.. code-block:: python

   # tests/test_new_feature.py
   import pytest
   import numpy as np
   import nudca
   
   class TestNewFeature:
       def setup_method(self):
           """Setup for each test method."""
           self.decay_db = nudca.load_decay_database()
       
       def test_basic_functionality(self):
           """Test basic functionality."""
           result = new_feature_function()
           assert result is not None
       
       def test_edge_cases(self):
           """Test edge cases."""
           with pytest.raises(ValueError):
               invalid_input_function()

Code Quality
============

Linting and Formatting
----------------------

**Black (code formatting):**

.. code-block:: bash

   # Format all Python files
   black nudca/ tests/
   
   # Check what would be formatted
   black --check nudca/ tests/

**Flake8 (linting):**

.. code-block:: bash

   # Run linter
   flake8 nudca/ tests/
   
   # With specific configuration
   flake8 --config=setup.cfg nudca/

**MyPy (type checking):**

.. code-block:: bash

   # Type check
   mypy nudca/
   
   # Generate type coverage report
   mypy --html-report mypy-report nudca/

Pre-commit Hooks
----------------

Set up pre-commit hooks for automated code quality checks:

.. code-block:: bash

   # Install pre-commit
   pip install pre-commit
   
   # Install hooks
   pre-commit install
   
   # Run on all files
   pre-commit run --all-files

Configuration in `.pre-commit-config.yaml`:

.. code-block:: yaml

   repos:
     - repo: https://github.com/psf/black
       rev: 22.3.0
       hooks:
         - id: black
     - repo: https://github.com/pycqa/flake8
       rev: 4.0.1
       hooks:
         - id: flake8
     - repo: https://github.com/pre-commit/mirrors-mypy
       rev: v0.942
       hooks:
         - id: mypy

Performance Profiling
=====================

Profiling Tools
---------------

**cProfile:**

.. code-block:: python

   import cProfile
   import nudca
   
   def profile_decay_calculation():
       db = nudca.load_decay_database()
       matrix = nudca.load_decay_matrix()
       calc = nudca.RadioactiveDecay({'U238': 1.0}, db, matrix)
       calc.decay_process([1e6, 1e9, 1e12])
   
   cProfile.run('profile_decay_calculation()', 'profile_output.prof')

**line_profiler:**

.. code-block:: bash

   # Install
   pip install line_profiler
   
   # Add @profile decorator to functions
   # Run profiler
   kernprof -l -v script.py

**memory_profiler:**

.. code-block:: bash

   # Install
   pip install memory_profiler
   
   # Profile memory usage
   mprof run script.py
   mprof plot

Benchmarking
------------

**pytest-benchmark:**

.. code-block:: python

   def test_decay_calculation_performance(benchmark):
       db = nudca.load_decay_database()
       matrix = nudca.load_decay_matrix()
       
       def decay_calculation():
           calc = nudca.RadioactiveDecay({'U238': 1.0}, db, matrix)
           return calc.decay_process([1e6, 1e9, 1e12])
       
       result = benchmark(decay_calculation)
       assert result is not None

Run benchmarks:

.. code-block:: bash

   pytest tests/performance/ --benchmark-only

Packaging and Distribution
==========================

Building Packages
-----------------

**Source distribution:**

.. code-block:: bash

   python setup.py sdist

**Wheel distribution:**

.. code-block:: bash

   pip install wheel
   python setup.py bdist_wheel

**Using build (recommended):**

.. code-block:: bash

   pip install build
   python -m build

Local Installation Testing
--------------------------

.. code-block:: bash

   # Install from local wheel
   pip install dist/nudca-*.whl
   
   # Test in clean environment
   python -c "import nudca; print(nudca.__version__)"

Release Process
===============

Version Management
------------------

NuDCA uses semantic versioning (MAJOR.MINOR.PATCH):

.. code-block:: bash

   # Update version in setup.py and nudca/__init__.py
   # Create git tag
   git tag -a v0.2.0 -m "Release version 0.2.0"
   git push origin v0.2.0

Release Checklist
-----------------

1. **Update version numbers**
2. **Update CHANGELOG.md**
3. **Run full test suite**
4. **Build documentation**
5. **Create release branch**
6. **Tag release**
7. **Build and upload packages**
8. **Update GitHub release notes**

Continuous Integration
======================

GitHub Actions
--------------

NuDCA uses GitHub Actions for CI/CD. The workflow includes:

* **Testing** on multiple Python versions and operating systems
* **Code quality** checks (linting, formatting)
* **Documentation** building
* **Package** building and testing

Local CI Simulation
-------------------

.. code-block:: bash

   # Install act (GitHub Actions local runner)
   # Run CI locally
   act

Troubleshooting
===============

Common Build Issues
-------------------

**Import errors:**

.. code-block:: bash

   # Ensure proper Python path
   export PYTHONPATH="${PYTHONPATH}:$(pwd)"

**Missing dependencies:**

.. code-block:: bash

   # Install all dependencies
   pip install -e ".[all]"

**Permission errors (Windows):**

.. code-block:: bash

   # Run as administrator or use WSL2

**Memory issues:**

.. code-block:: bash

   # Increase available memory
   export CPPFLAGS="-O1"  # Reduce optimization for compilation

Platform-Specific Notes
-----------------------

**macOS with Apple Silicon:**

.. code-block:: bash

   # Use conda-forge for better ARM64 support
   conda install -c conda-forge numpy scipy

**Windows with Visual Studio:**

.. code-block:: bash

   # Ensure Visual Studio Build Tools are installed
   # Set environment variables if needed
   set DISTUTILS_USE_SDK=1

Getting Help
============

* **Documentation**: This guide and the :doc:`user_guide`
* **Issues**: `GitHub Issues <https://github.com/nudca-code/NuDCA/issues>`_
* **Discussions**: `GitHub Discussions <https://github.com/nudca-code/NuDCA/discussions>`_
* **Email**: chohonche@163.com

For contributing guidelines, see :doc:`development`. 
User Guide
==========

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin: 2em 0;">
       <h2 style="color: white; margin-top: 0;">📚 NuDCA User Guide</h2>
       <p style="font-size: 1.1em; margin-bottom: 0;">In-depth information on nuclear decay modeling concepts, astrophysical applications, and computational methods</p>
   </div>

This guide provides comprehensive explanations of NuDCA's key concepts and methodologies. It assumes you have a basic understanding of nuclear physics and are familiar with Python programming. For a quick introduction, see the :doc:`quickstart`.

**Contents:**

This comprehensive user guide covers all aspects of nuclear decay calculations and astrophysical modeling with NuDCA. The content is organized into logical sections that build upon each other, starting with fundamental concepts and progressing to advanced applications.

Core Concepts
=============

Nuclear Decay Fundamentals
---------------------------

**Radioactive Decay Law**

The fundamental equation governing radioactive decay is:

.. math::

   \frac{dN}{dt} = -\lambda N

where :math:`N` is the number of radioactive nuclei, :math:`\lambda` is the decay constant, and :math:`t` is time.

For decay chains with multiple isotopes, this extends to a system of coupled differential equations:

.. math::

   \frac{d\mathbf{N}}{dt} = \mathbf{A} \mathbf{N}

where :math:`\mathbf{N}` is the vector of nuclide abundances and :math:`\mathbf{A}` is the decay matrix.

**Half-Life and Decay Constants**

The relationship between half-life (:math:`t_{1/2}`) and decay constant (:math:`\lambda`) is:

.. math::

   t_{1/2} = \frac{\ln(2)}{\lambda}

NuDCA handles half-life data from various nuclear databases and automatically converts to decay constants for calculations.

**Branching Ratios**

When a parent nucleus can decay through multiple channels, each channel has a probability (branching ratio) :math:`b_i` where:

.. math::

   \sum_i b_i = 1

The partial decay constant for channel :math:`i` is:

.. math::

   \lambda_i = b_i \lambda_{total}

Matrix-Based Approach
---------------------

**Decay Matrix Construction**

NuDCA uses sparse matrix methods to efficiently handle large decay networks. The decay matrix :math:`\mathbf{A}` has the structure:

* Diagonal elements: :math:`A_{ii} = -\lambda_i` (decay out of nuclide :math:`i`)
* Off-diagonal elements: :math:`A_{ij} = \lambda_j b_{j \to i}` (production of :math:`i` from :math:`j`)

**Matrix Exponential Solution**

The solution to the linear system is given by the matrix exponential:

.. math::

   \mathbf{N}(t) = \exp(\mathbf{A}t) \mathbf{N}(0)

NuDCA implements efficient algorithms for computing matrix exponentials of sparse matrices.

**Numerical Stability**

For systems with widely varying time scales (stiff equations), NuDCA uses:

* Adaptive time stepping
* Krylov subspace methods
* Precomputed matrix decompositions

Energy Deposition and Heating
------------------------------

**Decay Energy Types**

NuDCA tracks multiple types of decay energy:

* **Electromagnetic (EM)**: Gamma rays, X-rays, conversion electrons
* **Light Particles (LP)**: Beta particles, positrons, neutrinos
* **Heavy Particles (HP)**: Alpha particles, recoil nuclei
* **Neutrinos**: Usually escape without depositing energy locally

**Thermalization Efficiency**

The fraction of decay energy that thermalizes depends on:

* Energy of the radiation
* Density and composition of the medium
* Time after the explosion

NuDCA provides models for thermalization efficiency based on:

.. math::

   \epsilon(t) = \frac{1}{1 + (t/t_0)^{-n}}

where :math:`t_0` and :math:`n` are material-dependent parameters.

Nuclear Data Management
=======================

**ENDF-B-VIII.1 Integration**

NuDCA integrates nuclear data from the Evaluated Nuclear Data File (ENDF-B-VIII.1), providing:

* Half-lives for ~3000 nuclides
* Decay modes and branching ratios
* Decay energies by particle type
* Fission fragment yields

**Data Validation and Quality Control**

The package includes tools for:

* Cross-checking data consistency
* Identifying missing decay paths
* Validating energy conservation
* Comparing with experimental measurements

**Custom Nuclear Data**

Users can incorporate custom nuclear data by:

.. code-block:: python

   from nudca.utils import DecayDatabaseManager
   
   # Load custom data
   manager = DecayDatabaseManager('custom_data')
   custom_db = manager.load_custom_database('my_data.json')

Astrophysical Applications
===========================

r-Process Nucleosynthesis
--------------------------

**Rapid Neutron Capture**

The r-process occurs in extremely neutron-rich environments where neutron capture rates exceed beta-decay rates. NuDCA models the subsequent decay of neutron-rich isotopes back to stability.

**Key Physics:**

* Initial abundance distribution from r-process simulations
* Beta-decay waiting points
* Fission cycling in superheavy elements
* Final stable abundance patterns

**Typical Workflow:**

.. code-block:: python

   # Load r-process abundance pattern
   inputer = nudca.Inputer(decay_database)
   r_abundances = inputer.read_abundance_from_file('r_process_yields.csv')
   
   # Calculate post-r-process decay
   decay_calc = nudca.RadioactiveDecay(
       initial_abundance=r_abundances,
       decay_database=decay_database,
       decay_matrix=decay_matrix
   )
   
   # Evolution over astrophysical timescales
   times = np.logspace(0, 17, 100)  # 1 s to age of universe
   final_abundances = decay_calc.decay_process(times)

Kilonova Modeling
-----------------

**Energy Deposition**

Kilonovae are powered by radioactive decay of r-process elements. The heating rate per unit mass is:

.. math::

   \dot{q}(t) = \sum_i \lambda_i N_i(t) Q_i \epsilon_i(t)

where :math:`Q_i` is the decay energy and :math:`\epsilon_i(t)` is the thermalization efficiency.

**Light Curve Calculation**

The bolometric luminosity is related to the heating rate by:

.. math::

   L_{bol}(t) = \int \rho(r) \dot{q}(r,t) 4\pi r^2 dr

**Multi-Zone Models**

For realistic kilonova modeling, NuDCA supports:

* Layered ejecta with different compositions
* Velocity-dependent abundance distributions
* Time-dependent opacity calculations

.. code-block:: python

   from nudca.kilonovae import KNeLightCurve
   
   # Create multi-zone kilonova model
   kn_model = KNeLightCurve(
       mass_ejecta=0.05,      # Solar masses
       vel_ejecta=0.2,        # Fraction of c
       opacity_model='Tanaka2020'
   )
   
   # Calculate light curves in multiple bands
   for band in ['u', 'g', 'r', 'i', 'z']:
       times, magnitudes = kn_model.light_curve(
           times=np.linspace(0.1, 30, 100),  # days
           band=band
       )

Best Practices
==============

Performance Optimization
-------------------------

**Memory Management**

For large decay networks:

.. code-block:: python

   # Use sparse matrices
   import scipy.sparse as sp
   
   # Precompute matrix factorizations
   matrix_builder = nudca.MatrixBuilder(decay_database)
   P, P_inv = matrix_builder.build_decay_matrix()
   
   # Cache results for repeated calculations
   decay_cache = {}

**Computational Efficiency**

* Use vectorized operations for time arrays
* Leverage NumPy broadcasting for abundance calculations
* Consider parallel processing for parameter studies

**Numerical Accuracy**

* Choose appropriate time stepping based on shortest half-life
* Monitor conservation laws (mass, energy)
* Use double precision for long-time evolution

Data Analysis Workflows
-----------------------

**Abundance Evolution Studies**

.. code-block:: python

   # Set up parameter study
   initial_compositions = [
       {'lanthanides': 0.1, 'actinides': 0.01},
       {'lanthanides': 0.3, 'actinides': 0.03},
       {'lanthanides': 0.5, 'actinides': 0.05}
   ]
   
   results = []
   for composition in initial_compositions:
       # Run decay calculation
       calc = nudca.RadioactiveDecay(composition, db, matrix)
       result = calc.decay_process(times)
       results.append(result)
   
   # Analyze systematic trends
   nudca.plot.abundance_evolution(times, results)

**Uncertainty Quantification**

* Monte Carlo sampling of nuclear data uncertainties
* Sensitivity analysis for key parameters
* Error propagation through decay calculations

Advanced Topics
===============

Custom Physics Models
----------------------

**Modified Decay Rates**

For extreme environments, decay rates may be modified by:

* High magnetic fields (electron capture rates)
* High temperatures (bound-state beta decay)
* Dense plasmas (screening effects)

**Fission Fragment Distributions**

For superheavy elements, custom fission models can be implemented:

.. code-block:: python

   def custom_fission_yields(Z, A, excitation_energy):
       """Custom fission fragment yield model."""
       # Implementation of specific fission model
       return fragment_yields

Integration with Other Codes
-----------------------------

**Hydrodynamics Coupling**

NuDCA can be integrated with hydrodynamics codes for self-consistent kilonova modeling:

.. code-block:: python

   # Interface with hydro code
   for timestep in hydro_evolution:
       # Update density and temperature
       density = hydro_data.density[timestep]
       temperature = hydro_data.temperature[timestep]
       
       # Calculate heating rates
       heating = decay_calc.heating_rates(time, density, temperature)
       
       # Pass back to hydro code
       hydro_data.heating_source[timestep] = heating

**Radiative Transfer**

For photon transport calculations:

.. code-block:: python

   from nudca.kilonovae import KNOpacity
   
   # Calculate wavelength-dependent opacities
   opacity_calc = KNOpacity()
   kappa_lambda = opacity_calc.opacity_spectrum(
       composition=abundances,
       temperature=temperature,
       density=density,
       wavelengths=wavelength_grid
   )

Troubleshooting
===============

Common Issues
-------------

**Convergence Problems**

* Check for extremely short half-lives requiring small time steps
* Verify nuclear data completeness
* Consider using implicit time stepping methods

**Memory Issues**

* Use sparse matrix formats consistently
* Clear intermediate results when possible
* Consider matrix-free methods for very large problems

**Accuracy Concerns**

* Monitor conservation laws throughout calculation
* Compare with analytical solutions for simple cases
* Use higher-order time integration schemes when needed

**Performance Issues**

* Profile code to identify bottlenecks
* Use compiled extensions for inner loops
* Consider approximate methods for preliminary studies

Getting Help
------------

* Check the :doc:`api` for detailed function documentation
* Browse :doc:`examples` for complete working examples
* Search :doc:`tutorials` for step-by-step guides
* Report issues on `GitHub <https://github.com/nudca-code/NuDCA/issues>`_
* Ask questions on `Discussions <https://github.com/nudca-code/NuDCA/discussions>`_ 
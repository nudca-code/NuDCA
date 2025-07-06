API Reference
=============

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin: 2em 0;">
       <h2 style="color: white; margin-top: 0;">🔍 Complete API Reference</h2>
       <p style="font-size: 1.1em; margin-bottom: 0;">Comprehensive documentation for all NuDCA functions, classes, and modules</p>
   </div>

This reference manual details functions, modules, and objects included in NuDCA, describing what they are and what they do. 
For learning how to use NuDCA, see the :doc:`/user/user_guide` and :doc:`quickstart`.

.. raw:: html

   <div style="background: #f8f9fa; padding: 1.5em; border-radius: 8px; margin: 2em 0;">
       <h3 style="color: #667eea; margin-top: 0;">📚 Organization</h3>
       <p style="color: #4a5568; margin-bottom: 1em;">NuDCA is organized into several modules, each serving specific aspects of nuclear decay calculations and astrophysical modeling:</p>
       <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1em;">
           <div style="background: #e8f5e8; padding: 1em; border-radius: 6px; border-left: 4px solid #48bb78;">
               <h4 style="color: #2f855a; margin: 0;">🔬 Core</h4>
               <p style="color: #4a5568; font-size: 0.9em; margin: 0.5em 0 0;">Nuclear decay calculations and databases</p>
           </div>
           <div style="background: #fff5e6; padding: 1em; border-radius: 6px; border-left: 4px solid #ed8936;">
               <h4 style="color: #dd6b20; margin: 0;">🌟 Kilonova</h4>
               <p style="color: #4a5568; font-size: 0.9em; margin: 0.5em 0 0;">Astrophysical modeling and analysis</p>
           </div>
           <div style="background: #e6f3ff; padding: 1em; border-radius: 6px; border-left: 4px solid #4299e1;">
               <h4 style="color: #2b6cb0; margin: 0;">📊 Visualization</h4>
               <p style="color: #4a5568; font-size: 0.9em; margin: 0.5em 0 0;">Plotting and data visualization</p>
           </div>
           <div style="background: #f0f0ff; padding: 1em; border-radius: 6px; border-left: 4px solid #9f7aea;">
               <h4 style="color: #6b46c1; margin: 0;">🛠️ Utilities</h4>
               <p style="color: #4a5568; font-size: 0.9em; margin: 0.5em 0 0;">Helper functions and I/O operations</p>
           </div>
       </div>
   </div>

.. raw:: html

   <div style="background: #fff3cd; padding: 1.5em; border-radius: 8px; margin: 2em 0; border-left: 4px solid #ffc107;">
       <h3 style="color: #856404; margin-top: 0;">💡 Usage Notes</h3>
       <ul style="color: #856404; margin-bottom: 0;">
           <li>All functions assume NumPy arrays for numerical inputs unless otherwise specified</li>
           <li>Nuclear data follows ENDF-B-VIII.1 conventions and nomenclature</li>
           <li>Time units are in seconds throughout the package</li>
           <li>Abundances are typically mass fractions or number densities</li>
           <li>Energy units are in MeV for decay energies and erg/s for heating rates</li>
           <li>Most functions support both single values and arrays for time-dependent calculations</li>
       </ul>
   </div>

Core Classes and Functions
==========================

.. currentmodule:: nudca

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   RadioactiveDecay
   DecayDatabase
   DecayMatrix
   DecayDiagram
   MatrixBuilder
   Nuclide

Top-level functions
-------------------

.. autosummary::
   :toctree: generated/
   :template: function.rst
   
   load_decay_database
   load_decay_matrix
   serialize_decay_database
   serialize_decay_matrix

.. _decay_core:

Decay Core (`nudca.decay_core`)
===============================

.. currentmodule:: nudca.decay_core

The decay_core module provides the fundamental classes for nuclear decay calculations.

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   RadioactiveDecay

.. _decay_network:

Decay Network (`nudca.decay_network`)
=====================================

.. currentmodule:: nudca.decay_network

Database and matrix management for decay networks.

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   DecayDatabase
   DecayMatrix
   MatrixBuilder

.. _kilonovae:

Kilonova Analysis (`nudca.kilonovae`)
=====================================

Tools for modeling kilonova phenomena and r-process nucleosynthesis.

Geometry and Profiles
---------------------

.. currentmodule:: nudca.kilonovae.geometry

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   Geometry
   DensityProfile
   VelocityProfile

Heating Rates
-------------

.. currentmodule:: nudca.kilonovae.heating_rate

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   RadioactiveHeatingRate
   EffectiveHeatingRate
   ThermalizationEfficiency

Light Curves
------------

.. currentmodule:: nudca.kilonovae.lightcurve

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   KNeLightCurve

Opacity
-------

.. currentmodule:: nudca.kilonovae.opacity

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   KNOpacity

.. _io:

Input/Output (`nudca.io`)
=========================

.. currentmodule:: nudca.io

Data import, export, and abundance calculations.

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   Inputer
   Outputer
   AbundanceEstimator

.. _plot:

Visualization (`nudca.plot`)
============================

.. currentmodule:: nudca.plot

Plotting functions for nuclear data visualization.

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   DecayDiagram

.. _nuclide:

Nuclide (`nudca.nuclide`)
=========================

.. currentmodule:: nudca.nuclide

Classes and functions for handling nuclear species.

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   Nuclide

.. autosummary::
   :toctree: generated/
   :template: exception.rst
   
   NuclideStrError

.. _utils:

Utilities (`nudca.utils`)
=========================

.. currentmodule:: nudca.utils

Helper functions and utility classes.

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   DecayDatabaseManager
   HalfLifeColorMap

.. _constants:

Constants (`nudca.constants`)
=============================

.. currentmodule:: nudca.constants

Physical constants and conversion factors used throughout NuDCA.

.. autosummary::
   :toctree: generated/
   :template: attribute.rst
   
   nudca.constants.PI
   nudca.constants.E_MASS
   nudca.constants.N_MASS
   nudca.constants.P_MASS
   nudca.constants.ATOMIC_MASS
   nudca.constants.NA_CGS
   nudca.constants.DAY_CGS
   nudca.constants.YEAR_CGS
   nudca.constants.EV_CGS
   nudca.constants.G_CGS
   nudca.constants.C_CGS
   nudca.constants.M_SUN_CGS
   nudca.constants.H_CGS
   nudca.constants.K_B_CGS
   nudca.constants.SIGMA_SB_CGS

Version Information
===================

.. autodata:: nudca.__version__
   :no-value:

.. raw:: html

   <div style="background: #d4edda; border: 1px solid #c3e6cb; padding: 1em; border-radius: 8px; margin: 2em 0;">
       <h3 style="margin-top: 0; color: #155724;">📦 Package Information</h3>
       <p>NuDCA follows semantic versioning. Check the <a href="changelog.html">changelog</a> for version history and updates.</p>
   </div> 
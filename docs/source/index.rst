.. NuDCA documentation master file, created by
   sphinx-quickstart on Sat Jul  5 14:39:28 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _mainpage:


.. toctree::
   :maxdepth: 1
   :hidden:

   Installation <installation>
   User Guide <user/index>
   Examples <examples/index>
   Tutorials <tutorials/index>
   API <api/index>
   Physics Background <physics/index>
   Changelog <changelog>

.. raw:: html

   <style>
   .hero-section {
       text-align: center;
       padding: 4em 2em;
       background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
       color: white;
       border-radius: 12px;
       margin-bottom: 3em;
       box-shadow: 0 8px 32px rgba(0,0,0,0.1);
   }
   .hero-section h1 {
       font-size: 3.5em;
       margin-bottom: 0.3em;
       color: white;
       font-weight: 700;
   }
   .hero-section .subtitle {
       font-size: 1.4em;
       color: rgba(255,255,255,0.9);
       margin-bottom: 0.5em;
   }
   .hero-section .description {
       font-size: 1.1em;
       color: rgba(255,255,255,0.8);
       max-width: 600px;
       margin: 0 auto 2em;
   }
   .hero-section .version-badge {
       background: rgba(255,255,255,0.2);
       color: white;
       padding: 0.5em 1em;
       border-radius: 20px;
       font-size: 0.9em;
       margin-bottom: 1em;
   }
   .section-grid {
       display: grid;
       grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
       gap: 2em;
       margin: 3em 0;
   }
   .section-card {
       background: white;
       padding: 2em;
       border-radius: 12px;
       border: 1px solid #e2e8f0;
       transition: all 0.3s ease;
       box-shadow: 0 2px 8px rgba(0,0,0,0.05);
       text-align: center;
   }
   .section-card:hover {
       transform: translateY(-4px);
       box-shadow: 0 8px 25px rgba(0,0,0,0.1);
   }
   .section-card h3 {
       color: #667eea;
       margin-top: 0;
       font-size: 1.4em;
   }
   .section-card p {
       color: #4a5568;
       line-height: 1.6;
       margin-bottom: 1.5em;
   }
   .section-card .btn {
       background: #667eea;
       color: white;
       padding: 0.8em 1.5em;
       border-radius: 6px;
       text-decoration: none;
       font-weight: 500;
       transition: all 0.2s ease;
       display: inline-block;
   }
   .section-card .btn:hover {
       background: #764ba2;
       color: white;
       text-decoration: none;
       transform: translateY(-1px);
   }
   .useful-links {
       background: #f8f9fa;
       padding: 2em;
       border-radius: 8px;
       margin: 3em 0;
   }
   .useful-links h3 {
       color: #667eea;
       margin-top: 0;
   }
   .links-grid {
       display: grid;
       grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
       gap: 1em;
       margin-top: 1em;
   }
   .links-grid a {
       color: #667eea;
       text-decoration: none;
       font-weight: 500;
   }
   .links-grid a:hover {
       color: #764ba2;
       text-decoration: underline;
   }
   </style>

   <div class="hero-section">
       <div class="version-badge">Version: 0.1.0</div>
       <h1>NuDCA documentation</h1>
       <p class="subtitle">Nuclear Decay Chains in Astrophysics</p>
       <p class="description">NuDCA is the fundamental package for nuclear decay calculations and kilonova modeling in Python. It provides comprehensive tools for radioactive decay chain analysis, r-process nucleosynthesis simulations, and astrophysical phenomena modeling.</p>
   </div>

.. raw:: html

   <div class="useful-links">
       <h3>🔗 Useful links</h3>
       <div class="links-grid">
           <a href="installation.html">Installation</a>
           <a href="https://github.com/nudca-code/NuDCA">Source Repository</a>
           <a href="https://github.com/nudca-code/NuDCA/issues">Issue Tracker</a>
           <a href="https://github.com/nudca-code/NuDCA/discussions">Q&A Support</a>
           <a href="mailto:chohonchen@163.com">Mailing List</a>
           <a href="changelog.html">Release Notes</a>
       </div>
   </div>

.. raw:: html

   <div class="section-grid">
       <div class="section-card">
           <h3>🚀 Getting started</h3>
           <p>New to NuDCA? Check out the Quick Start Guide. It contains an introduction to NuDCA's main concepts and links to additional tutorials.</p>
           <a href="quickstart.html" class="btn">To the quick start guide</a>
       </div>
       
       <div class="section-card">
           <h3>📚 User guide</h3>
           <p>The user guide provides in-depth information on the key concepts of nuclear decay modeling with useful background information and explanation.</p>
           <a href="user/index.html" class="btn">To the user guide</a>
       </div>
       
       <div class="section-card">
           <h3>📖 API reference</h3>
           <p>The reference guide contains a detailed description of the functions, modules, and objects included in NuDCA. It assumes you understand the key concepts.</p>
           <a href="api/index.html" class="btn">To the reference guide</a>
       </div>
       
       <div class="section-card">
           <h3>🤝 Contributor's guide</h3>
           <p>Want to add to the codebase? Can help add translation or improve documentation? The contributing guidelines will guide you through the process.</p>
           <a href="development.html" class="btn">To the contributor's guide</a>
       </div>
   </div>


Scientific Applications
=======================

NuDCA enables cutting-edge research in nuclear astrophysics and related fields:

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 2em; margin: 2em 0;">
       <div style="background: #f0f8ff; padding: 1.5em; border-radius: 8px; border-left: 4px solid #4299e1;">
           <h4 style="color: #2b6cb0; margin-top: 0;">🌟 Kilonova AT2017gfo</h4>
           <p style="color: #2c5282;">NuDCA models the nuclear decay chains powering kilonova events, providing insights into r-process element synthesis and energy deposition rates in neutron star mergers.</p>
       </div>
       
       <div style="background: #f0fff4; padding: 1.5em; border-radius: 8px; border-left: 4px solid #48bb78;">
           <h4 style="color: #22543d; margin-top: 0;">⚛️ Nuclear Data Validation</h4>
           <p style="color: #2f855a;">Researchers use NuDCA to validate nuclear decay data from ENDF-B-VIII.1, ensuring consistency in half-lives, branching ratios, and decay energies across databases.</p>
       </div>
       
       <div style="background: #fef5e7; padding: 1.5em; border-radius: 8px; border-left: 4px solid #ed8936;">
           <h4 style="color: #c05621; margin-top: 0;">🔬 r-Process Chains</h4>
           <p style="color: #dd6b20;">Detailed studies of rapid neutron capture processes in extreme astrophysical environments, helping understand the origin of heavy elements in the universe.</p>
       </div>
       
       <div style="background: #faf5ff; padding: 1.5em; border-radius: 8px; border-left: 4px solid #9f7aea;">
           <h4 style="color: #553c9a; margin-top: 0;">📊 Kilonova Modeling</h4>
           <p style="color: #6b46c1;">Modeling the emission of light from the decay of r-process elements in kilonovae, providing insights into the properties of these events.</p>
       </div>
   </div>

Package Overview
================

NuDCA provides comprehensive tools for nuclear decay calculations and astrophysical modeling:

**Core Features:**

* **Nuclear Decay Database**: Complete nuclear data management with ENDF-B-VIII.1 support
* **Matrix-Based Calculations**: Efficient sparse matrix operations for complex decay networks  
* **Kilonova Modeling**: Specialized tools for r-process nucleosynthesis and light curve analysis
* **Visualization**: Advanced plotting capabilities for nuclear charts and decay chains
* **High Performance**: Optimized algorithms leveraging NumPy and SciPy

**Scientific Domains:**

* Nuclear Physics and Radioactive Decay
* Astrophysics and Neutron Star Mergers
* Kilonova and Transient Phenomena
* Nuclear Data Analysis and Validation

Citation
========

If you use NuDCA in your research, please cite:

.. code-block:: bibtex

   @software{nudca,
     author = {XXX},
     title = {NuDCA: A Numerical Code for Nuclear Decay Chains in Astrophysics},
     year = {2025},
     publisher = {GitHub},
     url = {https://github.com/nudca-code/NuDCA}
   }



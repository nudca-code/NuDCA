Quick Start
===========

Get up and running with NuDCA in minutes!

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin: 2em 0; text-align: center;">
       <h2 style="color: white; margin-top: 0; font-size: 2.5em;">🚀 Ready in 3 Steps</h2>
       <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1.5em; margin: 2em 0;">
           <div style="background: rgba(255,255,255,0.1); padding: 1.5em; border-radius: 8px;">
               <h3 style="color: white; margin: 0 0 0.5em 0; font-size: 1.5em;">1️⃣</h3>
               <p style="color: rgba(255,255,255,0.9); margin: 0;">Install NuDCA</p>
           </div>
           <div style="background: rgba(255,255,255,0.1); padding: 1.5em; border-radius: 8px;">
               <h3 style="color: white; margin: 0 0 0.5em 0; font-size: 1.5em;">2️⃣</h3>
               <p style="color: rgba(255,255,255,0.9); margin: 0;">Load nuclear data</p>
           </div>
           <div style="background: rgba(255,255,255,0.1); padding: 1.5em; border-radius: 8px;">
               <h3 style="color: white; margin: 0 0 0.5em 0; font-size: 1.5em;">3️⃣</h3>
               <p style="color: rgba(255,255,255,0.9); margin: 0;">Start calculating!</p>
           </div>
       </div>
   </div>

Installation
------------

.. code-block:: bash

   pip install nudca

Basic Example
-------------

Here's a complete working example to get you started:

.. code-block:: python

   import nudca
   import numpy as np
   import matplotlib.pyplot as plt

   # Load nuclear decay data
   decay_db = nudca.load_decay_database()
   decay_matrix = nudca.load_decay_matrix()

   # Define initial composition (e.g., pure U-238)
   initial = {'U238': 1.0}

   # Create decay calculator
   calculator = nudca.RadioactiveDecay(
       initial_abundance=initial,
       decay_database=decay_db,
       decay_matrix=decay_matrix
   )

   # Calculate decay over time
   times = np.logspace(0, 17, 100)  # 1 second to ~3 billion years
   abundances = calculator.decay_process(times)

   # Plot results
   plt.figure(figsize=(10, 6))
   for nuclide, abundance in abundances.items():
       if max(abundance) > 0.01:  # Only plot significant abundances
           plt.loglog(times, abundance, label=nuclide)
   
   plt.xlabel('Time (s)')
   plt.ylabel('Abundance')
   plt.title('Nuclear Decay Chain of U-238')
   plt.legend()
   plt.grid(True, alpha=0.3)
   plt.show()

Key Concepts
------------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1.5em; margin: 2em 0;">
       <div style="background: #f8f9fa; padding: 1.5em; border-radius: 8px; border-left: 4px solid #667eea;">
           <h3 style="margin-top: 0; color: #667eea;">📊 Decay Database</h3>
           <p>Contains nuclear data including half-lives, decay modes, branching ratios, and decay energies from ENDF-B-VIII.1.</p>
           <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px;">
   db = nudca.load_decay_database()
   t_half = db.get_nuclide_half_life('Co60')
   print(f"Co-60 half-life: {t_half}")</pre>
       </div>
       
       <div style="background: #f8f9fa; padding: 1.5em; border-radius: 8px; border-left: 4px solid #764ba2;">
           <h3 style="margin-top: 0; color: #764ba2;">🔗 Decay Matrix</h3>
           <p>Efficient sparse matrix representation of decay networks for fast computation of complex decay chains.</p>
           <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px;">
   matrix = nudca.load_decay_matrix()
   # Pre-computed for optimal performance</pre>
       </div>
       
       <div style="background: #f8f9fa; padding: 1.5em; border-radius: 8px; border-left: 4px solid #f093fb;">
           <h3 style="margin-top: 0; color: #f093fb;">⚡ Heating Rates</h3>
           <p>Calculate energy deposition from radioactive decay, crucial for kilonova modeling.</p>
           <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px;">
   heating = calculator.decay_heating_rates(
       times, energy_type='EM'
   )</pre>
       </div>
   </div>

Common Use Cases
----------------

**1. Single Nuclide Decay**

.. code-block:: python

   # Decay of a single radioactive isotope
   calc = nudca.RadioactiveDecay(
       {'Co60': 1e6},  # 1 million atoms
       decay_db, decay_matrix
   )
   times = np.linspace(0, 20*365*24*3600, 100)  # 20 years
   result = calc.decay_process(times)

**2. Decay Chain Analysis**

.. code-block:: python

   # Visualize decay chains
   diagram = nudca.DecayDiagram(decay_db)
   diagram.plot_decay_chains('U238', max_generations=10)

**3. Kilonova Light Curves**

.. code-block:: python

   from nudca.kilonovae import KNeLightCurve
   
   # Create kilonova model
   kn = KNeLightCurve(
       mass_ejecta=0.05,  # Solar masses
       vel_ejecta=0.2     # Fraction of c
   )
   
   # Calculate light curve
   times = np.logspace(-1, 2, 100)  # 0.1 to 100 days
   t, mag = kn(times, band='g')

Tips & Tricks
-------------

.. raw:: html

   <div style="background: #e6f3ff; padding: 1.5em; border-radius: 8px; margin: 2em 0;">
       <h3 style="margin-top: 0; color: #0066cc;">💡 Pro Tips</h3>
       <ul>
           <li><strong>Performance:</strong> Use vectorized operations with NumPy arrays for time points</li>
           <li><strong>Memory:</strong> For large decay networks, work with sparse matrices</li>
           <li><strong>Visualization:</strong> Use log scales for time and abundance plots</li>
           <li><strong>Accuracy:</strong> The matrix method ensures high precision even for stiff equations</li>
       </ul>
   </div>

Next Steps
----------

Ready to dive deeper? Check out:

- :doc:`tutorials` - Comprehensive tutorials on specific topics
- :doc:`examples` - Jupyter notebooks with real-world examples  
- :doc:`api` - Complete API reference
- `GitHub Examples <https://github.com/nudca-code/NuDCA/tree/main/examples>`_ - More code examples

Try NuDCA Interactive Examples
================================

.. raw:: html

   <div style="background: #f8f9fa; padding: 2em; border-radius: 8px; margin: 2em 0;">
       <h3 style="color: #667eea; margin-top: 0;">🧪 Interactive Nuclear Physics</h3>
       <p style="color: #4a5568; font-size: 1.1em;">
           Explore nuclear decay chains and kilonova phenomena with these interactive examples.
           Each example demonstrates a core concept in nuclear astrophysics.
       </p>
   </div>

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 2em; margin: 2em 0;">
       <div style="background: white; padding: 2em; border-radius: 8px; border: 1px solid #e2e8f0; box-shadow: 0 2px 8px rgba(0,0,0,0.05);">
           <h4 style="color: #667eea; margin-top: 0;">🔬 Single Isotope Decay</h4>
           <p style="color: #4a5568; margin-bottom: 1em;">Watch how a single radioactive isotope decays over time following exponential decay law.</p>
           <div style="background: #f6f8fa; padding: 1em; border-radius: 4px; font-family: monospace; font-size: 0.9em; border: 1px solid #e2e8f0;">
               <code style="color: #24292e;">
               import nudca<br>
               import numpy as np<br><br>
               # Create Co-60 decay simulation<br>
               calc = nudca.RadioactiveDecay(<br>
               &nbsp;&nbsp;&nbsp;&nbsp;{'Co60': 1e6},  # 1 million atoms<br>
               &nbsp;&nbsp;&nbsp;&nbsp;nudca.load_decay_database(),<br>
               &nbsp;&nbsp;&nbsp;&nbsp;nudca.load_decay_matrix()<br>
               )<br><br>
               # 20 years of decay<br>
               times = np.linspace(0, 20*365*24*3600, 100)<br>
               result = calc.decay_process(times)
               </code>
           </div>
       </div>
       
       <div style="background: white; padding: 2em; border-radius: 8px; border: 1px solid #e2e8f0; box-shadow: 0 2px 8px rgba(0,0,0,0.05);">
           <h4 style="color: #764ba2; margin-top: 0;">⚛️ Decay Chain Analysis</h4>
           <p style="color: #4a5568; margin-bottom: 1em;">Explore complex decay chains like the uranium series with branching ratios and multiple pathways.</p>
           <div style="background: #f6f8fa; padding: 1em; border-radius: 4px; font-family: monospace; font-size: 0.9em; border: 1px solid #e2e8f0;">
               <code style="color: #24292e;">
               import nudca<br><br>
               # Analyze U-238 decay chain<br>
               db = nudca.load_decay_database()<br>
               diagram = nudca.DecayDiagram(db)<br><br>
               # Visualize decay pathways<br>
               diagram.plot_decay_chains(<br>
               &nbsp;&nbsp;&nbsp;&nbsp;'U238', <br>
               &nbsp;&nbsp;&nbsp;&nbsp;max_generations=10<br>
               )
               </code>
           </div>
       </div>
       
       <div style="background: white; padding: 2em; border-radius: 8px; border: 1px solid #e2e8f0; box-shadow: 0 2px 8px rgba(0,0,0,0.05);">
           <h4 style="color: #f093fb; margin-top: 0;">🌟 Kilonova Modeling</h4>
           <p style="color: #4a5568; margin-bottom: 1em;">Simulate the light curves of kilonova events powered by r-process nucleosynthesis.</p>
           <div style="background: #f6f8fa; padding: 1em; border-radius: 4px; font-family: monospace; font-size: 0.9em; border: 1px solid #e2e8f0;">
               <code style="color: #24292e;">
               from nudca.kilonovae import KNeLightCurve<br><br>
               # Create kilonova model<br>
               kn = KNeLightCurve(<br>
               &nbsp;&nbsp;&nbsp;&nbsp;mass_ejecta=0.05,  # Solar masses<br>
               &nbsp;&nbsp;&nbsp;&nbsp;vel_ejecta=0.2     # Fraction of c<br>
               )<br><br>
               # Calculate light curve<br>
               times = np.logspace(-1, 2, 100)<br>
               t, mag = kn(times, band='g')
               </code>
           </div>
       </div>
   </div>

Performance Tips
----------------

.. raw:: html

   <div style="background: #e6f3ff; padding: 1.5em; border-radius: 8px; margin: 2em 0; border-left: 4px solid #4299e1;">
       <h3 style="margin-top: 0; color: #2b6cb0;">⚡ Optimization Strategies</h3>
       <ul style="color: #2c5282; margin-bottom: 0;">
           <li><strong>Vectorization:</strong> Use NumPy arrays for time points to leverage vectorized operations</li>
           <li><strong>Sparse Matrices:</strong> For large decay networks, the sparse matrix format saves memory and computation time</li>
           <li><strong>Precomputed Data:</strong> Load decay databases and matrices once at the beginning of your session</li>
           <li><strong>Time Selection:</strong> Use logarithmic time spacing for decay calculations over many orders of magnitude</li>
           <li><strong>Filtering:</strong> Focus on abundant isotopes (>1% of total) for visualization to reduce clutter</li>
       </ul>
   </div>

Learning Path
-------------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1.5em; margin: 2em 0;">
       <div style="background: #f0fff4; padding: 1.5em; border-radius: 8px; border: 1px solid #9ae6b4;">
           <h4 style="color: #22543d; margin-top: 0;">🎯 Beginner</h4>
           <ul style="color: #2f855a; margin-bottom: 0;">
               <li>Single isotope decay</li>
               <li>Half-life calculations</li>
               <li>Basic plotting</li>
               <li>Nuclear chart exploration</li>
           </ul>
       </div>
       
       <div style="background: #fef5e7; padding: 1.5em; border-radius: 8px; border: 1px solid #fbb6ce;">
           <h4 style="color: #c53030; margin-top: 0;">🔥 Intermediate</h4>
           <ul style="color: #e53e3e; margin-bottom: 0;">
               <li>Complex decay chains</li>
               <li>Branching ratios</li>
               <li>Energy deposition</li>
               <li>Abundance evolution</li>
           </ul>
       </div>
       
       <div style="background: #f0f4ff; padding: 1.5em; border-radius: 8px; border: 1px solid #90cdf4;">
           <h4 style="color: #2c5282; margin-top: 0;">🚀 Advanced</h4>
           <ul style="color: #3182ce; margin-bottom: 0;">
               <li>Kilonova modeling</li>
               <li>r-process nucleosynthesis</li>
               <li>Custom opacity calculations</li>
               <li>Multi-zone simulations</li>
           </ul>
       </div>
   </div>

Community & Support
-------------------

.. raw:: html

   <div style="background: white; padding: 2em; border-radius: 8px; border: 1px solid #e2e8f0; box-shadow: 0 2px 8px rgba(0,0,0,0.05); margin: 2em 0;">
       <h3 style="color: #667eea; margin-top: 0;">🤝 Join the NuDCA Community</h3>
       <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1.5em;">
           <div style="text-align: center;">
               <h4 style="color: #4a5568; margin-bottom: 0.5em;">📧 Email Support</h4>
               <p style="color: #718096; margin-bottom: 0.5em;">Get help with specific problems</p>
               <a href="mailto:chohonche@163.com" style="color: #667eea; text-decoration: none; font-weight: 500;">chohonche@163.com</a>
           </div>
           
           <div style="text-align: center;">
               <h4 style="color: #4a5568; margin-bottom: 0.5em;">💬 GitHub Issues</h4>
               <p style="color: #718096; margin-bottom: 0.5em;">Report bugs or request features</p>
               <a href="https://github.com/nudca-code/NuDCA/issues" style="color: #667eea; text-decoration: none; font-weight: 500;">Open an Issue</a>
           </div>
           
           <div style="text-align: center;">
               <h4 style="color: #4a5568; margin-bottom: 0.5em;">📚 Documentation</h4>
               <p style="color: #718096; margin-bottom: 0.5em;">Complete reference guide</p>
               <a href="api.html" style="color: #667eea; text-decoration: none; font-weight: 500;">API Reference</a>
           </div>
           
           <div style="text-align: center;">
               <h4 style="color: #4a5568; margin-bottom: 0.5em;">🎓 Tutorials</h4>
               <p style="color: #718096; margin-bottom: 0.5em;">Step-by-step guides</p>
               <a href="tutorials.html" style="color: #667eea; text-decoration: none; font-weight: 500;">Learn More</a>
           </div>
       </div>
   </div> 
Examples
========

.. raw:: html

   <div style="text-align: center; margin: 2em 0;">
       <h2 style="color: #667eea;">Learn by Example</h2>
       <p style="font-size: 1.2em; color: #666;">Explore real-world applications of NuDCA through interactive notebooks and code snippets</p>
   </div>

Example Gallery
---------------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1.5em; margin: 2em 0;">
       
       <div style="background: #f8f9fa; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
           <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 1em; color: white;">
               <h3 style="margin: 0;">🌟 AT2017gfo Kilonova</h3>
           </div>
           <div style="padding: 1.5em;">
               <p>Analysis of the famous neutron star merger event GW170817 and its electromagnetic counterpart.</p>
               <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px; font-size: 0.9em;">examples/AT2017gfo.ipynb</pre>
               <ul style="margin: 0.5em 0;">
                   <li>Light curve fitting</li>
                   <li>Ejecta mass estimation</li>
                   <li>r-process nucleosynthesis</li>
               </ul>
           </div>
       </div>
       
       <div style="background: #f8f9fa; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
           <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); padding: 1em; color: white;">
               <h3 style="margin: 0;">⚛️ Decay Networks</h3>
           </div>
           <div style="padding: 1.5em;">
               <p>Build and analyze complex nuclear decay networks with branching ratios.</p>
               <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px; font-size: 0.9em;">examples/decay_network.ipynb</pre>
               <ul style="margin: 0.5em 0;">
                   <li>Matrix construction</li>
                   <li>Chain visualization</li>
                   <li>Equilibrium analysis</li>
               </ul>
           </div>
       </div>
       
       <div style="background: #f8f9fa; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
           <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); padding: 1em; color: white;">
               <h3 style="margin: 0;">☢️ Spontaneous Fission</h3>
           </div>
           <div style="padding: 1.5em;">
               <p>Model spontaneous fission fragments and their subsequent decay chains.</p>
               <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px; font-size: 0.9em;">examples/SF_Fragment.ipynb</pre>
               <ul style="margin: 0.5em 0;">
                   <li>Fragment yields</li>
                   <li>Mass distributions</li>
                   <li>Decay heat calculations</li>
               </ul>
           </div>
       </div>
       
       <div style="background: #f8f9fa; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
           <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); padding: 1em; color: white;">
               <h3 style="margin: 0;">💫 Kilonova Suite</h3>
           </div>
           <div style="padding: 1.5em;">
               <p>Complete kilonova modeling toolkit with multiple specialized notebooks.</p>
               <pre style="background: #e9ecef; padding: 0.5em; border-radius: 4px; font-size: 0.9em;">examples/kilonovae/</pre>
               <ul style="margin: 0.5em 0;">
                   <li><code>geometry.ipynb</code> - Ejecta models</li>
                   <li><code>heating_rate.ipynb</code> - Energy deposition</li>
                   <li><code>lightcurve.ipynb</code> - Photometry</li>
                   <li><code>opacity.ipynb</code> - Radiative transfer</li>
               </ul>
           </div>
       </div>
       
   </div>

Quick Code Examples
-------------------

**Basic Decay Calculation** 🧮

.. code-block:: python

   import nudca
   import numpy as np
   
   # Simple decay of Co-60
   db = nudca.load_decay_database()
   matrix = nudca.load_decay_matrix()
   
   decay = nudca.RadioactiveDecay(
       {'Co60': 1.0},  # Start with 1 atom
       db, matrix
   )
   
   # Calculate over 10 half-lives
   t_half = db.get_nuclide_half_life('Co60', units='s')
   times = np.linspace(0, 10 * t_half, 100)
   
   abundances = decay.decay_process(times)

**Kilonova Light Curve** 🌟

.. code-block:: python

   from nudca.kilonovae import KNeLightCurve
   import matplotlib.pyplot as plt
   
   # Create kilonova model
   kn = KNeLightCurve(
       lightcurve_type='Magnitude',
       mass_ejecta=0.05,     # M_sun
       vel_ejecta=0.2,       # c
       opacity_type='Tanaka'
   )
   
   # Generate light curve
   times = np.logspace(-1, 2, 100)  # 0.1 to 100 days
   t, mag = kn(times, band='g')
   
   # Plot
   plt.figure(figsize=(8, 6))
   plt.plot(t, mag, 'b-', linewidth=2)
   plt.xlabel('Time (days)')
   plt.ylabel('Magnitude')
   plt.gca().invert_yaxis()
   plt.grid(True, alpha=0.3)
   plt.title('Kilonova Light Curve (g-band)')

**Nuclear Chart Visualization** 📊

.. code-block:: python

   from nudca import DecayDiagram
   
   # Create decay diagram plotter
   diagram = DecayDiagram(db)
   
   # Plot section of nuclear chart
   diagram.plot_nuclear_chart(
       Z_range=(80, 100),
       N_range=(120, 150),
       property='half_life',
       log_scale=True
   )

**r-Process Abundance Pattern** 🌊

.. code-block:: python

   from nudca.io import AbundanceEstimator
   import pandas as pd
   
   # Load solar r-process pattern
   solar_r = pd.read_excel('data/solar_r_abundance_pattern.xlsx')
   
   # Estimate initial abundances
   estimator = AbundanceEstimator(db)
   initial = estimator.initial_abundances_rProcess(
       dict(zip(solar_r['Nuclide'], solar_r['Abundance']))
   )
   
   # Calculate decay evolution
   decay = nudca.RadioactiveDecay(initial, db, matrix)
   times = np.logspace(6, 10, 50)  # 10^6 to 10^10 seconds
   final = decay.decay_process(times)

Running the Examples
--------------------

.. raw:: html

   <div style="background: #e6f3ff; padding: 1.5em; border-radius: 8px; margin: 2em 0;">
       <h3 style="margin-top: 0; color: #0066cc;">📝 How to Run</h3>
       <ol>
           <li><strong>Clone the repository:</strong><br>
               <code style="background: white; padding: 0.5em; border-radius: 4px; display: block; margin: 0.5em 0;">git clone https://github.com/nudca-code/NuDCA.git<br>cd NuDCA/examples</code>
           </li>
           <li><strong>Install Jupyter:</strong><br>
               <code style="background: white; padding: 0.5em; border-radius: 4px; display: block; margin: 0.5em 0;">pip install jupyter</code>
           </li>
           <li><strong>Launch notebooks:</strong><br>
               <code style="background: white; padding: 0.5em; border-radius: 4px; display: block; margin: 0.5em 0;">jupyter notebook</code>
           </li>
       </ol>
   </div>

Example Data
------------

The ``data/`` directory contains various datasets:

.. raw:: html

   <table style="width: 100%; margin: 1em 0;">
       <tr style="background: #f8f9fa;">
           <th style="padding: 0.5em; text-align: left;">Dataset</th>
           <th style="padding: 0.5em; text-align: left;">Description</th>
       </tr>
       <tr>
           <td style="padding: 0.5em;"><code>AT2017gfo/</code></td>
           <td style="padding: 0.5em;">Observational data from GW170817 kilonova</td>
       </tr>
       <tr style="background: #f8f9fa;">
           <td style="padding: 0.5em;"><code>Solar_Abundance/</code></td>
           <td style="padding: 0.5em;">Solar system abundance patterns</td>
       </tr>
       <tr>
           <td style="padding: 0.5em;"><code>ENDF-B-VIII.1_*</code></td>
           <td style="padding: 0.5em;">Nuclear decay data from ENDF database</td>
       </tr>
       <tr style="background: #f8f9fa;">
           <td style="padding: 0.5em;"><code>SkyNet_Y_filtered.xlsx</code></td>
           <td style="padding: 0.5em;">Nucleosynthesis yields from SkyNet</td>
       </tr>
   </table>

Contributing Examples
---------------------

We welcome community contributions! To add your own example:

1. Create a well-documented Jupyter notebook
2. Include clear explanations and visualizations
3. Add any required data files
4. Submit a pull request on GitHub

.. raw:: html

   <div style="text-align: center; margin: 3em 0;">
       <a href="https://github.com/nudca-code/NuDCA/tree/main/examples" 
          style="background: #667eea; color: white; padding: 1em 2em; border-radius: 6px; text-decoration: none; display: inline-block;">
          View All Examples on GitHub →
       </a>
   </div> 
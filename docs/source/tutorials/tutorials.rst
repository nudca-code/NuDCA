Tutorials
=========

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin-bottom: 2em; text-align: center;">
       <h2 style="color: white; margin: 0;">Master NuDCA Step by Step</h2>
       <p style="margin: 0.5em 0 0 0;">From basics to advanced kilonova modeling</p>
   </div>

Tutorial 1: Nuclear Decay Fundamentals
--------------------------------------

.. raw:: html

   <div style="display: grid; grid-template-columns: 1fr 2fr; gap: 2em; align-items: start; margin: 2em 0;">
       <div style="background: #f8f9fa; padding: 1.5em; border-radius: 8px;">
           <h3 style="margin-top: 0; color: #667eea;">📚 Key Concepts</h3>
           <ul style="list-style: none; padding: 0;">
               <li>✓ Decay constants: λ = ln(2)/T<sub>1/2</sub></li>
               <li>✓ Bateman equations</li>
               <li>✓ Branching ratios</li>
               <li>✓ Matrix methods</li>
           </ul>
       </div>
       <div>

**Single Isotope Decay**

.. code-block:: python

   import nudca
   import numpy as np
   import matplotlib.pyplot as plt

   # Load data
   db = nudca.load_decay_database()
   matrix = nudca.load_decay_matrix()

   # Decay of Co-60 (t₁/₂ = 5.27 years)
   calc = nudca.RadioactiveDecay({'Co60': 1e6}, db, matrix)
   
   # Time points over 5 half-lives
   t_half = 5.27 * 365.25 * 24 * 3600  # Convert to seconds
   times = np.linspace(0, 5 * t_half, 100)
   
   # Calculate
   abundances = calc.decay_process(times)
   
   # Plot
   plt.figure(figsize=(8, 5))
   plt.semilogy(times / (365.25*24*3600), abundances['Co60'], 'b-', lw=2)
   plt.xlabel('Time (years)')
   plt.ylabel('Number of atoms')
   plt.title('Exponential Decay of Co-60')
   plt.grid(True, alpha=0.3)

.. raw:: html

       </div>
   </div>

**Decay Chain Example**

.. code-block:: python

   # U-238 decay chain
   initial = {'U238': 1.0}
   calc = nudca.RadioactiveDecay(initial, db, matrix)
   
   # Evolution over geological time
   times = np.logspace(6, 17, 200)  # 10^6 to 10^17 seconds
   result = calc.decay_process(times)
   
   # Plot major isotopes
   fig, ax = plt.subplots(figsize=(10, 6))
   for nuclide in ['U238', 'U234', 'Th230', 'Ra226', 'Pb206']:
       if nuclide in result:
           ax.loglog(times/3.15e7, result[nuclide], label=nuclide, lw=2)
   
   ax.set_xlabel('Time (years)')
   ax.set_ylabel('Relative abundance')
   ax.set_title('U-238 Decay Chain Evolution')
   ax.legend()
   ax.grid(True, alpha=0.3)

Tutorial 2: Kilonova Modeling
-----------------------------

.. raw:: html

   <div style="background: #fff3cd; border-left: 4px solid #ffc107; padding: 1em; margin: 1em 0;">
       <strong>🌟 Kilonovae:</strong> Electromagnetic transients powered by r-process radioactive decay following neutron star mergers.
   </div>

**Step 1: Define Ejecta Properties**

.. code-block:: python

   from nudca.kilonovae import KNeLightCurve
   
   # Typical NS merger parameters
   params = {
       'mass_ejecta': 0.05,      # Solar masses
       'vel_ejecta': 0.2,        # Fraction of c
       'opacity_type': 'Tanaka',  # Lanthanide-rich opacity
       'kappa_const': 10.0       # cm²/g
   }
   
   kn = KNeLightCurve(**params)

**Step 2: Calculate Light Curves**

.. code-block:: python

   # Multi-band photometry
   times = np.logspace(-1, 2, 100)  # 0.1 to 100 days
   bands = ['u', 'g', 'r', 'i', 'z', 'J', 'H', 'K']
   
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
   
   # Magnitude evolution
   for band in bands:
       t, mag = kn(times, band=band)
       ax1.plot(t, mag, label=f'{band}-band', lw=2)
   
   ax1.set_xlabel('Time (days)')
   ax1.set_ylabel('Absolute Magnitude')
   ax1.set_title('Kilonova Light Curves')
   ax1.invert_yaxis()
   ax1.legend(loc='lower right')
   ax1.grid(True, alpha=0.3)
   ax1.set_xscale('log')
   
   # Color evolution
   t, g_mag = kn(times, band='g')
   t, r_mag = kn(times, band='r')
   t, i_mag = kn(times, band='i')
   
   ax2.plot(t, g_mag - r_mag, 'b-', label='g-r', lw=2)
   ax2.plot(t, r_mag - i_mag, 'r-', label='r-i', lw=2)
   ax2.set_xlabel('Time (days)')
   ax2.set_ylabel('Color')
   ax2.set_title('Color Evolution')
   ax2.legend()
   ax2.grid(True, alpha=0.3)
   ax2.set_xscale('log')

Tutorial 3: Advanced Analysis
-----------------------------

**Heating Rate Calculations**

.. code-block:: python

   from nudca.kilonovae import RadioactiveHeatingRate
   
   # r-process abundance pattern
   initial_abundance = {
       'U238': 1e-8, 'Th232': 3e-8, 'Eu153': 1e-6,
       'Gd156': 2e-6, 'Pt195': 5e-7, 'Au197': 3e-7
   }
   
   # Calculate heating
   heating = RadioactiveHeatingRate(initial_abundance, db, matrix)
   times = np.logspace(3, 8, 100)  # 10^3 to 10^8 seconds
   
   # Different energy channels
   heat_em = heating.calculate(times, energy_type='EM')
   heat_total = heating.calculate(times, energy_type='Total')
   
   # Thermalization efficiency
   from nudca.kilonovae import ThermalizationEfficiency
   therm = ThermalizationEfficiency()
   f_th = therm(times, mass_ejecta=0.05, vel_ejecta=0.2)
   
   # Effective heating
   heat_eff = heat_total * f_th
   
   # Plot
   plt.figure(figsize=(10, 6))
   plt.loglog(times/86400, heat_total, 'b-', label='Total', lw=2)
   plt.loglog(times/86400, heat_em, 'r--', label='EM only', lw=2)
   plt.loglog(times/86400, heat_eff, 'g:', label='Thermalized', lw=3)
   plt.xlabel('Time (days)')
   plt.ylabel('Heating rate (erg/s/g)')
   plt.title('Radioactive Heating in Kilonovae')
   plt.legend()
   plt.grid(True, alpha=0.3)

**Nuclear Chart Visualization**

.. code-block:: python

   from nudca import DecayDiagram
   
   diagram = DecayDiagram(db)
   
   # Plot r-process path region
   fig, ax = plt.subplots(figsize=(12, 8))
   diagram.plot_nuclear_chart(
       Z_range=(50, 85),
       N_range=(80, 140),
       property='half_life',
       log_scale=True,
       ax=ax
   )
   
   # Overlay r-process path
   # (simplified illustration)
   N_rprocess = np.arange(82, 130)
   Z_rprocess = N_rprocess / 1.5  # N/Z ~ 1.5 for r-process
   ax.plot(N_rprocess, Z_rprocess, 'r-', lw=3, label='r-process path')
   ax.legend()

Best Practices
--------------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1.5em; margin: 2em 0;">
       <div style="background: #e6f3ff; padding: 1.5em; border-radius: 8px;">
           <h3 style="margin-top: 0; color: #0066cc;">⚡ Performance</h3>
           <ul>
               <li>Use vectorized operations</li>
               <li>Leverage sparse matrices</li>
               <li>Pre-compute decay matrices</li>
               <li>Cache frequently used data</li>
           </ul>
       </div>
       <div style="background: #e8f5e9; padding: 1.5em; border-radius: 8px;">
           <h3 style="margin-top: 0; color: #2e7d32;">📊 Visualization</h3>
           <ul>
               <li>Use log scales for time/abundance</li>
               <li>Choose appropriate color maps</li>
               <li>Add error bars when available</li>
               <li>Label axes with units</li>
           </ul>
       </div>
       <div style="background: #fff3e0; padding: 1.5em; border-radius: 8px;">
           <h3 style="margin-top: 0; color: #ef6c00;">🔬 Physics</h3>
           <ul>
               <li>Check unit consistency</li>
               <li>Validate against known results</li>
               <li>Consider numerical precision</li>
               <li>Document assumptions</li>
           </ul>
       </div>
   </div>

Next Steps
----------

Ready for more? Explore:

- 📓 :doc:`Example Notebooks <examples>` - Complete working examples
- 🔧 :doc:`API Reference <api>` - Detailed function documentation
- 🌟 `Research Applications <https://github.com/nudca-code/NuDCA/wiki>`_ - Real-world use cases

.. raw:: html

   <div style="text-align: center; margin: 3em 0;">
       <p style="font-size: 1.2em; color: #666;">
           Questions? Join our <a href="https://github.com/nudca-code/NuDCA/discussions">community discussions</a> 
           or <a href="mailto:chohonche@163.com">contact us</a>
       </p>
   </div> 
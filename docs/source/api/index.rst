
########################
NuDCA API reference
########################

:Release: |version|
:Date: |today|

.. raw:: html

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2em; border-radius: 8px; margin: 2em 0;">
       <h2 style="color: white; margin-top: 0;">🔍 Complete API Reference</h2>
       <p style="font-size: 1.1em; margin-bottom: 0;">Comprehensive documentation for all NuDCA functions, classes, and modules</p>
   </div>

This reference manual details functions, modules, and objects included in NuDCA, describing what they are and what they do. 
For learning how to use NuDCA, see the :doc:`/user/user_guide` and :doc:`/user/quickstart`.

.. raw:: html

   <div style="background: #f8f9fa; padding: 1.5em; border-radius: 8px; margin: 2em 0;">
       <h3 style="color: #667eea; margin-top: 0;">📚 Organization</h3>
       <p style="color: #4a5568; margin-bottom: 1em;">NuDCA is organized into several modules, each serving specific aspects of nuclear decay calculations and astrophysical modeling:</p>
       <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1em;">
           <div style="background: #e8f5e8; padding: 1em; border-radius: 6px; border-left: 4px solid #48bb78;">
               <h4 style="color: #2f855a; margin: 0;">🔬 Radioactive Decay</h4>
               <p style="color: #4a5568; font-size: 0.9em; margin: 0.5em 0 0;">Nuclear decay calculations and databases</p>
           </div>
           <div style="background: #fff5e6; padding: 1em; border-radius: 6px; border-left: 4px solid #ed8936;">
               <h4 style="color: #dd6b20; margin: 0;">🌟 Kilonovae</h4>
               <p style="color: #4a5568; font-size: 0.9em; margin: 0.5em 0 0;">Astrophysical modeling and analysis</p>
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



Python API
==========

.. toctree::
   :maxdepth: 1

   module_structure




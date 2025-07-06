Installation
============

.. raw:: html

   <div style="text-align: center; margin: 2em 0;">
       <div style="display: inline-block; background: #f8f9fa; padding: 1.5em 3em; border-radius: 8px; border: 2px solid #667eea;">
           <h2 style="margin: 0; color: #667eea;">pip install nudca</h2>
           <p style="margin: 0.5em 0 0 0; color: #666;">That's it! You're ready to go.</p>
       </div>
   </div>

Requirements
------------

.. raw:: html

   <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1em; margin: 1em 0;">
       <div>
           <h3>✅ Required</h3>
           <ul>
               <li>Python ≥ 3.8</li>
               <li>NumPy ≥ 1.20.0</li>
               <li>SciPy ≥ 1.7.0</li>
               <li>Pandas ≥ 1.3.0</li>
               <li>Matplotlib ≥ 3.4.0</li>
               <li>NetworkX</li>
               <li>Numba</li>
               <li>Astropy</li>
           </ul>
       </div>
       <div>
           <h3>📦 Optional</h3>
           <ul>
               <li>Jupyter - for notebooks</li>
               <li>h5py - for HDF5 support</li>
           </ul>
       </div>
   </div>

Installation Methods
--------------------

**Standard Installation**

.. code-block:: bash

   pip install nudca

**Development Installation**

.. code-block:: bash

   git clone https://github.com/nudca-code/NuDCA.git
   cd NuDCA
   pip install -e ".[dev]"

**Conda Installation** (coming soon)

.. code-block:: bash

   conda install -c conda-forge nudca

Verify Installation
-------------------

.. code-block:: python

   import nudca
   print(f"NuDCA version: {nudca.__version__}")
   
   # Test basic functionality
   db = nudca.load_decay_database()
   print(f"Loaded {len(db.nuclides)} nuclides")

Troubleshooting
---------------

.. raw:: html

   <details style="margin: 1em 0;">
   <summary style="cursor: pointer; padding: 0.5em; background: #f8f9fa; border-radius: 4px;">
       <strong>🔧 Common Issues (click to expand)</strong>
   </summary>
   <div style="padding: 1em; border-left: 3px solid #667eea; margin-top: 0.5em;">

**ImportError: No module named 'nudca'**

Make sure NuDCA is installed in your current Python environment:

.. code-block:: bash

   pip list | grep nudca

**Version Conflicts**

Create a fresh virtual environment:

.. code-block:: bash

   python -m venv nudca_env
   source nudca_env/bin/activate  # Windows: nudca_env\Scripts\activate
   pip install nudca

**Missing Dependencies**

Install all requirements:

.. code-block:: bash

   pip install numpy scipy pandas matplotlib networkx numba astropy

.. raw:: html

   </div>
   </details>

Platform-Specific Notes
-----------------------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 1em; margin: 2em 0;">
       <div style="text-align: center; padding: 1em; background: #f8f9fa; border-radius: 8px;">
           <h3>🐧 Linux</h3>
           <p>Fully supported. Tested on Ubuntu, Fedora, and CentOS.</p>
       </div>
       <div style="text-align: center; padding: 1em; background: #f8f9fa; border-radius: 8px;">
           <h3>🍎 macOS</h3>
           <p>Fully supported. Both Intel and Apple Silicon.</p>
       </div>
       <div style="text-align: center; padding: 1em; background: #f8f9fa; border-radius: 8px;">
           <h3>🪟 Windows</h3>
           <p>Fully supported. Use Anaconda for easier setup.</p>
       </div>
   </div>

Next Steps
----------

✨ **Installation complete!** Now you can:

- 📖 Follow the :ref:`Quick Start Guide <mainpage>`
- 🎯 Try the :ref:`Examples <examples/index>`
- 🔬 Explore the :ref:`API Reference <api/index>`

Need Help?
----------

- 📋 Check `GitHub Issues <https://github.com/nudca-code/NuDCA/issues>`_
- 📧 Contact: chohonchen@163.com 
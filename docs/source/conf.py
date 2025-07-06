# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# Add the package directory to Python path
package_dir = Path(__file__).parent.parent.parent
sys.path.insert(0, str(package_dir))
sys.path.insert(0, str(package_dir / 'nudca'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'NuDCA'
copyright = '2025, NuDCA Team'
author = 'Qiuhong Chen & Menghua Chen'
language = 'en'
version = '0.1.0'
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx_copybutton',
    'sphinx_design',
    'sphinx_tags',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_book_theme'
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

# Custom CSS files are now defined below in the comprehensive section

# Theme options
html_theme_options = {
    "github_url": "https://github.com/nudca-code/NuDCA",
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "logo": {
        "text": "NuDCA",
        "alt_text": "NuDCA - Nuclear Decay Chains in Astrophysics",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/nudca-code/NuDCA",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/nudca",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
        {
            "name": "Discussions",
            "url": "https://github.com/nudca-code/NuDCA/discussions",
            "icon": "fa-solid fa-comments",
            "type": "fontawesome",
        },
    ],
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links", "theme-switcher"],
    "navbar_persistent": ["search-button"],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink"],
}

# HTML title
html_title = "NuDCA Documentation"

# Favicon
html_favicon = None  # Add a favicon.ico file if available

# HTML context for GitHub integration
html_context = {
    "github_user": "nudca-code",
    "github_repo": "NuDCA", 
    "github_version": "main",
    "doc_path": "docs/source",
}

# -- Extension configuration -------------------------------------------------

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'no-index': False,
}

# Prevent duplicate warnings
add_module_names = False

# Autosummary settings
autosummary_generate = True

# Copy button settings
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = "\\"

# MathJax settings for better math rendering
mathjax3_config = {
    'tex': {
        'packages': {'[+]': ['ams', 'color']},
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
    },
}


# TODO extension settings
todo_include_todos = True
todo_emit_warnings = True

# Doctest settings
doctest_global_setup = '''
import numpy as np
import nudca
'''

# MyST parser settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

# Custom CSS and JS files
html_css_files = [
    'custom.css',
]

# html_js_files = [
#     'custom.js',
# ]

# Additional template paths
html_extra_path = []

# Social cards and metadata
html_meta = {
    'description': 'NuDCA: Nuclear Decay Chains in Astrophysics - A Python package for nuclear decay calculations and kilonova modeling',
    'keywords': 'nuclear physics, astrophysics, decay chains, kilonova, r-process, nucleosynthesis, python',
    'author': 'NuDCA Team',
    'viewport': 'width=device-width, initial-scale=1',
    'theme-color': '#667eea',
}

# OpenGraph metadata
html_meta.update({
    'og:title': 'NuDCA Documentation',
    'og:description': 'Nuclear Decay Chains in Astrophysics - A Python package for nuclear decay calculations and kilonova modeling',
    'og:type': 'website',
    'og:url': 'https://nudca-code.github.io/nudca-main/',
    'og:image': 'https://nudca-code.github.io/nudca-main/_static/nudca-logo.png',
    'twitter:card': 'summary_large_image',
    'twitter:title': 'NuDCA Documentation',
    'twitter:description': 'Nuclear Decay Chains in Astrophysics - A Python package for nuclear decay calculations and kilonova modeling',
})

# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------
project = 'Verkko-fillet'
author = 'Juhyun Kim'
release = '0.1'  # The full version, including alpha/beta/rc tags

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',      # Include documentation from docstrings
    'sphinx.ext.napoleon',     # Support for Google-style and NumPy-style docstrings
    'sphinx.ext.viewcode',     # Add links to highlighted source code
]

#templates_path = ['_templates']  # Path to custom HTML templates
exclude_patterns = []  # List of patterns to ignore when looking for source files

# -- Options for HTML output -------------------------------------------------
html_theme = 'alabaster'  # Choose a theme (default is 'alabaster')
#html_static_path = ['_static']  # Path to custom static files (e.g., CSS)

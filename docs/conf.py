from __future__ import annotations
import warnings
warnings.filterwarnings('ignore')
import sys
import os
from datetime import datetime
from functools import partial
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING
from docutils import nodes


# Ensure Biopython is imported correctly
try:
    from Bio import BiopythonWarning
except ImportError:
    print("Biopython is not installed. Please install it using 'pip install biopython'.")
    sys.exit(1)

warnings.simplefilter('ignore', BiopythonWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="Bio")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # Code that generates warnings
    import verkkofillet

from packaging.version import Version
import plotly.io as pio


pio.renderers.default = "json"

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "extensions"))

# -- General configuration ------------------------------------------------

nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = "4.0"  # Nicer param docs
suppress_warnings = [
    "myst.header",  # https://github.com/executablebooks/MyST-Parser/issues/262
]

# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------
project = 'Verkko-fillet'
author = 'Juhyun Kim'
release = '0.1'  # The full version, including alpha/beta/rc tags
repository_url = "https://github.com/jjuhyunkim/verkko-fillet"
copyright = f"{datetime.now():%Y}, the verkko-fillet development team"

# -- General configuration ---------------------------------------------------
def linkcode_resolve(domain, info):
    if domain != 'py':
        return None
    if not info['module']:
        return None
    filename = info['module'].replace('.', '/')
    return "https://github.com/jjuhyunkim/verkko-fillet/%s.py" % filename


extensions = [
    #"myst_parser",
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "sphinxcontrib.bibtex",
    "matplotlib.sphinxext.plot_directive",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "sphinx.ext.linkcode",
    "sphinx_design",
    "sphinx_tabs.tabs",
]

# source_suffix = {
#     '.rst': 'restructuredtext',
#     '.txt': 'restructuredtext',
#     '.md': 'markdown',
# }

# myst_commonmark_only = False

# Bibliography settings
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

# default settings
templates_path = ["_templates"]
# master_doc = "index"
root_doc = "index"
default_role = "literal"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "release-notes/[!i]*.md",
]

# -- Options for jupyter notebooks ----
nb_execution_mode = "off"  # Prevents notebook execution
nb_output_mime_renderers = {
    "application/vnd.plotly.v1+json": "json"
}

# -- Options for HTML output ----------------------------------------------
html_static_path = ["_static"]
html_show_sphinx = False
html_logo = "_static/verkko-fillet_logo.png"
html_title = "verkko-fillet"
html_theme = "sphinx_book_theme"


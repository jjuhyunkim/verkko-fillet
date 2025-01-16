from __future__ import annotations

import sys
from datetime import datetime
from functools import partial
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

from docutils import nodes
from packaging.version import Version

import os
import sys
sys.path.insert(0, os.path.abspath('.'))
HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]

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
extensions = [
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

def linkcode_resolve(domain, info):
    if domain != 'py':
        return None
    if not info['module']:
        return None
    filename = info['module'].replace('.', '/')
    return "https://github.com/jjuhyunkim/verkko-fillet/%s.py" % filename

# Bibliography settings
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

# default settings
templates_path = ["_templates"]
master_doc = "index"
default_role = "literal"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    # exclude all 0.x.y.md files, but not index.md
    "release-notes/[!i]*.md",
]

# -- Options for jupyter notebooks ----
nb_execution_mode = "off"  # Prevents notebook execution


# -- Options for HTML output ----------------------------------------------

# The theme is sphinx-book-theme, with patches for readthedocs-sphinx-search

html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
}
html_static_path = ["_static"]
html_show_sphinx = False
html_logo = "_static/verkko-fillet_logo.png"
html_title = "verkko-fillet"


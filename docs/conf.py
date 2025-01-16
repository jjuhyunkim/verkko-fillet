from __future__ import annotations

import sys
from datetime import datetime
from functools import partial
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

import matplotlib  # noqa
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
    "git_ref",  # needs to be before scanpydoc.rtd_github_links
    "scanpydoc",  # needs to be before sphinx.ext.linkcode
    "sphinx.ext.linkcode",
    "sphinx_design",
    "sphinx_tabs.tabs",
    "sphinx_search.extension",
    "sphinxext.opengraph",
    *[p.stem for p in (HERE / "extensions").glob("*.py") if p.stem not in {"git_ref"}],
]

#templates_path = ['_templates']  # Path to custom HTML templates
exclude_patterns = []  # List of patterns to ignore when looking for source files

# -- Options for HTML output -------------------------------------------------
# html_theme = 'sphinx_rtd_theme'
html_theme = 'alabaster'  # Choose a theme (default is 'alabaster')
#html_static_path = ['_static']  # Path to custom static files (e.g., CSS)


from __future__ import annotations
"""Post-verkko cleaning and gap filling in Python."""

import sys 

from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl
from ._run_shell import run_shell
from ._default_func import check_user_input,print_directory_tree,addHistory

sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})

import sys

from packaging.version import Version
import importlib.metadata

print("verkko-fillet version:" + importlib.metadata.version("verkkofillet"))


__all__ = [
    "pp",
    "pl",
    "tl",
    'run_shell',
    'check_user_input',
    'print_directory_tree',
    'addHistory',
]
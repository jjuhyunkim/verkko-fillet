"""Post-verkko cleaning and gap filling in Python."""

import sys 

from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl
from ._run_shell import run_shell
from ._default_func import check_user_input,print_directory_tree,addHistory

sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})

__version__ = '0.1.0'
print(f"{__name__} version {__version__}")

__all__ = [
    "pp",
    "pl",
    "tl",
    'run_shell',
    'check_user_input',
    'print_directory_tree',
    'addHistory',
]
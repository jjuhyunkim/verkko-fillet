import warnings
# Ensure Biopython is imported correctly
try:
    from Bio import BiopythonWarning
except ImportError:
    print("Biopython is not installed. Please install it using 'pip install biopython'.")
    sys.exit(1)

warnings.simplefilter('ignore', BiopythonWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="Bio")
warnings.filterwarnings('ignore')

from ._graphAlign import graphIdx,graphAlign
from ._chrAssign import chrAssign,showPairwiseAlign,convertRefName
from ._kmer import mkMeryl,calQV
from ._asm_stats import getT2T
from ._run_rm_rDNA import rmrDNA
from ._insert_gap import insertGap
from ._modiFasta import flipContig,renameContig,sortContig,filterContigs
from ._detect_internal_telomere import detect_internal_telomere,runTrimming

__all__ = [
    'graphIdx',
    'graphAlign',
    'chrAssign',
    'mkMeryl',
    'calQV',
    'getT2T',
    'rmrDNA',
    'insertGap',
    'convertRefName',
    'showPairwiseAlign',
    'flipContig',
    'runTrimming',
    'detect_internal_telomere',
    'renameContig',
    'sortContig',
    'filterContigs',
]

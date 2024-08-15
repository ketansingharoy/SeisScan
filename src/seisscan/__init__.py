"""
A Python package for microearthquake detection and location.
"""


# read version from installed package
from importlib.metadata import version
__version__ = version("seisscan")


from seisscan.read_data.read_fdsn import read_fdsn, read_fdsn_inventory
from seisscan.read_data.read_example import read_example
from seisscan.imaging.waveform import prs
from seisscan.waveform_similarity.peak_cross_correlation import do_pcc
from seisscan.waveform_similarity.local_similarity import do_ls
from seisscan.backprojection.brightness import Brightness4, Stack
from seisscan.backprojection.backprojection import do_bp, prepare_traveltime_lookup_table

__all__ = [
    'read_example',
    'read_fdsn',
    'read_fdsn_inventory',
    'prs',
    'do_pcc',
    'do_ls',
    'prepare_traveltime_lookup_table',
    'do_bp',
    'Brightness4',
    'Stack'
]
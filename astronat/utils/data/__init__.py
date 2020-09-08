# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Data Utilities.

There are a number of x-matching utilities.

"""

__author__ = "Nathaniel Starkman"


__all__ = [
    # modules
    "core",
    "crossmatch",
    # functions
    "xmatch_data_graph",
    "indices_xmatch_fields",
    "xmatch_fields",
    "non_xmatched",
    "indices_xmatch_coords",
    "coord_has_duplicates",
    "xmatch_coords",
    # "indices_xfind_coords",  # TODO
    # "xfind_coords",  # TODO
    "xmatch",
    # conversions
    "xmatch_decorator",
    "Table_to_BaseCoordinateFrame",  # TODO Move to more general location
    "BasecoordinateFrame_to_SkyCoord",  # TODO move to more general location
]


##############################################################################
# IMPORTS

from utilipy.data_utils.crossmatch import (
    indices_xmatch_fields,
    non_xmatched,
    xmatch_fields,
)

from . import core, crossmatch
from .core import BasecoordinateFrame_to_SkyCoord  # conversion functions
from .core import Table_to_BaseCoordinateFrame, xmatch_data_graph
from .crossmatch import indices_xmatch_coords  # other
from .crossmatch import (
    coord_has_duplicates,
    xmatch,
    xmatch_coords,
    xmatch_decorator,
)

##############################################################################
# END

# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Data Utilities.

There are a number of x-matching utilities.

"""

__author__ = "Nathaniel Starkman"


__all__ = [
    # modules
    "core",
    "xmatch",
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
    "XMatch",
    # conversions
    "xmatch_decorator",
    "Table_to_BaseCoordinateFrame",  # TODO Move to more general location
    "BasecoordinateFrame_to_SkyCoord",  # TODO move to more general location
]


##############################################################################
# IMPORTS

# BUILT IN

# THIRD PARTY

from utilipy.data_utils.xmatch import (
    indices_xmatch_fields,
    xmatch_fields,
    non_xmatched,
)


# PROJECT-SPECIFIC

from . import core, xmatch

from .core import (
    xmatch_data_graph,
    # conversion functions
    Table_to_BaseCoordinateFrame,
    BasecoordinateFrame_to_SkyCoord,
)
from .xmatch import (
    indices_xmatch_coords,
    coord_has_duplicates,
    xmatch_coords,
    XMatch,
    # other
    xmatch_decorator,
)


##############################################################################
# END

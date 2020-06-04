# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : typing
# PROJECT : astronat
#
# ----------------------------------------------------------------------------

"""Typing library utilities."""

__author__ = "Nathaniel Starkman"
__credits__ = ["typing"]


__all__ = [
    "array_like",
    "ExpandedUnitType",
]


###############################################################################
# IMPORTS

# BUILT-IN

import typing as T
from typing import *

# CUSTOM

from utilipy.utils.typing import array_like


# THIRD PARTY

from astropy.units import Unit, Quantity


###############################################################################
# CODE
###############################################################################

#####################################################################
# Astropy

ExpandedUnitType = T.TypeVar("ExpandedUnitType", Unit, Quantity)
"""Expanded Unit Type.

Type variable for Astropy [astropy]_ units, including quantities.

References
----------
.. [astropy] Astropy Collaboration et al., 2018, AJ, 156, 123.

Examples
--------

    >>> x: ExpandedUnitType = 10 * u.km
    >>> isinstance(x, u.Unit)
    True

"""


###############################################################################
# __ALL__

if hasattr(T, "__all__"):
    __all__ += T.__all__
else:
    __all__ += list(dir(T))


###############################################################################
# END

# -*- coding: utf-8 -*-

"""Astropy Units. Extended.

provides full drop-in replacement for astropy units.

"""

__author__ = "Nathaniel Starkman"


__all__ = [
    # top-level
    "amuse",
    "composite",
    "convert",
    "core",
    "decorators",
    "full_amuse",
    # core
    "quantity_return_",
    # decorators
    "quantity_output",
    "quantity_io",
    "QuantityInputOutput",
]


#############################################################################
# IMPORTS

# THIRD-PARTY

from astropy.units import *  # for drop-in use, parts will be overridden
from astropy.units import (  # explicit imports for mypy compatibility
    # units
    rad,
    deg,
    m,
    AU,
    pc,
    mag,
    # functions
    get_physical_type,
)


# PROJECT-SPECIFIC

# units to add
# here first so that other stuff can access them
from .amuse import *
from .composite import *

# more stuff
from .convert import from_amuse, hms_str_to_unit
from .core import quantity_return_
from .decorators import quantity_output, quantity_io
from .decorators import QuantityInputOutput

# Import modules into top-level directory
from . import amuse, composite, convert, core, decorators, full_amuse


##############################################################################
# CODE
##############################################################################


##############################################################################
# END

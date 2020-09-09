# -*- coding: utf-8 -*-

"""Astropy Units. Extended.

provides full drop-in replacement for astropy units.

"""

__author__ = "Nathaniel Starkman"


#############################################################################
# IMPORTS

# for drop-in use, parts will be overridden;
# explicit imports for mypy compatibility; units; functions
from astropy import units
from astropy.units import *  # noqa
from astropy.units import AU, deg, get_physical_type, m, mag, pc, rad

# Import modules into top-level directory
from . import amuse, composite, convert, core, decorators, full_amuse

# units to add
# here first so that other stuff can access them
from .amuse import *  # noqa
from .composite import *  # noqa

# more stuff
from .convert import from_amuse, hms_str_to_unit
from .core import quantity_return_
from .decorators import QuantityInputOutput, quantity_io  # , quantity_output

##############################################################################
# PARAMETERS

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
    # convert
    "from_amuse",
    "hms_str_to_unit",
    # decorators
    # "quantity_output",
    "quantity_io",
    "QuantityInputOutput",
]

# __all__ += composite.__all__
# __all__ += amuse.__all__
# __all__ += units.__all__

##############################################################################
# END

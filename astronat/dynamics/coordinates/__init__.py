# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Coordinate Systems for use with Orbits and Streams.

..todo ::

    move some of this to astronat package

"""

__credits__ = ["Astropy", "Jo Bovy", "Jeremy Webb"]


##############################################################################
# IMPORTS

from . import attributes, builtin_frames, ic, representations
from .builtin_frames import *
from .ic import PhaseSpacePosition
from .representations import *

# __all__
__all__ = [
    # modules
    "builtin_frames",
    "attributes",
    "representations",
    "ic",
    # specific
    "PhaseSpacePosition",
]
__all__ += builtin_frames.__all__
__all__ += representations.__all__
# "OrbitSkyOffsetRepresentation",
# "OrbitSkyOffsetUnitRepresentation",

##############################################################################
# END

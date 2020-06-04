# -*- coding: utf-8 -*-

"""Dynamics (Mostly Galactic).

Built for Galpy and Gala inter-operability to permit developing a single code
and testing on two pipelines.

"""

__author__ = "Nathaniel Starkman"
__credits__ = [
    "Jo Bovy & galpy contributors",
    "Adrian Price-Whelan & Gala contributors",
]

__all__ = [
    # modules
    "coordinates",
    "common",
    "potential",
    # functions
    "PhaseSpacePosition",
]


##############################################################################
# IMPORTS

# BUILT IN


# THIRD PARTY


# PROJECT-SPECIFIC

# modules
from . import coordinates, common, potential

# functions
from .coordinates import PhaseSpacePosition


##############################################################################
# END

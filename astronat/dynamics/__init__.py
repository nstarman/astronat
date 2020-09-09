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
    # "vesc",
    # Types
    "PotentialType",
    "GalpyBasePotentialType",
    "GalpyPotentialType",
    "GalaPotentialType",
    # transformations
    "galpy_Orbit_to_SkyCoord",
    "galpy_Orbit_to_gala_Orbit",
    "gala_Orbit_to_galpy_Orbit",
    # functions
    "PhaseSpacePosition",
    "GalpyCompositePotential",
]


##############################################################################
# IMPORTS

# modules
from . import common, coordinates, potential
from .common import GalaPotentialType, PotentialType

# functions
from .coordinates import PhaseSpacePosition
from .potential import (
    GalpyBasePotentialType,
    GalpyCompositePotential,
    GalpyPotentialType,
)

# from . import vesc


##############################################################################
# END

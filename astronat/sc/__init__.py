# -*- coding: utf-8 -*-

"""astro functions where the arguments are SkyCoords."""

__author__ = "Nathaniel Starkman"


__all__ = [
    # modules
    "core",
    # distance modulus
    "apparent_to_absolute_magnitude",
    "absolute_to_apparent_magnitude",
    "distanceModulus_magnitude",
    "distanceModulus_distance",
    "distanceModulus",
    # parallax
    "parallax_angle",
    "parallax_distance",
    "parallax",
    # angular separation
    "max_angular_separation",
]


#############################################################################
# IMPORTS

from .core import (
    # distance modulus
    apparent_to_absolute_magnitude,
    absolute_to_apparent_magnitude,
    distanceModulus_magnitude,
    distanceModulus_distance,
    distanceModulus,
    # parallax
    parallax_angle,
    parallax_distance,
    parallax,
    # angular separation
    max_angular_separation,
)

# module
from . import core


#############################################################################
# END

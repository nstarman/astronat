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

from . import core

from .core import (  # distance modulus; parallax; angular separation
    absolute_to_apparent_magnitude,
    apparent_to_absolute_magnitude,
    distanceModulus,
    distanceModulus_distance,
    distanceModulus_magnitude,
    max_angular_separation,
    parallax,
    parallax_angle,
    parallax_distance,
)

#############################################################################
# END

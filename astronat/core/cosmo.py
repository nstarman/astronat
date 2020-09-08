# -*- coding: utf-8 -*-

"""Cosmology related functions.

These functions are implemented here in minimal form for low overhead.
For equivalent functions with the full set of bells and whistles, see
:mod:`~astronat.common`.

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "scale_factor",
    "hubble_distance",
    "luminosity_distance",
    "comoving_transverse_distance",
    "angular_diameter_distance",
    "comoving_distance",
    "Hubble_parameterization",
]


##############################################################################
# IMPORTS

import typing as T

import numpy as np

from .distance import (
    angular_diameter_distance,
    comoving_distance,
    comoving_transverse_distance,
    hubble_distance,
    luminosity_distance,
    scale_factor,
)

###############################################################################
# CODE
###############################################################################

#####################################################################
# Magnitudes


def Hubble_parameterization(
    z: T.Sequence, omM: float, omR: float, omL: float, omK: float = 0.0
) -> T.Sequence:
    r"""Hubble parameterization.

    ::

        E(z) = H(z) / H0
             = sqrt(omL + omK(1+z)^2 + omM(1+z)^3 + omR(1+z)^4)

    Parameters
    ----------
    DM: array_like
        Distance Modulus value, in [magnitude].
        If Quantity, pass as ``DM.to_value('mag')``

    Returns
    -------
    d: scalar, array
        Luminosity distance value in parsecs.

    References
    ----------
    https://en.wikipedia.org/wiki/Distance_measures_(cosmology)

    """
    return np.sqrt(
        omL
        + omK * (1 + z) ** 2.0
        + omM * (1 + z) ** 3.0
        + omR * (1 + z) ** 4.0
    )


# /def


###############################################################################
# END

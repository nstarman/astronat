# -*- coding: utf-8 -*-

"""Distance related functions.

These functions are implemented here in minimal form for low overhead.
For equivalent functions with the full set of bells and whistles, see
:mod:`~astronat.common`.

"""

__author__ = "Nathaniel Starkman"


##############################################################################
# IMPORTS

import typing as T

import numba
import numpy as np

##############################################################################
# PARAMETERS

__all__: T.List[str] = [
    "apparent_to_absolute_magnitude",
    "absolute_to_apparent_magnitude",
]


###############################################################################
# CODE
###############################################################################

#####################################################################
# Magnitudes


@numba.njit
def apparent_to_absolute_magnitude(m: T.Sequence, d: T.Sequence) -> T.Sequence:
    """Calculate the absolute magnitude.

    M = m - 5 log10(d) + 5

    Parameters
    ----------
    m: array_like
        apparent magnitude
    d: array_like
        Luminosity distance to object in [parsec].
        If Quantity, pass as ``d.to_value('pc')``

    Returns
    -------
    M: array_like
        absolute magnitudes

    """
    return m - 5.0 * np.log10(d) + 5.0


# /def


@numba.njit
def absolute_to_apparent_magnitude(M: T.Sequence, d: T.Sequence) -> T.Sequence:
    """Calculate the apparent magnitude.

        m = M + 5 log10(d) - 5

    Parameters
    ----------
    M: array_like
        absolute magnitude
    d: array_like
        Luminosity distance to object in [parsec].
        If Quantity, pass as ``d.to_value('pc')``

    Returns
    -------
    m: array_like
        apparent magnitudes

    """
    return M + 5.0 * np.log10(d) - 5.0


# /def


###############################################################################
# END

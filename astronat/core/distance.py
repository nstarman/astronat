# -*- coding: utf-8 -*-

"""Distance related functions.

These functions are implemented here in minimal form for fast compute.
For equivalent functions with the full set of bells and whistles, see
:mod:`~astronat.common`.

The units are natural units.

"""

##############################################################################
# IMPORTS

import typing as T

import numba
import numpy as np

from ..constants import default_values as constants

##############################################################################
# PARAMETERS

__all__: T.List[str] = [
    "scale_factor",
    "hubble_distance",
    "luminosity_distance",
    "comoving_transverse_distance",
    "angular_diameter_distance",
    "comoving_distance",
    "distance_to_modulus",
    "modulus_to_distance",
    "parallax_angle",
    "parallax_distance",
    "max_angular_separation",
]


_AU_to_pc = constants.AU_to_pc
_m_to_pc = 3.2407792700054e-17

pi = np.pi

# c_in_kms = constants.c / 1000.0

H0 = 70  # TODO make a configurable parameter
"""Hubble constant. Value in [km / s / Mpc]."""


###############################################################################
# CODE
###############################################################################


@numba.njit
def scale_factor(z: float) -> float:
    """Scale factor from redshift.

    Parameters
    ----------
    z : float
        Redshift.

    Returns
    -------
    a : float
        Scale factor.

    """
    return 1 / (1 + z)


# /def


# -------------------------------------------------------------------


@numba.njit
def hubble_distance(H0: float = H0) -> float:
    """Hubble distance.

    Parameters
    ----------
    H0 : float
        Hubble constant. value in [km / s / Mpc].

    Returns
    -------
    dH : float
        Hubble distance. Value in [parsec]

    """
    return constants.c / H0 / 1e9


# /def

# -------------------------------------------------------------------


@numba.njit
def luminosity_distance(L: T.Sequence, F: T.Sequence) -> T.Sequence:
    """Luminosity distance from luminosity and flux.

        dL^2 = L / 4 pi F

    Since distances can only be measured observationally, The luminosity
    distance is the fundamental distance measure in Astronomy.

    Parameters
    ----------
    L : array_like
        Luminosity. value in [W]
    F : array_like
        Flux. value in [W / m^2]
        If Quantity, pass as ``distance.to_value('pc')`

    Returns
    -------
    dL: scalar, array
        Luminosity distance in [parsec]

    """
    return np.sqrt(L / (4 * pi * F)) * _m_to_pc


# /def

# -------------------------------------------------------------------


@numba.njit
def comoving_transverse_distance(dL: T.Sequence, z: T.Sequence) -> T.Sequence:
    """Comoving Transverse Distance.

        D_L = (1 + z) D_M

    Parameters
    ----------
    dL : array_like
        Luminosity distance.
    z : array_like
        Redshift.

    Returns
    -------
    dM: scalar, array
        Comoving transverse distance. value in units of dL.

    References
    ----------
    https://en.wikipedia.org/wiki/Luminosity_distance

    """
    return dL / (1 + z)


# /def

# -------------------------------------------------------------------


@numba.njit
def angular_diameter_distance(dL: T.Sequence, z: T.Sequence) -> T.Sequence:
    """Angular diameter distance (Etherington's distance).

        D_L = (1 + z)^2 D_A

    The Etherington's distance-duality equation is the relationship between
    the luminosity distance of standard candles and the angular diameter
    distance.[1] The equation is as follows:  :math:`d_{L}=(1+z)^{2}d_{A}}`,
    where `dL` is the luminosity distance and dA the angular-diameter
    distance.

    Parameters
    ----------
    dL : array_like
        Luminosity distance.
    z : array_like
        Redshift.

    Returns
    -------
    dA: scalar, array
        Angular diameter distance. value in units of dL.

    References
    ----------
    https://en.wikipedia.org/wiki/Luminosity_distance

    """
    return dL / (1 + z) ** 2


etherington_distance = angular_diameter_distance
# /def


# -------------------------------------------------------------------


@numba.njit
def comoving_distance(dM: T.Sequence, H0: float = H0, omK: float = 0.0):
    """Comoving distance from comoving transverse distance.

    To use the base distance, the luminosity distance,
    ::

        dM = comoving_transverse_distance(dL, z)

    Parameters
    ----------
    dM : array_like
        Luminosity distance. Value in [parsec].
    H0 ; float
        Hubble constant. Value in [km / s / Mpc]
    omK : float
        Omega curvature

    Returns
    -------
    dC : T.Sequence
        Comoving distance. Value in [parsec]

    """
    if omK == 0.0:
        dC = dM

    elif omK > 0.0:
        dC = (hubble_distance(H0) / np.sqrt(omK)) * np.arcsinh(
            dM * np.sqrt(omK) / hubble_distance(H0)
        )

    else:  # omK < 0
        dC = (hubble_distance(H0) / np.sqrt(-omK)) * np.arcsin(
            dM * np.sqrt(-omK) / hubble_distance(H0)
        )

    return dC


# /def


#####################################################################
# Distance Modulus


@numba.njit
def distance_to_modulus(d: T.Sequence) -> T.Sequence:
    """Distance Modulus from the distance.

        DM = 5 log10(d / 10 pc)

    Parameters
    ----------
    d: array_like
        Luminosity distance value in [parsec].
        If Quantity, pass as ``distance.to_value('pc')`

    Returns
    -------
    DM: scalar, array
        Distance Modulus value in magnitudes

    """
    return 5.0 * np.log10(d) - 5.0


# /def


@numba.njit
def modulus_to_distance(DM: T.Sequence) -> T.Sequence:
    r"""Distance Modulus.

        d = 10^((DM + 5) / 5)

    Parameters
    ----------
    DM: array_like
        Distance Modulus value, in [magnitude].
        If Quantity, pass as ``DM.to_value('mag')``

    Returns
    -------
    d: scalar, array
        Luminosity distance value in parsecs.

    """
    return np.power(10.0, np.divide(DM, 5.0) + 1)


# /def


###############################################################################
# Parallax


@numba.njit
def parallax_angle(d: T.Sequence) -> T.Sequence:
    """Compute parallax angle from distance.

        arctan(1 [AU] / d [pc])

    Parameters
    ----------
    d: array_like
        Distance value in [parsec].
        If Quantity, pass as ``d.to_value('pc')``

    Returns
    -------
    angle: array_like
        Parallax angle in [radians].

    """
    return np.arctan(_AU_to_pc / d)


# /def


@numba.njit
def parallax_distance(p: T.Sequence) -> T.Sequence:
    """Compute parallax distance from angle.

        1 [pc] / tan(p)

    Parameters
    ----------
    p: array_like
       Parallax angle value in [radians].

    Returns
    -------
    d: array_like
        Distance value in [parsec].

    """
    return _AU_to_pc / np.tan(p)


# /def


###############################################################################
# Angular Separation


@numba.njit
def max_angular_separation(d: T.Sequence, doff: T.Sequence) -> T.Sequence:
    """Maximum angular separation.

    Parameters
    ----------
    d: array_like
        Ddistance to original coordinate
    doff: array_like
        Distance offset from original coordinate
        Value in same units as `d`.

    Returns
    -------
    angle: array_like
        maximum angular separation, value in [radian].

    """
    return np.fabs(np.arctan(np.divide(doff, d)))


# /def


###############################################################################
# END

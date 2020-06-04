# -*- coding: utf-8 -*-

"""Core set of astronomy functions.

These functions are implemented here in minimal form for low overhead.
For equivalent functions with the full set of bells and whistles, see
:mod:`~astronat.common`.

"""

__author__ = "Nathaniel Starkman"


##############################################################################
# IMPORTS

# BUILT-IN

import typing as T


# THIRD PARTY

import numpy as np
import numba

# PROJECT-SPECIFIC


from ..constants import default_values as constants


##############################################################################
# PARAMETERS

__all__: T.List[str] = [
    # "distance_to_modulus",
    # "modulus_to_distance",
    # "apparent_to_absolute_magnitude",
    # "absolute_to_apparent_magnitude",
]


# _AU_to_pc = constants.AU_to_pc
# pi = np.pi


###############################################################################
# CODE
###############################################################################

#####################################################################
# General


@numba.njit
def lens_power(f: T.Sequence) -> T.Sequence:
    """Lens power from focal length.

    Parameters
    ----------
    f: array_like
        Focal length value in [meters].

    Returns
    -------
    P: array_like
        Lens power, value in [1 / meter].

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return 1 / f


# /def


@numba.njit
def lateral_magnification(x1: float, x2: float) -> float:
    """Lateral magnification.

    ::

        m = -x2 / x1 = y2 / y1

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return -x2 / x1


# /def


@numba.njit
def angular_magnification(theta1: float, theta2: float) -> float:
    """Lateral magnification.

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return theta2 / theta1


# /def


#####################################################################
# Radiometry


@numba.njit
def radiosity(E: T.Sequence, M: T.Sequence) -> T.Sequence:
    """Lateral magnification.

    Returns
    -------
    radiosity : array_like
        Value in [W / m^2]

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return E + M


# /def


#####################################################################
# Luminal ElectroMagnetic Waves


@numba.njit
def dielectric_mean_energy_density(
    E: T.Sequence, B: T.Sequence, permittivity: float, permeability: float
) -> T.Sequence:
    """Dielectric Energy Density.

    .. todo::

        have defaults for permittivity and permeability.

    Returns
    -------
    radiosity : array_like
        Value in [W / m^2]

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return (permittivity * E ** 2 + permeability * B ** 2) / 2


# /def


# @numba.njit
# def irradiance():


@numba.njit
def doppler_wavelength(lambda0: T.Sequence, v: float):

    return lambda0 * np.sqrt((constants.c - v) / constants.c + v)


# /def


@numba.njit
def doppler_velocity(Lambda, lambda0: T.Sequence):

    return np.fabs(Lambda - lambda0) * constants.c / lambda0


# /def


@numba.njit
def Cherenkov_angle(n, v):

    return np.arccos(constants.c / n / v)


# /def


#####################################################################
# Geometric Optics


@numba.njit
def refractive_index_from_material(
    permittivity: T.Sequence, permeability: T.Sequence
) -> T.Sequence:
    """Refractive index from material properties.

    Parameters
    ----------
    permittivity : array_like
        Value in [Farad / meter].

        The absolute permittivity, often simply called permittivity and
        denoted by the Greek letter epsilon, is a measure of the electric
        polarizability of a dielectric. A material with high permittivity
        polarizes more in response to an applied electric field than a
        material with low permittivity, thereby storing more energy in the
        electric field. In electrostatics, the permittivity plays an important
        role in determining the capacitance of a capacitor [#]_.

    permeability : array_like
        Value in [Newton / Ampere squared].

        permeability is the measure of the resistance of a material against
        the formation of a magnetic field, otherwise known as distributed
        inductance in transmission line theory. Hence, it is the degree of
        magnetization that a material obtains in response to an applied
        magnetic field. Magnetic permeability is typically represented by the
        (italicized) Greek letter mu [#]_.

    Returns
    -------
    n: array_like
        Refractive index. [unitless].

    References
    ----------
    .. [#] https://en.wikipedia.org/wiki/Permittivity
    .. [#] https://en.wikipedia.org/wiki/Permeability_(electromagnetism)

    """

    return np.sqrt(permittivity * permeability)


# /def


@numba.njit
def critical_angle(n1: T.Sequence, n2: T.Sequence) -> T.Sequence:
    """Critical angle from refractive indices.

    Parameters
    ----------
    n1, n2: array_like
        Refractive index. [unitless].

    Returns
    -------
    theta_c: array_like
        critical angle, in [radians]

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return np.arcsin(n2 / n1)


# /def


@numba.njit
def critical_angle_from_velocities(
    v1: T.Sequence, v2: T.Sequence
) -> T.Sequence:
    """Critical angle from light velocity in media.

    A monochromatic beam will have different velocities in different media.

    Parameters
    ----------
    v1, v2: array_like
        Light velocity. Inversely proportional to refractive index.
        Values from same units.

    Returns
    -------
    theta_c: array_like
        critical angle, in [radians].

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return critical_angle(v2, v1)  # Inverse propto to n.


# /def


@numba.njit
def critical_angle_from_wavelength(
    lambda1: T.Sequence, lambda2: T.Sequence
) -> T.Sequence:
    """Critical angle from wavelengths in media.

    A monochromatic beam will have its wavelength altered in different media.

    Parameters
    ----------
    lambda1, lambda2: array_like
        The wavelength in medias 1 and 2, respectively.
        Inversely proportional to refractive index.
        Values from same units.

    Returns
    -------
    theta_c: array_like
        critical angle, in [radians].

    References
    ----------
    https://en.wikipedia.org/wiki/List_of_photonics_equations

    """
    return critical_angle(lambda2, lambda1)  # Inverse propto to n.


# /def


###############################################################################
# END

# -*- coding: utf-8 -*-

"""Modeling with :mod:`~astropy.modeling`."""

__author__ = "Nathaniel Starkman"


__all__ = [
    # functions
    "cartesian_to_spherical",
    "cartesian_to_coslatspherical",
    # "_wrap_model_output",
]


##############################################################################
# IMPORTS

import astropy.coordinates as coord
import astropy.units as u
from astropy.modeling import FittableModel

##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################

# class _wrap_model_output:
#     """Expose wrapped value as "item", for non-Sequence Model returns."""

#     def __init__(self, value):
#         self.value = value

#     def item(self):
#         return self.value


# # /class


#####################################################################


class cartesian_to_spherical(FittableModel):
    """Model Converting Cartesian to Spherical Coordinates.

    .. todo::

        For speed, use the actual mechanics of transformation, not the Astropy
        representation classes.


    """

    # inputs
    n_inputs = 6

    # outputs
    n_outputs = 6
    outputs = ("lon", "lat", "distance", "d_lon", "d_lat", "d_distance")
    return_units = {
        "lon": u.deg,
        "lat": u.deg,
        "distance": u.kpc,
        "d_lon": u.mas / u.yr,
        "d_lat": u.mas / u.yr,
        "d_distance": u.km / u.s,
    }

    @staticmethod
    def evaluate(
        x: u.kpc,
        y: u.kpc,
        z: u.kpc,
        v_x: u.km / u.s,
        v_y: u.km / u.s,
        v_z: u.km / u.s,
    ):
        # TODO using actual mechanics of transformation
        d_xyz = coord.CartesianDifferential(d_x=v_x, d_y=v_y, d_z=v_z)
        xyz = coord.CartesianRepresentation(x=x, y=y, z=z, differentials=d_xyz)

        sph = xyz.represent_as(
            coord.SphericalRepresentation,
            differential_class=coord.SphericalDifferential,
        )
        d_sph = sph.differentials["s"]

        return (
            sph.lon,
            sph.lat,
            sph.distance,
            d_sph.d_lon,
            d_sph.d_lat,
            d_sph.d_distance,
        )

    # /def


# /class


# -------------------------------------------------------------------


class cartesian_to_coslatspherical(FittableModel):
    """Model Converting Cartesian to Spherical Coordinates.

    .. todo::

        For speed, use the actual mechanics of transformation, not the Astropy
        representation classes.

    """

    # inputs
    n_inputs = 6

    # outputs
    n_outputs = 6
    outputs = ("lon", "lat", "distance", "d_lon", "d_lat", "d_distance")
    return_units = {
        "lon": u.deg,
        "lat": u.deg,
        "distance": u.kpc,
        "d_lon_coslat": u.mas / u.yr,
        "d_lat": u.mas / u.yr,
        "d_distance": u.km / u.s,
    }

    @staticmethod
    def evaluate(
        x: u.kpc,
        y: u.kpc,
        z: u.kpc,
        v_x: u.km / u.s,
        v_y: u.km / u.s,
        v_z: u.km / u.s,
    ):
        # TODO using actual mechanics of transformation
        d_xyz = coord.CartesianDifferential(d_x=v_x, d_y=v_y, d_z=v_z)
        xyz = coord.CartesianRepresentation(x=x, y=y, z=z, differentials=d_xyz)

        sph = xyz.represent_as(
            coord.SphericalRepresentation,
            differential_class=coord.SphericalCosLatDifferential,
        )
        d_sph = sph.differentials["s"]

        return (
            sph.lon,
            sph.lat,
            sph.distance,
            d_sph.d_lon_coslat,
            d_sph.d_lat,
            d_sph.d_distance,
        )

    # /def


# /class


##############################################################################
# END

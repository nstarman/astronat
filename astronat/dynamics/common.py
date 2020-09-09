# -*- coding: utf-8 -*-

"""Common Code. Defines usefule Types and conversion functions."""

__all__ = [
    # Types
    "PotentialType",
    "GalpyBasePotentialType",
    "GalpyPotentialType",
    "GalaPotentialType",
    # conversions
    "galpy_Orbit_to_SkyCoord",
    "galpy_Orbit_to_gala_Orbit",
    "gala_Orbit_to_galpy_Orbit",
]


##############################################################################
# IMPORTS

import typing as T

import gala.potential as galapot
from astropy import coordinates as coords
from gala.dynamics import Orbit as galaOrbit
from galpy.orbit import Orbit as galpyOrbit
from utilipy.data_utils import DataTransform
from utilipy.utils.typing import FrameOptionsType

from ..common import data_graph
from .potential import GalpyBasePotentialType, GalpyPotentialType

##############################################################################
# PARAMETERS

GalaPotentialType = T.TypeVar(
    "GalaPotentialType",
    galapot.PotentialBase,
    galapot.CompositePotential,
    galapot.CCompositePotential,
)

PotentialType = T.TypeVar(
    "PotentialType", GalaPotentialType, GalpyPotentialType
)


##############################################################################
# CODE
##############################################################################


@data_graph.register(DataTransform, galpyOrbit, coords.SkyCoord)
def galpy_Orbit_to_SkyCoord(data, frame: FrameOptionsType = coords.ICRS, **kw):
    """Transform galpy Orbit to Astropy SkyCoord.

    .. todo::

        - galpy obs=[X,Y,Z], ro= distance, vo= velocity

    Parameters
    ----------
    data : :class:`~galpy.orbit.Orbit` instance
    frame : str or `BaseCoordinateFrame` or `SkyCoord`, optional
        frame to which to transform the SkyCoord

    Returns
    -------
    :class:`~astropy.coordinates.SkyCoord`
        in frame `frame`

    """
    sc = data.SkyCoord(data.time())
    sc = sc.transform_to(frame)

    return sc


# /def


# -------------------------------------------------------------------


@data_graph.register(DataTransform, galpyOrbit, galaOrbit)
def galpy_Orbit_to_gala_Orbit(data, frame: FrameOptionsType = "icrs", **kw):
    """Transform galpy Orbit to gala Orbit.

    .. todo::

        - galpy obs=[X,Y,Z], ro= distance, vo= velocity

    Parameters
    ----------
    data : :class:`~galpy.orbit.Orbit` instance
    frame : str or `BaseCoordinateFrame` or `SkyCoord`, optional
        frame to which to transform the SkyCoord

    Returns
    -------
    :class:`~gala.dynamics.Orbit`
        in frame `frame`

    """
    raise Exception("Not yet implemented")

    # return sc


# /def


@data_graph.register(DataTransform, galaOrbit, galpyOrbit)
def gala_Orbit_to_galpy_Orbit(data, frame: FrameOptionsType = "icrs", **kw):
    """Transform gala Orbit to galpy Orbit.

    .. todo::

        - galpy obs=[X,Y,Z], ro= distance, vo= velocity

    Parameters
    ----------
    data : :class:`~gala.dynamics.Orbit` instance
    frame : str or `BaseCoordinateFrame` or `SkyCoord`, optional
        frame to which to transform the SkyCoord

    Returns
    -------
    :class:`~galpy.orbit.Orbit`
        in frame `frame`

    """
    raise Exception("Not yet implemented")

    # return sc


# /def


##############################################################################
# END

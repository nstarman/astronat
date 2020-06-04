# -*- coding: utf-8 -*-

"""Data-type conversion functions, registered into data TransformGraph."""

__author__ = "Nathaniel Starkman"


# __all__ = [
#     ""
# ]


##############################################################################
# IMPORTS

# BUILT-IN

# THIRD PARTY

from astropy import coordinates as coords

from utilipy.data_utils import DataTransform
from utilipy.utils.typing import FrameOptionsType


# PROJECT-SPECIFIC

from ..common import data_graph


##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


# TODO move to utilipy itself?
@data_graph.register(
    DataTransform, coords.SkyCoord, coords.BaseCoordinateFrame
)
def SkyCoord_to_Frame(data, frame: FrameOptionsType = None, **kw):
    """Transform Astropy SkyCoord to Coordinate Frame.

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
    if frame is not None:
        data = data.transform_to(frame)

    return data.frame


# /def


# -------------------------------------------------------------------


# TODO move to utilipy itself?
@data_graph.register(
    DataTransform, coords.BaseCoordinateFrame, coords.BaseCoordinateFrame,
)
def Frame_to_Frame(data, frame: FrameOptionsType = None, **kw):
    """Transform Astropy SkyCoord to Coordinate Frame.

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
    if frame is not None:
        data = data.represent_as(frame)

    return data


# /def


##############################################################################
# END

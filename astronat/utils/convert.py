# -*- coding: utf-8 -*-

"""Data-type conversion functions, registered into data TransformGraph."""

__author__ = "Nathaniel Starkman"


# __all__ = [
#     ""
# ]


##############################################################################
# IMPORTS

import astropy.coordinates as coord
import numpy as np
from astropy.table import Table
from utilipy.data_utils import DataTransform
from utilipy.utils.typing import FrameOptionsType, TableType

from ..common import data_graph

##############################################################################
# CODE
##############################################################################


# TODO move to utilipy itself?
@data_graph.register(DataTransform, coord.SkyCoord, coord.BaseCoordinateFrame)
def SkyCoord_to_Frame(data, frame: FrameOptionsType = None, **kw):
    """Transform Astropy SkyCoord to Coordinate Frame.

    Parameters
    ----------
    data : `SkyCoord`
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
    DataTransform, coord.BaseCoordinateFrame, coord.BaseCoordinateFrame,
)
def Frame_to_Frame(
    data, frame: FrameOptionsType = None, **kw
) -> FrameOptionsType:
    """Transform Astropy SkyCoord to Coordinate Frame.

    Parameters
    ----------
    data : `BaseCoordinateFrame`
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

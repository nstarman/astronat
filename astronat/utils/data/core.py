# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Data-Graph.

.. todo::

    copy `data_graph`

"""

__all__ = [
    "xmatch_data_graph",
]


##############################################################################
# IMPORTS

# THIRD PARTY

from astropy.coordinates import (
    SkyCoord,
    BaseCoordinateFrame as BCFrame,
)
from astropy.table import Table
from astropy.time import Time

from utilipy.data_utils.xfm import (
    data_graph as xmatch_data_graph,
    DataTransform,
)


##############################################################################
# DATA TRANFSFORM REGISTRY OVERRRIDES


@xmatch_data_graph.register(DataTransform, Table, SkyCoord)
def Table_to_BaseCoordinateFrame(data, obstime=None):
    """`~Table` to `BaseCoordinateFrame`.

    .. todo::

        - first try to determine if a SkyCoord is embedded in the table
        - more robust method of determining epoch

    """
    frame = SkyCoord.guess_from_table(data)

    if obstime is None:

        if "obstime" in data.dtype.fields:
            frame.obstime = Time(data["obstime"])
        elif "epoch" in data.dtype.fields:
            frame.obstime = Time(data["epoch"])
        elif "ref_epoch" in data.dtype.fields:
            frame.obstime = Time(data["epoch"])

        elif "obstime" in data.meta:
            frame.obstime = Time(data.meta["obstime"])
        elif "epoch" in data.meta:
            frame.obstime = Time(data.meta["epoch"])
        elif "ref_epoch" in data.meta:
            frame.obstime = Time(data.meta["epoch"])

        else:
            print("Could not determine epoch")

    else:
        frame.obstime = obstime

    return frame


# /def


@xmatch_data_graph.register(DataTransform, BCFrame, SkyCoord)
def BasecoordinateFrame_to_SkyCoord(data, obstime=None):
    """Transform from BaseCoordinateFrame to SkyCoord."""
    out = SkyCoord(data)
    if hasattr(out, "obstime"):
        if obstime is not None:
            out.obstime = obstime
    else:  # TODO do if obstime is None?
        out.obstime = obstime

    return out


# /def


##############################################################################
# END

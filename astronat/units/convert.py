# -*- coding: utf-8 -*-

"""Unit Conversions."""

__author__ = "Nathaniel Starkman"


__all__ = ["from_amuse", "hms_str_to_unit"]


##############################################################################
# IMPORTS

import typing as T

import astropy.coordinates as coord
import astropy.units as u
import numpy as np

from . import Unit

# from .decorators import quantity_output

###############################################################################
# CODE
###############################################################################


def from_amuse(quantity):
    """Convert AMUSE quantity to astropy quantity.

    only amuse quantities need to be converted.
    floats and astropy quantities are left as is.

    Parameters
    ----------
    quantity: AMUSE quantity or array-like

    Returns
    -------
    quantity: astropy quantity or array-like

    Notes
    -----
    requires that amuse units are string represented in astropy

    """
    try:
        quantity.value_in
    except Exception:
        out = quantity
    else:
        out: u.Quantity
        out = quantity.value_in(quantity.unit) * Unit(str(quantity.unit))

    return out


# /def


###############################################################################


# @quantity_output(unit=u.deg)
def hms_str_to_unit(hms: T.Sequence[str], empty_to_=np.NaN, copy=False):
    """Change a hms string to Astropy deg units.

    Parameters
    ----------
    hms : str
        hour-minute-second array of strings

    Returns
    -------
    Angle

    """
    return coord.Angle(hms, unit="hourangle", copy=copy)


# /def


###############################################################################
# END

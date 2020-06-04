# -*- coding: utf-8 -*-

"""Unit Conversions."""

__author__ = "Nathaniel Starkman"


__all__ = ["from_amuse", "hms_str_to_unit"]


##############################################################################
# IMPORTS

# BUILT-IN

import typing as T


# THIRD PARTY

import astropy.coordinates as coord
import astropy.units as u


# PROJECT-SPECIFIC

from . import Unit
from .decorators import quantity_output


###############################################################################
# PARAMETERS

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


@quantity_output(unit=u.deg)
def hms_str_to_unit(hms: T.Sequence[str]):
    """Change a hms string to Astropy deg units.

    Parameters
    ----------
    hms : str
        hour-minute-second array of strings

    """
    return coord.Angle(hms, unit="hourangle")


# /def


###############################################################################
# END

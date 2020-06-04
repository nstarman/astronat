# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : unit decorators
# PROJECT : astronat
#
# ----------------------------------------------------------------------------

"""Conversions."""

__author__ = "Nathaniel Starkman"


__all__ = ["from_amuse"]


##############################################################################
# IMPORTS

# PROJECT-SPECIFIC

from . import Unit


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
        out = quantity.value_in(quantity.unit) * Unit(str(quantity.unit))

    return out


# /def


###############################################################################


###############################################################################
# END

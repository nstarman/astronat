# -*- coding: utf-8 -*-

"""Tests for :mod:`~astronat.utils.table.utils`."""

__all__ = [
    "test_rename_columns",
    "test_cast_columns",
]


##############################################################################
# IMPORTS

import astropy.units as u
from astropy.table import QTable

from astronat.utils.table import utils

##############################################################################
# PARAMETERS

tbl = QTable([[2.0, 5.0], ["x", "y"]], names=("a", "b"))

##############################################################################
# CODE
##############################################################################


def test_rename_columns():
    """Test `~astronat.utils.table.utils.rename_columns`."""
    utils.rename_columns(table=tbl, rename={"a": "A"})

    assert tbl.colnames == ["A", "b"]

    utils.rename_columns(table=tbl, rename={"A": "a"})

    assert tbl.colnames == ["a", "b"]


# /def

# -------------------------------------------------------------------


def test_cast_columns():
    """Test `~astronat.utils.table.utils.rename_columns`."""
    utils.cast_columns(table=tbl, recast={"a": lambda x: x * u.km})

    assert all(tbl["a"] == [2.0, 5.0] * u.km)


# /def


##############################################################################
# END

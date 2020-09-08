# -*- coding: utf-8 -*-

"""Test `~astronat.phot.MegaCam_PanSTARRS.data`."""

__all__ = ["test_MegaCamGen1_from_PS1"]


##############################################################################
# IMPORTS

from astropy.table import QTable

from .. import data

##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


def test_MegaCamGen1_from_PS1():
    """Test `~astronat.phot.MegaCam_PanSTARRS.data.read_MegaCamGen1_from_PS1`.

    .. todo::

        finish

    """
    df = data.read_MegaCamGen1_from_PS1()

    assert isinstance(df, QTable)

    return


# /def


# ------------------------------------------------------------------------


##############################################################################
# END

# -*- coding: utf-8 -*-

"""Combining relevant functions."""

__author__ = "Nathaniel Starkman"
__credits__ = [
    "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html"
]


__all__ = [
    # Mixed_MegaCamGen3_PS1
    "mixed",
    # MegaCamGen3_from_PS1
    "U_MP9302",
    "u_MC",
    "G_MP9402",
    "g_MC",
    "R_MP9602",
    "r_MC",
    "I_MP9703",
    "i_MC",
    "Z_MP9901",
    "z_MC",
    "GRI_MP9605",
    "gri_MC",
]


#############################################################################
# IMPORTS

from . import Mixed_MegaCamGen3_PS1 as mixed
from .MegaCamGen3_from_PS1 import G_MP9402 as g_MC
from .MegaCamGen3_from_PS1 import GRI_MP9605 as gri_MC
from .MegaCamGen3_from_PS1 import I_MP9703 as i_MC
from .MegaCamGen3_from_PS1 import R_MP9602 as r_MC
from .MegaCamGen3_from_PS1 import U_MP9302 as u_MC
from .MegaCamGen3_from_PS1 import Z_MP9901 as z_MC

# from .PS1_from_MegaCamGen3 import


#############################################################################
# END

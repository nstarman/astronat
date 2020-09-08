# -*- coding: utf-8 -*-

"""MegaCam Generation 1 and panSTARRS generation 1.

    The conversions are sourced from [1]_

Modules
-------
Mixed_MegaCamGen1_PS1
MegaCamGen1_from_PS1:
PS1_from_MegaCamGen1:

References
----------
.. [1] http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html

"""

__author__ = "Nathaniel Starkman"
__credits__ = [
    "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html"
]


__all__ = [
    # Mixed_MegaCamGen1_PS1
    "mixed",
    # MegaCamGen1_from_PS1
    "U_MP9301",
    "u_MC",
    "G_MP9401",
    "g_MC",
    "R_MP9601",
    "r_MC",
    "I_MP9701",
    "i_MC",
    "Z_MP9801",
    "z_MC",
    "umg_MC",
    "umr_MC",
    "umi_MC",
    "umz_MC",
    "gmr_MC",
    "gmi_MC",
    "gmz_MC",
    "rmi_MC",
    "rmz_MC",
    "imz_MC",
    # PS1_from_MegaCamGen1
    "g_PS",
    "r_PS",
    "i_PS",
    "z_PS",
    "gmr_PS",
    "gmi_PS",
    "gmz_PS",
    "rmi_PS",
    "rmz_PS",
    "imz_PS",
]


#############################################################################
# IMPORTS

from . import Mixed_MegaCamGen1_PS1 as mixed
from .MegaCamGen1_from_PS1 import G_MP9401 as g_MC
from .MegaCamGen1_from_PS1 import I_MP9701 as i_MC
from .MegaCamGen1_from_PS1 import R_MP9601 as r_MC
from .MegaCamGen1_from_PS1 import U_MP9301 as u_MC
from .MegaCamGen1_from_PS1 import Z_MP9801 as z_MC
from .MegaCamGen1_from_PS1 import GmI as gmi_MC
from .MegaCamGen1_from_PS1 import GmR as gmr_MC
from .MegaCamGen1_from_PS1 import GmZ as gmz_MC
from .MegaCamGen1_from_PS1 import ImZ as imz_MC
from .MegaCamGen1_from_PS1 import RmI as rmi_MC
from .MegaCamGen1_from_PS1 import RmZ as rmz_MC
from .MegaCamGen1_from_PS1 import UmG as umg_MC
from .MegaCamGen1_from_PS1 import UmI as umi_MC
from .MegaCamGen1_from_PS1 import UmR as umr_MC
from .MegaCamGen1_from_PS1 import UmZ as umz_MC
from .PS1_from_MegaCamGen1 import G as g_PS
from .PS1_from_MegaCamGen1 import GmI as gmi_PS
from .PS1_from_MegaCamGen1 import GmR as gmr_PS
from .PS1_from_MegaCamGen1 import GmZ as gmz_PS
from .PS1_from_MegaCamGen1 import I_band as i_PS
from .PS1_from_MegaCamGen1 import ImZ as imz_PS
from .PS1_from_MegaCamGen1 import R as r_PS
from .PS1_from_MegaCamGen1 import RmI as rmi_PS
from .PS1_from_MegaCamGen1 import RmZ as rmz_PS
from .PS1_from_MegaCamGen1 import Z as z_PS

#############################################################################
# END

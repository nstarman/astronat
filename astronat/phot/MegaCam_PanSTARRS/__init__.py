# -*- coding: utf-8 -*-

"""MegaCam_Panstarrs Filter Transformations."""

__author__ = "Nathaniel Starkman"
__credits__ = [
    "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html"
]

__all__ = [
    "MCgen1_PS1",
    "MCgen2_PS1",
    "MCgen3_PS1",
    # General
    "UmG",
    "UmR",
    "UmI",
    "UmZ",
    "GmR",
    "GmI",
    "GmZ",
    "RmI",
    "RmZ",
    "ImZ",
]

# __all__ += [
#     # MegaCamGen1_from_PS1
#     "U_MP9301",
#     "G_MP9401",
#     "R_MP9601",
#     "I_MP9701",
#     "Z_MP9801",
#     "umg_MC1PS1",
#     "umr_MC1PS1",
#     "umi_MC1PS1",
#     "umz_MC1PS1",
#     "gmr_MC1PS1",
#     "gmi_MC1PS1",
#     "gmz_MC1PS1",
#     "rmi_MC1PS1",
#     "rmz_MC1PS1",
#     "imz_MC1PS1",
#     # PS1_from_MegaCamGen1
#     "g_PS1MC1",
#     "r_PS1MC1",
#     "i_PS1MC1",
#     "z_PS1MC1",
#     "gmr_PS1MC1",
#     "gmi_PS1MC1",
#     "gmz_PS1MC1",
#     "rmi_PS1MC1",
#     "rmz_PS1MC1",
#     "imz_PS1MC1",
#     # mixed
#     "mixed_MC1PS1",
#     "mixed_PS1MC1",
#     # MegaCamGen2_from_PS1
#     "I_MP9702",
#     # PS1_from_MegaCamGen2
#     # mixed
#     "mixed_MC2PS1",
#     "mixed_PS1MC2",
#     # MegaCamGen3_from_PS1
#     "U_MP9302",
#     "G_MP9402",
#     "R_MP9602",
#     "I_MP9703",
#     "Z_MP9901",
#     "GRI_MP9605",
#     # PS1_from_MegaCamGen3
#     # mixed
#     "mixed_MC3PS1",
#     "mixed_PS1MC3",
# ]


#############################################################################
# IMPORTS

from typing import Any

from astropy.table import Table
from typing_extensions import Literal

from astronat import units

from . import MCgen1_PS1, MCgen2_PS1, MCgen3_PS1

# Not in __all__
# TODO deprecate, but needed in below functions
from .MCgen1_PS1 import g_PS as g_PS1MC1
from .MCgen1_PS1 import gmi_MC as gmi_MC1PS1
from .MCgen1_PS1 import gmi_PS as gmi_PS1MC1
from .MCgen1_PS1 import gmr_MC as gmr_MC1PS1
from .MCgen1_PS1 import gmr_PS as gmr_PS1MC1
from .MCgen1_PS1 import gmz_MC as gmz_MC1PS1
from .MCgen1_PS1 import gmz_PS as gmz_PS1MC1
from .MCgen1_PS1 import i_PS as i_PS1MC1
from .MCgen1_PS1 import imz_MC as imz_MC1PS1
from .MCgen1_PS1 import imz_PS as imz_PS1MC1
from .MCgen1_PS1 import mixed as mixed_MC1PS1
from .MCgen1_PS1 import mixed as mixed_PS1MC1
from .MCgen1_PS1 import r_PS as r_PS1MC1
from .MCgen1_PS1 import rmi_MC as rmi_MC1PS1
from .MCgen1_PS1 import rmi_PS as rmi_PS1MC1
from .MCgen1_PS1 import rmz_MC as rmz_MC1PS1
from .MCgen1_PS1 import rmz_PS as rmz_PS1MC1
from .MCgen1_PS1 import umg_MC as umg_MC1PS1
from .MCgen1_PS1 import umi_MC as umi_MC1PS1
from .MCgen1_PS1 import umr_MC as umr_MC1PS1
from .MCgen1_PS1 import umz_MC as umz_MC1PS1
from .MCgen1_PS1 import z_PS as z_PS1MC1
from .MCgen1_PS1.MegaCamGen1_from_PS1 import (
    G_MP9401,
    I_MP9701,
    R_MP9601,
    U_MP9301,
    Z_MP9801,
)

# mixed
from .MCgen2_PS1 import mixed as mixed_MC2PS1
from .MCgen2_PS1 import mixed as mixed_PS1MC2
from .MCgen2_PS1.MegaCamGen2_from_PS1 import I_MP9702
from .MCgen3_PS1 import mixed as mixed_MC3PS1
from .MCgen3_PS1 import mixed as mixed_PS1MC3
from .MCgen3_PS1.MegaCamGen3_from_PS1 import (
    G_MP9402,
    GRI_MP9605,
    I_MP9703,
    R_MP9602,
    U_MP9302,
    Z_MP9901,
)

#############################################################################
# CODE
#############################################################################


def _b1mb2(ps: Any, band1: str, band2: str, **kw: Any) -> Any:
    """B1 - B2."""
    b1 = globals().get(band1)
    b2 = globals().get(band2)
    return b1(ps, **kw) - b2(ps, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def UmG(
    ps: Table,
    u_band: Literal["9301", "9302"] = "9302",
    g_band: Literal["9401", "9402"] = "9402",
    **kw
) -> units.mag:
    """U-G color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for u(g)_band functions

    u_band: {'9301', '9302'}
        R band to use
        (default '9302')

    g_band : {'9401', '9402'}
        G band to use
        (default '9402')

    kwargs
        passes to U & G-band functions

    Returns
    -------
    U-G color

    """
    return _b1mb2(ps, "U_MP" + u_band, "G_MP" + g_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def UmR(
    ps: Table,
    u_band: Literal["9301", "9302"] = "9302",
    r_band: Literal["9601", "9602"] = "9602",
    **kw
) -> units.mag:
    """U-R color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for u(r)_band functions

    u_band: {'9301', '9302'}
        R band to use
        (default '9302')

    r_band: {'9601', '9602'}
        R band to use
        (default '9602')

    kwargs
        passes to U & R-band functions

    Returns
    -------
    U-R color

    """
    return _b1mb2(ps, "U_MP" + u_band, "R_MP" + r_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def UmI(
    ps: Table,
    u_band: Literal["9301", "9302"] = "9302",
    i_band: Literal["9701", "9702", "9703"] = "9703",
    **kw
) -> units.mag:
    """U-I color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for u(i)_band functions

    u_band: {'9301', '9302'}
        U band to use
        (default '9302')

    i_band: {'9701', '9702', '9703'}
        I band to use
        (default '9703')

    kwargs
        passes to U & I-band functions

    Returns
    -------
    U-I color

    """
    return _b1mb2(ps, "U_MP" + u_band, "I_MP" + i_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def UmZ(
    ps: Table,
    u_band: Literal["9301", "9302"] = "9302",
    z_band: Literal["9801", "9901"] = "9901",
    **kw
) -> units.mag:
    """U-Z color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for u(z)_band functions

    u_band: {'9301', '9302'}
        U band to use
        (default '9302')

    z_band: {'9801', '9901'}
        Z band to use
        (default '9901')

    kwargs
        passes to U & Z-band functions

    Returns
    -------
    U-Z color

    """
    return _b1mb2(ps, "U_MP" + u_band, "Z_MP" + z_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def GmR(
    ps: Table,
    g_band: Literal["9401", "9402"] = "9402",
    r_band: Literal["9601", "9602"] = "9602",
    **kw
) -> units.mag:
    """G-R color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for g(r)_band functions

    g_band : {'9401', '9402'}
        G band to use
        (default '9402')

    r_band: {'9601', '9602'}
        R band to use
        (default '9602')

    kwargs
        passes to G & R-band functions

    Returns
    -------
    G-R color

    """
    return _b1mb2(ps, "G_MP" + g_band, "R_MP" + r_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def GmI(
    ps: Table,
    g_band: Literal["9401", "9402"] = "9402",
    i_band: Literal["9701", "9702", "9703"] = "9703",
    **kw
) -> units.mag:
    """G-I color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for g(i)_band functions

    g_band : {'9401', '9402'}
        G band to use
        (default '9402')

    i_band: {'9701', '9702', '9703'}
        I band to use
        (default '9703')

    kwargs
        passes to G & I-band functions

    Returns
    -------
    G-I color

    """
    return _b1mb2(ps, "G_MP" + g_band, "I_MP" + i_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def GmZ(
    ps: Table,
    g_band: Literal["9401", "9402"] = "9402",
    z_band: Literal["9801", "9901"] = "9901",
    **kw
) -> units.mag:
    """G-Z color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for g(z)_band functions

    g_band : {'9401', '9402'}
        G band to use
        (default '9402')

    z_band: {'9801', '9901'}
        Z band to use
        (default '9901')

    kwargs
        passes to G & Z-band functions

    Returns
    -------
    G-Z color

    """
    return _b1mb2(ps, "G_MP" + g_band, "Z_MP" + z_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def RmI(
    ps: Table,
    r_band: Literal["9601", "9602"] = "9602",
    i_band: Literal["9701", "9702", "9703"] = "9703",
    **kw
) -> units.mag:
    """R-I color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for r(i)_band functions

    r_band: {'9601', '9602'}
        R band to use
        (default '9602')

    i_band: {'9701', '9702', '9703'}
        I band to use
        (default '9703')

    kwargs
        passes to R & I-band functions

    Returns
    -------
    R-I color

    """
    return _b1mb2(ps, "R_MP" + r_band, "I_MP" + i_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def RmZ(
    ps: Table,
    r_band: Literal["9601", "9602"] = "9602",
    z_band: Literal["9801", "9901"] = "9901",
    **kw
) -> units.mag:
    """R-Z color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for r(z)_band functions

    r_band: {'9601', '9602'}
        R band to use
        (default '9602')

    z_band: {'9801', '9901'}
        Z band to use
        (default '9901')

    kwargs
        passes to R & Z-band functions

    Returns
    -------
    R-Z color

    """
    return _b1mb2(ps, "R_MP" + r_band, "Z_MP" + z_band, **kw)


# /def

# -------------------------------------------------------------------


@units.quantity_io()
def ImZ(
    ps: Table,
    i_band: Literal["9701", "9702", "9703"] = "9703",
    z_band: Literal["9801", "9901"] = "9901",
    **kw
) -> units.mag:
    """G-I color.

    Parameters
    ----------
    ps : astropy.table.Table
        need arguments for i(z)_band functions

    i_band: {'9701', '9702', '9703'}
        I band to use
        (default '9703')

    z_band: {'9801', '9901'}
        Z band to use
        (default '9901')

    kwargs
        passes to I & Z-band functions

    Returns
    -------
    I-Z color

    """
    return _b1mb2(ps, "I_MP" + i_band, "Z_MP" + z_band, **kw)


# /def

# -------------------------------------------------------------------

#############################################################################
# END

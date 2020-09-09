# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Coordinate Systems Utilities."""

__author__ = "Nathaniel Starkman"
__credits__ = ["Astropy", "Jo Bovy", "Jeremy Webb"]

__all__ = [
    # projection functions
    "catalog_match_sky_track",
    "catalog_match_sky_inverse_track",
    "make_catalog_match_sky_track_functions",
    # functions
    "catalog_match_track",
    "catalog_match_inverse_track",
    "make_catalog_match_track_functions",
]


##############################################################################
# IMPORTS

import functools

import astropy.units as u
import numpy as np
from astropy.coordinates import match_coordinates_3d, match_coordinates_sky

##############################################################################
# PARAMETERS

dmls = u.dimensionless_unscaled


##############################################################################
# CODE
##############################################################################

#####################################################################
# projected version


def catalog_match_sky_track(
    coords, catalog, affine_param, adj_sep_sgn=True,
):
    """Map Astropy SkyCoord to orbit catalog.

    Does a ``match_coordinates_sky`` using orbit as catalog
    to get the time index and 2d distance.

    .. todo::

        update position_angle calculation to work when `catalog` or `coords`
        is a BaseCoordinateFrame, not just a SkyCoord.
        This would allow the interpolation functions to not construct
        a SkyCoord every time.

    Parameters
    ----------
    coords : SkyCoord
    catalog : SkyCoord
    affine_param : Sequence
        same length as catalog
    adj_sep_sgn : bool, optional
        whether the separation is the absolute value and the sign (+/-)
        needs to be inferred using the orbit's direction. Default is True
        to match `method` default, "closest".

    Returns
    -------
    afn : Quantity
        the affine parameter
    sep2d : Angle
    distance : Distance
    _PA : Angle
        position angle from orbit to original coordinate
        hidden variable

    Notes
    -----
    To get the sign of the separation (+/-) uses the direction
    of motion of the orbit, with right/left positive/negative
    This is done by approximating the orbit vector as the
    direction from the closest point (from the time index)
    to the next point. Rather than fancy trig, we instead use
    the difference in position angles (PA) between the point
    and closest orbit point (PA1) compared to the PA between
    point and the next orbit point (PA2).

        sep = sign(PA2 - PA1) |sep|

    """
    idx, sep2d, _, = match_coordinates_sky(coords, catalog)
    afn = affine_param[idx]

    if adj_sep_sgn:  # need to correct the sign of sep2d
        pa1 = catalog[idx].position_angle(coords) * dmls
        pa2 = catalog[idx + 1].position_angle(coords) * dmls
        sep2d *= np.sign(pa2 - pa1)

    return afn, sep2d, coords.distance, pa1


# /def

# -------------------------------------------------------------------


def catalog_match_sky_inverse_track(
    coords, catalog, affine_param,
):
    """Reverse mapping of above.

    Parameters
    ----------
    coords : orbitskyoffset_frame
        time, sep, distance
    catalog : SkyCoord
    affine_param : Sequence
        same length as catalog
    afn_name : str
        :deprecated:
        name of the affine param as attribute in `coords`

    Returns
    -------
    lon, lat : Quantity
    distance : Quantity

    """
    # TODO checks that `time` in affine_param range
    idx = np.argmin(  # get index into non-abs value
        np.abs(  # so get actual minimum, not most negative
            np.subtract(
                affine_param[:, None],
                coords.afn
                # getattr(coords, afn_name)
            )  # for broadcasting
        ),
        axis=0,  # column-wise, b/c broadcasted
    )

    orbit_pos = catalog[idx]

    # need to know offset direction. Most should have _PA
    pa = coords.data._PA  # TODO better access
    # but for anything that doesn't, need to find
    isnan = np.isnan(pa)
    # approximate by orthogonal to tangent at orbit_pos
    # do by average mean{PA(0, 1), PA(-1, 0)}
    # TODO more accurate PA, do this by orthogonal to on-sky tangent vec.
    pa[isnan] = (
        catalog[idx[isnan] - 1].position_angle(catalog[idx[isnan] + 1])
        * dmls  # for convert quantity list to list quantity
    ) + 90 * u.deg  # rotate to be orthogonal

    # Now offset by `sep` in direction `pa`
    out = orbit_pos.directional_offset_by(
        pa, np.abs(coords.sep)  # need abs() b/c `adj_sep_sgn`
    ).represent_as("spherical")

    return out.lon, out.lat, coords.distance


# /def


# -------------------------------------------------------------------


def make_catalog_match_sky_track_functions(
    catalog, affine_param, afn_name: str, adj_sep_sgn: bool = True,
):
    """Make Catalog Match Track Functions.

    Parameters
    ----------
    catalog:SkyCoord
    affine_param: Sequence
    afn_name: str
    adj_sep_sgn : bool, optional
        whether the separation is the absolute value and the sign (+/-)
        needs to be inferred using the orbit's direction. Default is True
        to match `method` default, "closest".

    Returns
    -------
    track_fn : Callable
        functools partial of catalog_match_sky_track with specified
        `catalog`, `affine_param`, `adj_sep_sgn`

    inverse_track_fn : Callable
        functools partial of catalog_match_sky_inverse_track with
        specified `catalog`, `affine_param`, `afn_name`

    """
    track_fn = functools.partial(
        catalog_match_sky_track,
        catalog=catalog,
        affine_param=affine_param,
        adj_sep_sgn=adj_sep_sgn,
    )

    inverse_track_fn = functools.partial(
        catalog_match_sky_inverse_track,
        catalog=catalog,
        affine_param=affine_param,
    )

    return track_fn, inverse_track_fn


# /def


##############################################################################
# 3D version


def catalog_match_track(
    coords, catalog, affine_param, adj_sep_sgn=True,
):
    """Map Astropy SkyCoord to orbit catalog.

    Does a ``match_coordinates_3d`` using orbit as catalog
    to get the time index and 3d distance.

    Returns
    -------
    afn : Quantity
        the affine parameter
    x, y : Quantity
        distance in tangent plane

    """
    idx, _, sep3d, = match_coordinates_3d(coords, catalog)
    afn = affine_param[idx]

    # now project onto plane orthogonal to curve
    # TODO more rigorous tangent vector method
    cart = catalog.cartesian
    tvec = cart[idx + 1] - cart[idx - 1]  # tangent vector
    that = tvec / np.linalg.norm(tvec)

    # define vectors along curve from central point
    delta2 = cart[idx + 1] - cart[idx]
    delta1 = cart[idx - 1] - cart[idx]

    fac = np.divide(np.dot(that, delta1), np.dot(that, delta2))

    xvec = delta1 - fac * delta2
    xhat = xvec / np.linalg.norm(xvec)

    yvec = np.cross(that, xvec)
    yhat = yvec / np.linalg.norm(yvec)

    # basic projection
    # https://en.wikipedia.org/wiki/Vector_projection
    x = np.dot(sep3d, xhat)
    y = np.dot(sep3d, yhat)
    # and error in projection
    d_afn = np.dot(sep3d, that)

    return afn, x, y, d_afn


# /def

# -------------------------------------------------------------------


def catalog_match_inverse_track(
    coords, catalog, affine_param,
):
    """Reverse mapping of above.

    .. todo::

        not need to do PA by having a hidden variable in `coords`
        that tracks the PA when converted to this frame.

    Parameters
    ----------
    coords : orbitoffset_frame
        time, sep, distance

    """
    # TODO checks that `time` in affine_param range
    idx = np.argmin(  # get index into non-abs value
        np.abs(  # so get actual minimum, not most negative
            np.subtract(
                affine_param[:, None],
                coords.afn
                # getattr(coords, afn_name)
            )  # for broadcasting
        ),
        axis=0,  # column-wise, b/c broadcasted
    )

    orbit_pos = catalog[idx].cartesian

    # now get the offset by x and y.
    # need to get the xhat and yhat direction
    # TODO more rigorous tangent vector method
    cart = catalog.cartesian
    tvec = cart[idx + 1] - cart[idx - 1]  # tangent vector
    that = tvec / np.linalg.norm(tvec)

    # define vectors along curve from central point
    delta2 = cart[idx + 1] - cart[idx]
    delta1 = cart[idx - 1] - cart[idx]

    fac = np.divide(np.dot(that, delta1), np.dot(that, delta2))

    xvec = delta1 - fac * delta2
    xhat = xvec / np.linalg.norm(xvec)

    yvec = np.cross(that, xvec)
    yhat = yvec / np.linalg.norm(yvec)

    offset = (
        (coords.x * xhat)
        + (coords.y * yhat)
        + (coords._d_afn * that)  # fixes error in projection
    )

    # get total position
    pos = orbit_pos + offset

    return pos.x, pos.y, pos.z


# /def

# -------------------------------------------------------------------


def make_catalog_match_track_functions(
    catalog, affine_param, afn_name: str, adj_sep_sgn: bool = True,
):
    """Make Catalog Match Track Functions.

    Parameters
    ----------
    catalog:SkyCoord
    affine_param: Sequence
    afn_name: str
    adj_sep_sgn : bool, optional
        whether the separation is the absolute value and the sign (+/-)
        needs to be inferred using the orbit's direction. Default is True
        to match `method` default, "closest".

    Returns
    -------
    track_fn : Callable
        functools partial of catalog_match_track with specified
        `catalog`, `affine_param`, `adj_sep_sgn`

    inverse_track_fn : Callable
        functools partial of catalog_match_inverse_track with
        specified `catalog`, `affine_param`, `afn_name`

    """
    track_fn = functools.partial(
        catalog_match_track,
        catalog=catalog,
        affine_param=affine_param,
        adj_sep_sgn=adj_sep_sgn,
    )

    inverse_track_fn = functools.partial(
        catalog_match_inverse_track,
        catalog=catalog,
        affine_param=affine_param,
    )

    return track_fn, inverse_track_fn


# /def

##############################################################################
# END

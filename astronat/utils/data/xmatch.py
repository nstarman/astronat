# -*- coding: utf-8 -*-

"""Catalog x-matching.

Built from the extensive x-matching program in Jo Bovy's ``gaia_tools``,
but reconfigured to use the Astropy framework as much as possible.

.. todo::

    - make a registry class-singleton for adding conversions to Table.

        - numpy recarray in some manner (?)
        - astropy NDData
        - pandas array
        - astroquery TableList
        - my (Q)TableList and TablesList

    - astroquery-style constraints classes for more constraint options.
      like color matching, etc.

Notes
-----
Jo Bovy's ``gaia_tools`` includes an excellent CDS x-match interface.
Astropy has basic cross-matching, documented
`here <https://docs.astropy.org/en/stable/coordinates/matchsep.html>_`
astroML has a very basic module "crossmatch"

for HEALPix, see ligo.skymap.postprocess.crossmatch

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["Jo Bovy", "Henry Leung"]


__all__ = [
    # imported functions
    "indices_xmatch_fields",
    "xmatch_fields",
    "non_xmatched",
    # functions
    "indices_xmatch_coords",
    "coord_has_duplicates",
    "xmatch_coords",
    # "indices_xfind_coords",  # TODO
    # "xfind_coords",  # TODO
    "XMatch",
    # conversions
    "xmatch_decorator",
]


##############################################################################
# IMPORTS

# BUILT-IN

import typing as T


# THIRD PARTY

import astropy.coordinates as coord
from astropy.coordinates import (
    SkyCoord,
    BaseCoordinateFrame as BCFrame,
)
from astropy.table import Table
import astropy.units as u

import numpy as np

from utilipy.utils import inspect, functools
from utilipy.data_utils.xmatch import (
    indices_xmatch_fields,
    xmatch_fields,
    non_xmatched,
)


# PROJECT-SPECIFIC

from .core import xmatch_data_graph as xdg


##############################################################################
# CODE
##############################################################################


def xmatch_decorator(
    function: T.Callable = None,
    *,
    _doc_style: str = "numpy",
    _doc_fmt: T.Dict[str, T.Any] = {},
    **catalogs,
):
    """Decorator for correctly dealing with Table inputs.

    Manage Table inputs to x-match functions. Convert. Clean up. etc.

    .. todo::

        - actually write function
        - natively support SkyCoord, BaseCoordinateFrame, (Q)Table.
          (extend via __astropy_table__ method)
        - make a registry class-singleton for adding conversions.

            - numpy recarray in some manner (?)
            - astropy NDData
            - pandas array
            - astroquery TableList
            - my (Q)TableList and TablesList

        - (here ?) take obstime & epoch into account
        - a "return" argument in "catalogs" tp return the munged data

    Parameters
    ----------
    function : Callable or None, optional
        the function to be decoratored
        if None, then returns decorator to apply.
    **catalogs: dict
        catalog information, where keyword is the argument name in `function`
        the values are themselves dictionaries, containing the information
        regarding the specified catalog.

        Dictionary Options:

            - dtype: class, optional
                the output (sub)class for the munged catalog(s).
                Coordinate classes are different than photometric classes, etc
                and this decorator needs to know what it's dealing with.

                For coordinates,
                :class:`~astropy.coordinates.BaseCoordinateFrame`
                is treated as a sub-class hint and the frame is undisturbed
                (even if its actually a SkyCoord).
                Unless a specific frame is needed, this is the preferred input.

            - return_dtype : bool, optional
                a "return" argument in "catalogs" to return the munged data
                TODO!!!

    Returns
    -------
    wrapper : Callable
        wrapper for `function` that manage input catalog tables.
        includes the original function in a method `.__wrapped__`

    Other Parameters
    ----------------
    _doc_style: str or formatter, optional
        `function` docstring style. Parameter to `utilipy.wraps`.
    _doc_fmt: dict, optional
        `function` docstring format arguments. Parameter to `utilipy.wraps`.

    """
    # ------------------------------
    # need to make sure **catalogs contains correct information

    # allowed_keys = {"dtype"}  # TODO detect non-allowed

    for name, info in catalogs.items():
        if "dtype" not in info:
            raise ValueError(f"Specify an `dtype` in catalog {name}")

    # ------------------------------
    if function is None:  # allowing for optional arguments
        return functools.partial(
            xmatch_decorator,
            _doc_style=_doc_style,
            _doc_fmt=_doc_fmt,
            **catalogs,
        )

    sig = inspect.fuller_signature(function)

    @functools.wraps(function, _doc_style=_doc_style, _doc_fmt=_doc_fmt)
    def wrapper(*args, _skip_decorator=False, **kwargs):
        """Wrapper docstring.

        Other Parameters
        ----------------
        _skip_decorator : bool, optional
            Whether to skip the decorator.
            default {_skip_decorator}

        Notes
        -----
        The wrapper

        """
        # whether to skip decorator or keep going
        if _skip_decorator:
            return function(*args, **kwargs)
        # else:

        ba = sig.bind_partial_with_defaults(*args, **kwargs)

        obstime = ba.arguments.get("obstime", None)
        # print("obstime: ", obstime)

        # Catalog munging
        for name, info in catalogs.items():
            catalog = ba.arguments[name]

            # get obstime if not provided or determined
            if obstime is None and hasattr(catalog, "obstime"):
                obstime = getattr(catalog, "obstime")
                # print(obstime)  # TODO report this somehow

            # Step 1, convert to correct output type
            # TODO use registry for this
            cat_dtype = info["dtype"]  # output type

            # TODO support passing multiple types
            if cat_dtype is BCFrame:

                if isinstance(  # note, works if subclass
                    catalog, (BCFrame, SkyCoord)
                ):
                    pass
                elif isinstance(catalog, Table):
                    catalog = xdg.get_transform(Table, SkyCoord)(catalog)

                # elif:

                else:
                    raise TypeError(
                        (
                            f"Unknown conversion {type(catalog)} "
                            "-> BaseCoordinateFrame."
                        )
                    )

            elif isinstance(cat_dtype, SkyCoord):
                if obstime is None and hasattr(catalog, "obstime"):
                    obstime = getattr(catalog, "obstime")

                if isinstance(catalog, SkyCoord):  # note, works if subclass
                    pass
                elif isinstance(catalog, BCFrame):
                    catalog = xdg.get_transform(BCFrame, SkyCoord)(catalog)

                elif isinstance(catalog, Table):
                    catalog = xdg.get_transform(Table, SkyCoord)(catalog)

            # elif isinstance(cat_dtype, BCFrame):  # a specific Frame class

            else:
                raise Exception("TODO")

            # Step 2, evolve observations by `obstime`
            try:
                new_catalog = SkyCoord.apply_space_motion(catalog, obstime)
            except Exception as e:  # can't evolve
                # need to check that this is not a problem
                if hasattr(catalog, "obstime"):
                    if catalog.obstime != obstime:  # Uh oh!, it is
                        raise Exception(
                            (
                                f"epoch mismatch: catalog {name}'s "
                                "coordinates cannot be transformed "
                                f"to match the epoch {obstime} "
                                f"(Exception {e})."
                            )
                        )
                    else:
                        pass  # the epochs match

                # TODO, do as warning instead
                # print(f"Not adjusting catalog {name} ({e})")
                new_catalog = catalog

            # reassign
            ba.arguments[name] = new_catalog

        # /def

        # ----------------------------------

        return_ = function(*ba.args, **ba.kwargs)

        # replace_catalogs = {}  # TODO implement when have "return" working
        # for name, info in catalogs.items():
        #     if info.get("return", False):  # get munged data
        #         replace_catalogs[name] = ba.arguments[name]

        # if replace_catalogs:  # it's not empty
        #     return return_, replace_catalogs

        return return_

    # /def

    return wrapper


# /def


##############################################################################
# Coordinate Cross-Matching
# where everything is nice


@xmatch_decorator(
    catalog={"dtype": BCFrame}, other={"dtype": BCFrame},
)
@u.quantity_input(maxdist=("angle", "length"))
def indices_xmatch_coords(
    catalog, other, maxdist=1 * u.arcsec, obstime=None, nthneighbor: int = 1,
):
    """Basic 2-catalog x-match. Returns match indices.

    see https://docs.astropy.org/en/stable/coordinates/matchsep.html

    Parameters
    ----------
    catalog, other : SkyCoord or BaseCoordinateFrame
        `catalog` is the "catalogcoord", `other` are the "matchcoord".
        Note that this is in the opposite order as Astropy.
    maxdist : `~astropy.units.Angle` or `~astropy.coordinates.Distance`, optional
        The maximum separation to be considered a match.
        If Angle, does an on-sky x-match using
        :func:`~astropy.coordinates.match_coordinates_sky`.
        If Distance, does a 3d x-match using
        :func:`~astropy.coordinates.match_coordinates_3d`.
    obstime : Time, optional
        If provided, the "epoch" at which to x-match the coords.
        If None (default), will use the obstime of the first catalog to have
        an obstime. An absence of obstime in a catalog means the coordinates
        are *right now*, if this is not the case, ensure the catalog has this
        information.

    Returns
    -------
    catalog_idx : integer array
        indices into `catalog` for the x-match.
    info : dict
        Useful information.

            - sep2d : on-sky separation (Angle)
            - dist3d : 3D distance (Quantity)

    Other Parameters
    ----------------
    nthneighbor : int
        The nthneighbor to use in ``match_coordinates_``.
        TODO, rename to "_nth" but "quantity_input" can't handle underscores.

    See Also
    --------
    :func:`~astropy.coordinates.match_coordinates_sky`
    :func:`~astropy.coordinates.match_coordinates_3d`

    """
    if u.get_physical_type(maxdist.unit) == "angle":
        idx, sep2d, dist3d = coord.match_coordinates_sky(
            matchcoord=other, catalogcoord=catalog, nthneighbor=nthneighbor,
        )
        sep = sep2d  # separation constraints on this
    elif u.get_physical_type(maxdist.unit) == "length":
        idx, sep2d, dist3d = coord.match_coordinates_3d(
            matchcoord=other, catalogcoord=catalog, nthneighbor=nthneighbor,
        )
        sep = dist3d  # separation constraints on this

    # separation constraints
    midx = np.where(sep < maxdist)[0]

    if len(midx) == 0:  # no matches
        return False, False, {}
    elif len(midx) == 1:  # correct for case of only 1 match
        midx = midx[0]
        sel = ...
    else:
        sel = midx

    catalog_idx = idx[sel]
    other_idx = midx

    info = {"sep2d": sep2d[sel], "dist3d": dist3d[sel]}

    return catalog_idx, other_idx, info


# /def

# -------------------------------------------------------------------


def coord_has_duplicates(catalog, tol=0.01 * u.arcsec, dup_lim=0):
    """Check for duplicates in the catalog to within `tol`.

    Since catalog matching can be an expensive operation, this also
    returns the indices and info.

    Parameters
    ----------
    catalog : SkyCoord or BaseCoordinateFrame
    tol : `~astropy.units.Angle` or `~astropy.coordinates.Distance`, optional
        The maximum separation to be considered a match.
        If Angle, does an on-sky x-match using
        :func:`~astropy.coordinates.match_coordinates_sky`.
        If Distance, does a 3d x-match using
        :func:`~astropy.coordinates.match_coordinates_3d`.
    dup_lim : int, optional
        number of duplicates allowed before failing
        (see Return for how this is handled)

    Return
    ------
    result : bool
        * False = passed the `dup_lim`, there are NO duplicates
        * True = failed the `dup_lim`, there ARE duplicates
    cat_idx : integer array or None
        If passed, indices into catalogs for the x-match.
        If failed, then True, because no duplicates found.
    info : dict or None
        If passed, info dict.
        If failed, then None, because no duplicates found.

    See Also
    --------
    :func:`~astronat.utils.data.xmatch.indices_xmatch_coords`
    :func:`~astronat.utils.data.xmatch.xmatch_coords`

    """
    cat_idx, _, info = indices_xmatch_coords(
        catalog, catalog, maxdist=tol, _nthneighbor=2
    )

    raise Exception(
        (
            "TODO: I think I need to reverse this to "
            "get the indices of the non-duplicates"
        )
    )

    if np.sum(cat_idx) > dup_lim:
        return True, cat_idx, info
    return False, True, None


# /def


# -------------------------------------------------------------------


@u.quantity_input(maxdist=("angle", "length"))
def xmatch_coords(
    catalog, other, maxdist: u.Quantity = 1 * u.arcsec, obstime=None
):
    """Basic 2-catalog x-match.

    see https://docs.astropy.org/en/stable/coordinates/matchsep.html

    Parameters
    ----------
    catalog, other : SkyCoord or BaseCoordinateFrame
        `catalog` is the "catalogcoord", `other` are the "matchcoord".
        Note that this is in the opposite order as Astropy,.
    maxdist : `~astropy.units.Angle` or `~astropy.coordinates.Distance`, optional
        The maximum separation to be considered a match.
        If Angle, does an on-sky x-match using
        :func:`~astropy.coordinates.match_coordinates_sky`.
        If Distance, does a 3d x-match using
        :func:`~astropy.coordinates.match_coordinates_3d`.
    obstime : Time, optional
        If provided, the "epoch" at which to x-match the coords.
        If None (default), will use the obstime of the first catalog to have
        an obstime. An absence of obstime in a catalog means the coordinates
        are *right now*, if this is not the case, ensure the catalog has this
        information.

    Returns
    -------
    catalog_matches, other_matches : catalog input types
        the x-matches
    info : dict
        Useful information.

            - catalog_idx : indices into `catalog` for the x-match.
            - m1 : indices into `cat2` for the x-match.
            - sep2d : on-sky separation (Angle)
            - dist3d : 3D distance (Quantity)

    See Also
    --------
    :func:`~astronat.utils.data.xmatch.indices_xmatch_coords`

    """
    catalog_idx, other_idx, info = indices_xmatch_coords(
        catalog,
        other,
        maxdist=maxdist,
        obstime=obstime,
        _skip_decorator=False,
    )

    info.update({"catalog_idx": catalog_idx, "other_idx": other_idx})

    # TODO more succinct, try-except?
    if not np.shape(catalog):  # it's a scalar and cannot be indexed
        catalog_matches = np.array([catalog])[
            catalog_idx
        ]  # see if the index works on 1-elem list
    else:
        catalog_matches = catalog[catalog_idx]

    if not np.shape(other):  # it's a scalar and cannot be indexed
        other_matches = np.array([other])[
            other_idx
        ]  # see if the index works on 1-elem list
    else:
        other_matches = other[other_idx]

    return catalog_matches, other_matches, info


# /def


# -------------------------------------------------------------------


@xmatch_decorator(
    coords1={"dtype": BCFrame}, coords2={"dtype": BCFrame},
)
@u.quantity_input(maxdist=("angle", "length"))
def indices_xfind_coords(
    coords1, coords2, maxdist=1 * u.arcsec, obstime=None,
):
    """Indices of X-Find coords.

    Returns
    -------
    idx1 : integer array
    idx2 : integer array
        indices into `other` for all coordinates which have a "match"
        in `catalog`.
    info : dict
        Useful information  # TODO

    See Also
    --------
    :func:`~astropy.coordinates.search_around_sky`
    :func:`~astropy.coordinates.search_around_3d`

    """
    if u.get_physical_type(maxdist.unit) == "angle":
        idx1, idx2, sep2d, dist3d = coord.search_around_sky(
            coords1, coords2, seplimit=maxdist,
        )
    elif u.get_physical_type(maxdist.unit) == "length":
        idx1, idx2, sep2d, dist3d = coord.search_around_3d(
            coords1, coords2, distlimit=maxdist,
        )

    info = {"sep2d": sep2d, "dist3d": dist3d}

    return idx1, idx2, info


# /def


# -------------------------------------------------------------------


@u.quantity_input(maxdist=("angle", "length"))
def xfind_coords(
    coords1, coords2, maxdist: u.Quantity = 1 * u.arcsec, obstime=None,
):
    """X-Find Coords.

    Returns
    -------
    coords1_idx : integer array
        indices into `coords1` for all coordinates within maxdist
        off the coords in `coords2`.
    coords2_idx : integer array
        indices into `coords2` for all coordinates which have a "match"
        in `coords1`.
    info : dict
        Useful information  # TODO

    See Also
    --------
    :func:`~astronat.utils.data.xmatch.indices_xfind_coords`

    """
    coords1_idx, coords2_idx, info = indices_xfind_coords(
        coords1, coords2, maxdist=maxdist, obstime=obstime
    )

    info.update({"coords1_idx": coords1_idx, "coords2_idx": coords2_idx})

    coords1_matches = coords1[coords1_idx]
    coords2_matches = coords2[coords2_idx]

    return coords1_matches, coords2_matches, info


# /def


# -------------------------------------------------------------------


@u.quantity_input(maxdist=("angle", "length"))
def xmatch_many_coords(
    *catalogs, maxdist: u.Quantity = 1 * u.arcsec, obstime=None
):
    """Cross-match all the catalogs.

    see https://docs.astropy.org/en/stable/coordinates/matchsep.html

    Parameters
    ----------
    catalogs : SkyCoord or BaseCoordinateFrame
        `catalogs` is the "catalogcoord", `others` are the "matchcoord".
        Note that this is in the opposite order as Astropy,.
    maxdist : `~astropy.units.Angle` or `~astropy.coordinates.Distance`, optional
        The maximum separation to be considered a match.
        If Angle, does an on-sky x-match using
        :func:`~astropy.coordinates.match_coordinates_sky`.
        If Distance, does a 3d x-match using
        :func:`~astropy.coordinates.match_coordinates_3d`.
    obstime : Time, optional
        If provided, the "epoch" at which to x-match the coords.
        If None (default), will use the obstime of the first catalog to have
        an obstime. An absence of obstime in a catalog means the coordinates
        are *right now*, if this is not the case, ensure the catalog has this
        information.

    Returns
    -------
    catalog_matches, *others_matches : catalog input types
        the x-matches

    Notes
    -----
    .. todo::
        return an info dict, like xmatch_coords

    See Also
    --------
    :func:`~astronat.utils.data.xmatch.indices_xfind_coords`

    """
    raise ValueError("TODO")

    # IDEA 2
    # starting from last catalog (assumed smallest?), use the search_around_*
    # function to get all candidate matches between last and penultimate.
    # throw out all stars in ultimate and penultimate which do not have any
    # candidate matches.
    # This is the new ultimate catalog
    # repeat above with 3rd last (now penultimate) and new ultimate catalog
    # work back until reach 1st catalog

    # list of indices for each catalog
    idxs = [] * len(catalogs)
    # the indices that the catalog one index up cares about
    # needed to detect whan a down-index match orphans an up-index coord.
    down_idx = [] * len(catalogs)

    # start with an all-true array
    catalog_idx = True
    down_idx[-1] = True

    for i, j in [(-2, -1), (-3, -2), (-4, -3)]:  # TODO go all the way down

        catalog_idx, other_idx, _ = indices_xfind_coords(
            catalogs[i],
            catalogs[j][catalog_idx],
            maxdist=maxdist,
            obstime=obstime,
        )

        # store the indices for matched values in "matchcoord"
        idxs[j] = np.unique(other_idx)
        # TODO detect if any up "up-index" catalogs have orphaned values
        down_idx[i] = catalog_idx

        catalog_idx = np.unique(catalog_idx)  # remove duplicates

    # catalog_idx, other_idx, _ = xmatch_coords(
    #     catalog,
    #     other,
    #     maxdist=maxdist,
    #     obstime=obstime,
    #     _skip_decorator=False,
    # )

    # return catalog_matches, other_matches


# /def


##############################################################################
# General Cross-Match
# allowing for more complicated x-matching
# TODO, allow arbitrary match fields


# @xmatch_decorator(
#     cat1={"dtype": BaseCoordinateFrame}, cat2={"dtype": BaseCoordinateFrame},
# )
# @static_citation_decorator(citation="https://github.com/jobovy/gaia_tools/")
@u.quantity_input(maxdist=("angle", "length"))
def XMatch(
    cat1,
    cat2,
    maxdist=1 * u.arcsec,
    obstime=None,
    swap: bool = False,
    match_fields: T.Optional[dict] = None,
    # **kwargs
):
    """Cross-match two catalogs.

    .. todo::

        - Correcting velocities for solar reflex motion ?
        - extra match fields, with constraints
        - find duplicates before matching ?
        - find duplicates while matching. It can happen
        - combine many matches simultaneously
            the difficulty is matches with metrics, like xmatch_coords.


    Parameters
    ----------
    cat1, cat2: Any
        the two catalogs to crossmatch
        must be convertible to SkyCoords (see transformation registry)
    maxdist : `~astropy.units.Angle` or `~astropy.coordinates.Distance`, optional
        The maximum separation to be considered a match.
        If Angle, does an on-sky x-match using
        :func:`~astropy.coordinates.match_coordinates_sky`.
        If Distance, does a 3d x-match using
        :func:`~astropy.coordinates.match_coordinates_3d`.
    swap : bool, optional
        If True reverse roles of `cat1` and `cat2`, (default False).
        Important when one of the catalogs has duplicates
    match_fields : str, optional
        data tag on which to additionally cross-match, default None.
        this also works for any discrete-valued data column.
        TODO, allow this to be a list of fields
    **kwargs
        for decorator?

    See Also
    --------
    :func:`~astronat.utils.data.xmatch.xmatch_coords`
    :func:`~astronat.utils.data.xmatch.indices_xmatch_coords`

    References
    ----------
    https://github.com/jobovy/gaia_tools/

    """
    if swap:  # not default
        c1 = cat2
        c2 = cat1
    else:
        c1 = cat1
        c2 = cat2

    if match_fields is None:
        c1_matches, c2_matches, info = xmatch_coords(
            c1, c2, maxdist, obstime=obstime
        )

    else:  # TODO implement more of this in the decorator somehow?

        # TODO general cross-match function
        # and a function for combining cross-matches

        # TODO speed up
        # can partly speed up by extracting the SkyCoord to be
        # used by indices_xmatch_coords, so it doesn't do that every
        # time

        # The code in `# <-- # -->` is modified from `gaia_tools`
        if isinstance(match_fields, str):
            field = match_fields

            # <--
            # does this field exist?  # TODO not assume type(Table)
            if field not in cat1.colnames:
                raise KeyError(f"'{field}' not in catalog cat1")
            if field not in cat2.colnames:
                raise KeyError(f"'{field}' not in catalog cat2")

            # need to get the field
            values_c1 = c1[field]
            len_vc1 = len(values_c1)

            values_c2 = c2[field]
            len_vc2 = len(values_c2)

            uniques = np.unique(values_c1)

            idxc1 = []
            idxc2 = []
            info = {"sep2d": [], "dist3d": []}
            for unique in uniques:
                idx_1 = np.arange(len_vc1)[values_c1 == unique]
                idx_2 = np.arange(len_vc2)[values_c2 == unique]

                temp_idx1, temp_idx2, temp_info = indices_xmatch_coords(
                    c1[idx_1], c2[idx_2], maxdist=maxdist, obstime=obstime
                )

                idxc1.extend(idx_1[temp_idx1])
                idxc2.extend(idx_2[temp_idx2])
                info["sep2d"].extend(temp_info.get("sep2d", []))
                info["dist3d"].extend(temp_info.get("dist3d", []))

            # -->

            idxc1 = np.array(idxc1)
            idxc2 = np.array(idxc2)

            info.update({"m1": idxc1, "m2": idxc2})

            c1_matches = cat1[idxc1]
            c2_matches = cat2[idxc2]

        else:
            raise Exception("TODO")

    if swap:  # switch back
        cat1_matches = c2_matches
        cat2_matches = c1_matches
    else:
        cat1_matches = c1_matches
        cat2_matches = c2_matches

    return cat1_matches, cat2_matches, info


# /def


##############################################################################
# END

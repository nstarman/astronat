# -*- coding: utf-8 -*-

"""Orbit SkyOffset Frames.

Like `~astropy.coordinates.Galactocentric` needs ``galcen_coord`` and other
information to convert to/from heliocentric or geocentric coordinate systems,
the `OrbitPseudoFrame` needs information about the orbit (phase-space position
and integration time) and the potential in which the orbit was integrated.

The information name, and frame attribute data descriptor are listed:

    - origin : `~astropy.coordinates.attributes.CoordinateAttribute`
    - potential: `PotentialAttribute`
    - afn_bounds : `~astropy.coordinates.attributes.QuantityAttribute`

Notes
-----
This coordinate system should be used for visualization, probably not science
for two reasons: it is under active development, and its very difficult
to interpret what kinematics mean when the coordinate system itself is a
function of an affine parameter.

.. todo::

    - Look at FunctionTransformWithFiniteDifference for the kinematics

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["Jeremy Webb", "Jo Bovy"]


__all__ = [
    "OrbitSkyOffsetFrame",
]


##############################################################################
# IMPORTS


import typing as T

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord  # TODO not need
from astropy.coordinates import (
    concatenate,
    frame_transform_graph,
    match_coordinates_sky,
)
from astropy.coordinates.attributes import (
    CoordinateAttribute,
    QuantityAttribute,
)
from astropy.coordinates.baseframe import BaseCoordinateFrame
from astropy.coordinates.representation import (
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.coordinates.transformations import FunctionTransform
from scipy import interpolate
from typing_extensions import Literal

from ..attributes import PotentialAttribute
from ..representations import (
    OrbitSkyOffsetRepresentation,
    OrbitSkyOffsetUnitRepresentation,
)
from ..transformations import FunctionWithKwargTransform
from .utils import (
    catalog_match_sky_inverse_track as catalog_match_inverse_track,
)
from .utils import catalog_match_sky_track as catalog_match_track

# from astropy.coordinates.attributes import Attribute
# from astropy.coordinates.baseframe import RepresentationMapping

##############################################################################
# PARAMETERS

_orbitskyoffset_cache = {}
"""The cache in which OrbitSkyOffsets are stored."""


##############################################################################
# CODE
##############################################################################


def make_orbitskyoffset_cls(
    framecls,
    track_fn,
    inverse_track_fn=None,
    track_fn_kw={},
    inverse_track_fn_kw={},
):
    """Make Sky-Offset Class from an Orbit in a potential.

    Create a new class that is the orbit sky offset frame for a specific class
    of origin frame. If such a class has already been created for this frame,
    the same class will be returned.

    The new class will always have component names "afn", "sep", "distance".

    Parameters
    ----------
    framecls : subclass of `~astropy.coordinates.BaseCoordinateFrame`
        coordinate frame class. The class to create the SkyOffsetFrame.
    potential : Any
        The potential in which the track was integrated
    track_fn : Callable[[affine param array], coordinate array]
        Function mapping affine parameter to the coordinate array.
    inverse_track_fn : Callable[[coordinate array], affine param array]
        Function mapping coordinate array to affine parameter.
    track_fn_kw : dict
        keyword arguments into `track_fn`
    inverse_track_fn_kw : dict
        keyword arguments into `inverse_track_fn`

    Returns
    -------
    orbitskyoffsetframecls : class
        The class for the new orbit-skyoffset frame.

    Notes
    -----
    This function is necessary because Astropy's frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct orbit-skyoffset frame class.  This function
    generates just that class, as well as ensuring that only one example of
    such a class actually gets created in any given python session.

    For implementation details, see
    :mod:`~astropy/coordinates/builtin_frames/skyoffset`

    .. todo::

        - fix cache. it erroneously returns same class for everything
          see SkyOffsetFrame for reference. maybe need to set cache key
          as a hash of the framecls, origin, potential, and trackfunc?

    """
    # if framecls in _orbitskyoffset_cache:  # FIXME, all are the same
    #     return _orbitskyoffset_cache[framecls]

    # the class of a class object is the metaclass
    framemeta = framecls.__class__

    # -----------------------------------------------------

    class OrbitSkyOffsetMeta(framemeta):
        """Metaclass for Orbit Sky-Offsets.

        This metaclass renames the class to be "SkyOffset<framecls>" and also
        adjusts the frame specific representation info so that spherical names
        are always "lon" and "lat" (instead of e.g. "ra" and "dec").

        """

        def __new__(cls, name, bases, members):
            # Only 'origin' is needed here, to set the origin frame properly.
            members["origin"] = CoordinateAttribute(
                frame=framecls, default=None
            )

            # This has to be done because FrameMeta will set these attributes
            # to the defaults from BaseCoordinateFrame when it creates the base
            # OrbitSkyOffsetFrame class initially.
            members["_default_representation"] = OrbitSkyOffsetRepresentation
            members[
                "_default_differential"
            ] = framecls._default_differential  # TODO replace

            newname = name[:-5] if name.endswith("Frame") else name
            newname += framecls.__name__

            return super().__new__(cls, newname, bases, members)

        # /def

    # /class

    # We need this to handle the intermediate metaclass correctly, otherwise
    # we could just subclass OrbitSkyOffsetFrame.
    _OrbitSkyOffsetFramecls = OrbitSkyOffsetMeta(
        "OrbitSkyOffsetFrame",
        (OrbitSkyOffsetFrame, framecls),
        {"__doc__": OrbitSkyOffsetFrame.__doc__},
    )

    # -----------------------------------------------------
    # register frame transform graph to/from reference from/to OrbitSkyOffset

    @frame_transform_graph.transform(
        FunctionWithKwargTransform,
        framecls,
        _OrbitSkyOffsetFramecls,
        func_kwargs=track_fn_kw,  # the **kwargs
    )
    def reference_to_orbitskyoffset(
        reference_coord, orbitskyoffset_frame, **kwargs
    ):
        """Compute the transformation from reference to orbit frame.

        Notes
        -----
        unlike matrix transforms, the FunctionTransform method requires
        manually passing the frame attributes. This can be generalized using
        the ``frame_attribute`` property.

        """
        afn_name = kwargs.pop("afn_name", None)
        afn, sep2d, distance, _PA = track_fn(reference_coord, **kwargs)

        # now need to decide how to treat distances, if they are / not present
        # this might be a temporary solution, but rather than directly
        # implementing a _OrbitSkyOffsetFramecls, first create a representation
        # normal or unit (if no distances) and then create the frame from the
        # representation, specifying which representation type was used.
        rep = OrbitSkyOffsetRepresentation(
            afn=afn, sep=sep2d, distance=distance, _PA=_PA, afn_name=afn_name,
        )
        representation_type = OrbitSkyOffsetRepresentation

        if all(reference_coord.distance.to_value() == 1) and (
            reference_coord.distance.unit == u.dimensionless_unscaled
        ):
            rep = rep.represent_as(OrbitSkyOffsetUnitRepresentation)
            representation_type = OrbitSkyOffsetUnitRepresentation

        # /if

        return _OrbitSkyOffsetFramecls(
            rep,
            representation_type=representation_type,
            # **orbitskyoffset_frame.frame_attributes  # TODO
            origin=orbitskyoffset_frame.origin,
            potential=orbitskyoffset_frame.potential,
            afn_bound_tail=orbitskyoffset_frame.afn_bound_tail,  # FIXME, shape prob
            afn_bound_lead=orbitskyoffset_frame.afn_bound_lead,
        )

    # /def

    @frame_transform_graph.transform(
        FunctionWithKwargTransform,
        _OrbitSkyOffsetFramecls,
        framecls,
        func_kwargs=inverse_track_fn_kw,  # the **kwargs
    )
    def skyoffset_to_reference(
        orbitskyoffset_coord, reference_frame, **kwargs
    ):
        """Convert an sky offset frame coordinate to the reference frame."""
        lon, lat, distance = inverse_track_fn(orbitskyoffset_coord, **kwargs)
        rep = SphericalRepresentation(lon=lon, lat=lat, distance=distance)

        if not hasattr(orbitskyoffset_coord, "distance"):
            rep = rep.represent_as(UnitSphericalRepresentation)

        return framecls(
            rep, representation_type=reference_frame.representation_type
        )

    # /def

    # -----------------------------------------------------
    # register between OrbitSkyOffset frame transforms

    # TODO
    @frame_transform_graph.transform(
        FunctionTransform, _OrbitSkyOffsetFramecls, _OrbitSkyOffsetFramecls
    )
    def skyoffset_to_skyoffset(from_skyoffset_coord, to_skyoffset_frame):
        """Transform between two orbit-skyoffset frames.

        Parameters
        ----------
        from_skyoffset_coord
        to_skyoffset_frame

        Returns
        -------
        to_skyoffset_coord

        """
        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
        intermediate_from = from_skyoffset_coord.transform_to(
            from_skyoffset_coord.origin
        )
        intermediate_to = intermediate_from.transform_to(
            to_skyoffset_frame.origin
        )

        return intermediate_to.transform_to(to_skyoffset_frame)

    # /def

    # -----------------------------------------------------

    _orbitskyoffset_cache[framecls] = _OrbitSkyOffsetFramecls

    return _OrbitSkyOffsetFramecls


# /def


# -------------------------------------------------------------------


class OrbitSkyOffsetFrame(BaseCoordinateFrame):
    """Sky Offset Frame from an Orbit in a Potential.

    A frame which is relative to some specific position and oriented to match
    its frame.

    SkyOffsetFrames always have component names for spherical coordinates
    of ``lon``/``lat``, *not* the component names for the frame of ``origin``.

    This is useful for calculating offsets and dithers in the frame of the sky
    relative to an arbitrary position. Coordinates in this frame are both
    centered on the position specified by the ``origin`` coordinate, *and*
    they are oriented in the same manner as the ``origin`` frame.  E.g., if
    ``origin`` is `~astropy.coordinates.ICRS`, this object's ``lat`` will be
    pointed in the direction of Dec, while ``lon`` will point in the direction
    of RA.

    For more on skyoffset frames, see :ref:`astropy-skyoffset-frames`.

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data,
        or use the other keywords
    origin : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the origin of this frame. Note that this
        origin is used purely for on-sky location.  It can have a
        ``distance`` but it will not be used by this ``OrbitSkyOffsetFrame``.
    potential : `~galpy.potential.Potential` or list thereof.

    Notes
    -----
    ``OrbitSkyOffsetFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``OrbitSkyOffsetFrame``.
    Instead, distinct classes are created on-the-fly for whatever the frame
    class is of ``origin``.

    """

    origin = CoordinateAttribute(default=None, frame=None)
    potential = PotentialAttribute()
    afn_bound_tail = QuantityAttribute(default=-np.inf * u.Myr, unit=u.Myr)
    afn_bound_lead = QuantityAttribute(default=np.inf * u.Myr, unit=u.Myr)

    @property
    def afn_bounds(self):
        """Affine bounds (tail, lead)."""
        return u.Quantity([self.afn_bound_tail, self.afn_bound_lead])

    # /def

    def __new__(cls, *args, **kwargs):
        """OrbitSkyOffsetFrame."""
        # We don't want to call this method if we've already set up
        # an orbitskyoffset frame for this class.
        if not (
            issubclass(cls, OrbitSkyOffsetFrame)
            and cls is not OrbitSkyOffsetFrame
        ):
            # We get the origin argument, and handle it here.
            try:
                origin_frame = kwargs["origin"]
            except KeyError:
                raise TypeError(
                    "Can't initialize an OrbitSkyOffsetFrame "
                    "without origin= keyword."
                )

            try:
                track_fn = kwargs.pop("track_fn")
            except KeyError:
                raise TypeError(
                    "Can't initialize an OrbitSkyOffsetFrame "
                    "without origin= keyword."
                )

            inverse_track_fn = kwargs.pop("inverse_track_fn", None)

            track_fn_kw = kwargs.pop("track_fn_kw", {})
            inverse_track_fn_kw = kwargs.pop("inverse_track_fn_kw", {})

            # and the potential argument
            try:
                kwargs["potential"]
            except KeyError:
                raise TypeError(
                    "Can't initialize an OrbitSkyOffsetFrame "
                    "without potential= keyword."
                )

            # and the afn arguments
            try:
                afn_bounds = kwargs["afn_bounds"]
            except KeyError:
                raise TypeError(
                    "Can't initialize an OrbitSkyOffsetFrame "
                    "without afn_bounds= keyword."
                )
            else:
                if len(afn_bounds) != 2:
                    raise ValueError("`afn_bounds` must be len= 2.")

            if hasattr(origin_frame, "frame"):
                origin_frame = origin_frame.frame

            newcls = make_orbitskyoffset_cls(
                origin_frame.__class__,
                track_fn=track_fn,
                track_fn_kw=track_fn_kw,
                inverse_track_fn=inverse_track_fn,
                inverse_track_fn_kw=inverse_track_fn_kw,
            )

            return newcls.__new__(newcls, *args, **kwargs)

        # /if

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)

        return super().__new__(cls, *args, **kwargs)

    # /def

    @classmethod
    def from_galpy_orbit(
        cls,
        *args,
        orbit,
        orbit_bkw=None,
        frame="galactocentric",
        method: T.Union[
            Literal["closest"],
            Literal["linear"],
            Literal["cubic"],
            T.Callable[[T.Sequence], T.Sequence],
        ] = "closest",
        time_unit=None,
        **kwargs,
    ):
        """Create an Orbit Sky-Offset Frame from a galpy orbit.

        Parameters
        ----------
        orbit : `~galpy.orbit.Orbit`
            An integrated single orbit., keyword only
            the initial 4-vector and conjugate momenta are taken as the origin.
        orbit_bkw : `~galpy.orbit.Orbit`, optional, keyword only
            An integrated single orbit in the opposite time direction.
            This allows for a leading and trailing orbit from the origin point.
            Must have same origin as `orbit`, but be integrated in reverse.
        frame : str or BaseCoordinateFrame, optional, keyword only
            frame in which to represent the orbit.
            calls ``transform_to(frame)`` on `orbit`'s ICRS SkyCoord output,
            so be careful about things like Galactocentric defaults
        method : {"closest", "linear", "cubic"} or Callable, optional, keyword
            how to construct the affine parameter function mapping time to
            coordinate The orbit integration is precomputed at discrete time
            intervals and the path frame needs to be able to match coordinates
            to the closest point on the orbit. This can be done by treating
            the orbit as an immutable catalog (option "closest", default), a
            linearly interpolatable set of points (option "linear") with
            :class:`~scipy.interpolate.interp1d`, a cubic interpolation
            (option "cubic") with :class:`~scipy.interpolate.CubicSpline`,
            or any user-provided univariate function like those of "linear"
            and "cubic".

            .. todo::

                minimization set by "tol" parameter for closest point on curve
        time_unit : Quantity or Nont, optional, keyword only
            preferred time unit. None means no modification.
            Galpy defaults to Gyr.

        Raises
        ------
        ValueError
            if `orbit` is not integrated
            if `orbit_bkw` is not integrated (and not None)
            if the potential in `orbit` and `orbit_bkw` do not match.

        Notes
        -----
        Make sure that the orbit does not wrap and get close to itself.
        There are currently no checks that this is happening and this
        can lead to some very strange coordinate projections.

        .. todo::

            - allow affine parameter to be arc length
            - break out into a function that calls OrbitSkyOffsetFrame
              so not a classmethod. Keeps this class general.

        """
        # ----------------------------------------------------
        # Checks
        # if orbit_bkw is not None, check has potential and matches orbit
        # check times go in opposite directions

        try:  # check orbit has potential
            orbit._pot
        except AttributeError:  # it does not
            raise ValueError("`orbit` must be integrated.")
        else:  # it does
            # grab the potential, integration time, and origin
            potential = orbit._pot
            t_fwd = orbit.time(use_physical=True)
            origin = orbit.SkyCoord().transform_to(frame)

        if orbit_bkw is not None:
            try:
                orbit._pot
            except AttributeError:
                raise ValueError("`orbit_bkw` must be integrated.")

            if orbit._pot != potential:
                raise ValueError(
                    ("potential in `orbit` and `orbit_bkw` do not match.")
                )

            # check time "directions" are opposite
            # they "fwd" direction can be back in time. That is permitted
            # just not allowed to have the "bkw" direction be the same.
            time_bkw_sgn = np.sign(orbit_bkw.t[1] - orbit_bkw.t[0])
            t_fwd_sgn = np.sign(orbit.t[1] - orbit.t[0])

            if time_bkw_sgn == t_fwd_sgn:
                raise ValueError(
                    (
                        "`orbit` and `orbit_bkw` must be integrated"
                        "in opposite time directions"
                    )
                )

            # get back time, converting to correct units
            t_bkw = orbit_bkw.time(use_physical=True)[::-1] << t_fwd.unit

            # concatenate fwd and bkw orbits into single orbit catalog
            orbit_catalog = concatenate(
                [
                    orbit_bkw.SkyCoord(t_bkw).transform_to(frame).frame,
                    orbit.SkyCoord(t_fwd).transform_to(frame).frame,
                ]
            )
            orbit_time = np.concatenate((t_bkw, t_fwd))

            # create time bounds
            _u = t_fwd.unit if time_unit is None else time_unit
            if t_fwd[-1] > t_bkw[-1]:
                t_bnds = [t_bkw[0], t_fwd[-1]] << _u
            else:
                t_bnds = [t_fwd[0], t_bkw[-1]] << _u

        else:
            # create orbit catalog
            orbit_catalog = orbit.SkyCoord(t_fwd).transform_to(frame)
            orbit_time = t_fwd

            # create time bounds
            _u = t_fwd.unit if time_unit is None else time_unit
            if t_fwd[-1] > t_fwd[0]:
                t_bnds = [t_fwd[0], t_fwd[-1]] << _u
            else:
                t_bnds = [t_fwd[-1], t_fwd[0]] << _u

        # convert orbit time to `time_unit`, if specified
        if time_unit is not None:
            orbit_time <<= time_unit  # (in-place modification)
        else:  # time unit is not None
            time_unit = orbit_time.unit

        # ----------------------------------------------------
        # construct affine function

        track_fn_kw = kwargs.pop("track_fn_kw", {"afn_name": "time"})
        inverse_track_fn_kw = kwargs.pop("inverse_track_fn_kw", {})

        if isinstance(method, str):
            # now need to check it's one of the supported strings
            if method.lower() == "closest":
                # does a catalog match between the coordinates and the
                # points on the orbit from the orbit integration
                # track_fn = ("closest", orbit_catalog, orbit_time)
                track_fn = catalog_match_track
                inverse_track_fn = catalog_match_inverse_track
                track_fn_kw = {
                    "catalog": orbit_catalog,
                    "affine_param": orbit_time,
                    "adj_sep_sgn": True,
                    "afn_name": "time",
                }
                inverse_track_fn_kw = {
                    "catalog": orbit_catalog,
                    "affine_param": orbit_time,
                }

                _interpolation_flag = False

            else:
                # need to handle interpolated functions separately,
                # because requires a closest point "optimization"
                if method.lower() == "linear":
                    method = interpolate.interp1d
                elif method.lower() == "cubic":
                    method = interpolate.CubicSpline
                else:
                    raise ValueError(f"method {method} not known.")

                _interpolation_flag = True

        elif callable(method):
            _interpolation_flag = True

        else:
            raise ValueError(f"method {method} not known.")

        # /if

        if _interpolation_flag:

            # get affine parameter and data interpolation ready
            affine_param = orbit_time.to_value(time_unit)
            _data = orbit_catalog.data._values
            _data = _data.view(np.float64).reshape(_data.shape + (-1,))

            # construct interpolation
            _track_array_fn = method(affine_param, _data.T)

            # astropy coordinate object reconstruction information
            _cmpt = [
                (c, orbit_catalog.data._units[c],)  # in right order, but dict
                for c in orbit_catalog.data.components  # tuple = order
            ]
            _frame = orbit_catalog.frame.realize_frame(None)
            _rep = orbit_catalog.frame.get_representation_cls()

            # interpolation function as astropy, not numpy
            def _track_fn(affine_param):
                """_track_array_fn converted back into a coordinate object."""
                _oc = _track_array_fn(affine_param)  # evaluate interpolation
                rep = _rep(**{c: _oc[i] * u for i, (c, u) in enumerate(_cmpt)})
                catalog = SkyCoord(  # make catalog (TODO not SkyCoord)
                    _frame.realize_frame(rep)
                )
                return catalog

            # make actual track and inverse functions
            def track_fn(
                coords, tol=None, init_sampler: T.Union[float, int] = 1e4
            ):
                """Map coordinates to catalog projection.

                .. todo::

                    change defualt `tol` to something else

                Parameters
                ----------
                coords: SkyCoord
                tol : float or None, optional
                    If None (default), does catalog match but no further
                    minimization. The catalog is the evaluation of the "method"
                    function with affine parameter linearly sampled with
                    `init_sampler` points between "afn_bounds"
                init_sampler : int or float, optional
                    the number of points in ``np.linspace`` for an inital
                    sampling of the affine parameter.

                """
                _aff = np.linspace(  # affine parameter
                    *t_bnds, num=int(init_sampler)
                )[1:-2]
                catalog = _track_fn(_aff)

                if tol is None:
                    return catalog_match_track(
                        coords,
                        catalog=catalog,
                        affine_param=_aff,
                        adj_sep_sgn=True,
                    )
                else:  # TODO actual minimization
                    # initial guess
                    idx, sep2d, _, = match_coordinates_sky(coords, catalog)
                    raise ValueError("Not yet implemented")

                return idx

            # /def

            def inverse_track_fn(coords, **kw):
                """Map catalog projection to coordinates.

                .. todo::

                    this is a very generic function. put somewhere else.

                Parameters
                ----------
                coords: SkyCoord
                    in orbit projection frame

                """
                orbit_pos = _track_fn(coords.afn)

                # need to know offset direction. Most should have _PA
                pa = coords.data._PA
                # Now offset by `sep` in direction `pa`
                out = orbit_pos.directional_offset_by(
                    pa, np.abs(coords.sep)  # need abs() b/c `adj_sep_sgn`
                ).represent_as("spherical")

                return out.lon, out.lat, coords.distance

            # /def
        # /if

        # TODO correct construction with init to get the Attributes
        self = cls(
            *args,
            origin=origin,
            potential=potential,
            track_fn=track_fn,
            inverse_track_fn=inverse_track_fn,
            track_fn_kw=track_fn_kw,
            inverse_track_fn_kw=inverse_track_fn_kw,
            afn_bounds=t_bnds,
            **kwargs,
        )

        return self

    # /def

    # ---------------------------------------------------------------

    def __init__(self, *args, **kwargs):
        """Initialize OrbitSkyOffsetFrame."""
        # remove arguments into __new__ that are not supported in __init__.
        kwargs.pop("track_fn", None)
        kwargs.pop("inverse_track_fn", None)
        kwargs.pop("track_fn_kw", None)
        kwargs.pop("inverse_track_fn_kw", None)

        afn_bounds = kwargs.pop("afn_bounds", None)
        if afn_bounds is not None:
            kwargs["afn_bound_tail"], kwargs["afn_bound_lead"] = afn_bounds

        # initialize
        super().__init__(*args, **kwargs)

        if self.origin is not None and not self.origin.has_data:
            raise ValueError(
                "The origin supplied to OrbitSkyOffsetFrame has no data."
            )

    # /def


# /class


##############################################################################
# END

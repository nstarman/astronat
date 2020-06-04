# -*- coding: utf-8 -*-

"""Offset Frame.

Extending SkyOffset frames to 3D

https://github.com/astropy/astropy/blob/59df4c077c3a0530af4735c44be3db8f2083ad6e/astropy/coordinates/builtin_frames/skyoffset.py
https://docs.astropy.org/en/stable/coordinates/matchsep.html#astropy-skyoffset-frames
https://docs.astropy.org/en/stable/api/astropy.coordinates.builtin_frames.Galactocentric.html#astropy.coordinates.Galactocentric

.. todo::

    - use AffineTransform for the transformation

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["astropy"]

# __all__ = [
#     ""
# ]


##############################################################################
# IMPORTS

# THIRD PARTY

from astropy import units as u
from astropy.utils.decorators import format_doc
from astropy.coordinates.angles import Angle
from astropy.coordinates.matrix_utilities import (
    rotation_matrix,
    matrix_product,
    matrix_transpose,
)
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import (
    BaseCoordinateFrame,
    frame_transform_graph,
    base_doc,
)
from astropy.coordinates.attributes import (
    CoordinateAttribute,
    QuantityAttribute,
    DifferentialAttribute,
)
from astropy.coordinates.transformations import (
    AffineTransform,
    DynamicMatrixTransform,
    FunctionTransform,
)
from astropy.coordinates.errors import ConvertError
from astropy.coordinates import ICRS


# PROJECT-SPECIFIC

__all__ = ["OffsetFrame"]


##############################################################################
# PARAMETERS

_offset_cache = {}

##############################################################################
# CODE
##############################################################################


def make_offset_cls(framecls):
    """
    Create a new class that is the sky offset frame for a specific class of
    origin frame. If such a class has already been created for this frame, the
    same class will be returned.
    The new class will always have component names for spherical coordinates of
    ``lon``/``lat``.
    Parameters
    ----------
    framecls : coordinate frame class (i.e., subclass of `~astropy.coordinates.BaseCoordinateFrame`)
        The class to create the OffsetFrame of.
    Returns
    -------
    offsetframecls : class
        The class for the new offset frame.
    Notes
    -----
    This function is necessary because Astropy's frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct offset frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given python session.
    """

    if framecls in _offset_cache:
        return _offset_cache[framecls]

    # the class of a class object is the metaclass
    framemeta = framecls.__class__

    class OffsetMeta(framemeta):
        """
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
            # OffsetFrame class initially.
            members[
                "_default_representation"
            ] = framecls._default_representation
            members["_default_differential"] = framecls._default_differential

            newname = name[:-5] if name.endswith("Frame") else name
            newname += framecls.__name__

            return super().__new__(cls, newname, bases, members)

    # We need this to handle the intermediate metaclass correctly, otherwise we could
    # just subclass OffsetFrame.
    _OffsetFramecls = OffsetMeta(
        "OffsetFrame",
        (OffsetFrame, framecls),
        {"__doc__": OffsetFrame.__doc__},
    )

    @frame_transform_graph.transform(
        FunctionTransform, _OffsetFramecls, _OffsetFramecls
    )
    def offset_to_offset(from_offset_coord, to_offset_frame):
        """Transform between two offset frames."""

        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
        intermediate_from = from_offset_coord.transform_to(
            from_offset_coord.origin
        )
        intermediate_to = intermediate_from.transform_to(
            to_offset_frame.origin
        )
        return intermediate_to.transform_to(to_offset_frame)

    @frame_transform_graph.transform(
        DynamicMatrixTransform, framecls, _OffsetFramecls
    )
    def reference_to_offset(reference_frame, offset_frame):
        """Convert a reference coordinate to an sky offset frame."""

        # Define rotation matrices along the position angle vector, and
        # relative to the origin.
        origin = offset_frame.origin.spherical
        mat1 = rotation_matrix(-offset_frame.rotation, "x")
        mat2 = rotation_matrix(-origin.lat, "y")
        mat3 = rotation_matrix(origin.lon, "z")
        return matrix_product(mat1, mat2, mat3)

    @frame_transform_graph.transform(
        DynamicMatrixTransform, _OffsetFramecls, framecls
    )
    def offset_to_reference(offset_coord, reference_frame):
        """Convert an sky offset frame coordinate to the reference frame"""

        # use the forward transform, but just invert it
        R = reference_to_offset(reference_frame, offset_coord)
        # transpose is the inverse because R is a rotation matrix
        return matrix_transpose(R)

    _offset_cache[framecls] = _OffsetFramecls
    return _OffsetFramecls

    # def get_matrix_vectors(offset_frame, inverse=False):
    #     """
    #     Use the ``inverse`` argument to get the inverse transformation, matrix and
    #     offsets to go from OffsetFrame to ICRS.
    #     """
    #     # shorthand
    #     of = offset_frame

    #     # rotation matrix to align x(ICRS) with the vector to the Galactic center
    #     mat1 = rotation_matrix(-of.origin_coord.dec, "y")
    #     mat2 = rotation_matrix(of.origin_coord.ra, "z")
    #     # extra roll away from the Galactic x-z plane
    #     mat0 = rotation_matrix(of.get_roll0() - of.roll, "x")

    #     # construct transformation matrix and use it
    #     R = matrix_product(mat0, mat1, mat2)

    #     # Now need to translate by Sun-Galactic center distance around x' and
    #     # rotate about y' to account for tilt due to Sun's height above the plane
    #     translation = r.CartesianRepresentation(
    #         of.galcen_distance * [1.0, 0.0, 0.0]
    #     )
    #     z_d = of.z_sun / of.galcen_distance
    #     H = rotation_matrix(-np.arcsin(z_d), "y")

    #     # compute total matrices
    #     A = matrix_product(H, R)

    #     # Now we re-align the translation vector to account for the Sun's height
    #     # above the midplane
    #     offset = -translation.transform(H)

    #     if inverse:
    #         # the inverse of a rotation matrix is a transpose, which is much faster
    #         #   and more stable to compute
    #         A = matrix_transpose(A)
    #         offset = (-offset).transform(A)
    #         offset_v = r.CartesianDifferential.from_cartesian(
    #             (-of.galcen_v_sun).to_cartesian().transform(A)
    #         )
    #         offset = offset.with_differentials(offset_v)

    #     else:
    #         offset = offset.with_differentials(of.galcen_v_sun)

    #     return A, offset

    # def _check_coord_repr_diff_types(c):
    #     if isinstance(c.data, r.UnitSphericalRepresentation):
    #         raise ConvertError(
    #             "Transforming to/from an OffsetFrame "
    #             "requires a 3D coordinate, e.g. (angle, angle, "
    #             "distance) or (x, y, z)."
    #         )

    #     if "s" in c.data.differentials and isinstance(
    #         c.data.differentials["s"],
    #         (
    #             r.UnitSphericalDifferential,
    #             r.UnitSphericalCosLatDifferential,
    #             r.RadialDifferential,
    #         ),
    #     ):
    #         raise ConvertError(
    #             "Transforming to/from an OffsetFrame "
    #             "requires a 3D velocity, e.g., proper motion "
    #             "components and radial velocity."
    #         )

    # @frame_transform_graph.transform(AffineTransform, ICRS, OffsetFrame)
    # def icrs_to_galactocentric(icrs_coord, offset_frame):
    #     _check_coord_repr_diff_types(icrs_coord)
    #     return get_matrix_vectors(offset_frame)

    # @frame_transform_graph.transform(AffineTransform, OffsetFrame, ICRS)
    # def galactocentric_to_icrs(galactocentric_coord, icrs_frame):
    #     _check_coord_repr_diff_types(galactocentric_coord)
    #     return get_matrix_vectors(galactocentric_coord, inverse=True)


##############################################################################


class OffsetFrame(BaseCoordinateFrame):
    """Offset Frame.

    A frame which is relative to some specific position and oriented to match
    its frame. OffsetFrames always have component names for spherical
    coordinates of ``lon``/``lat``, *not* the component names for the frame of
    ``origin``. This is useful for calculating offsets and dithers in the
    frame of the sky relative to an arbitrary position. Coordinates in this
    frame are both centered on the position specified by the ``origin``
    coordinate, *and* they are oriented in the same manner as the ``origin``
    frame.  E.g., if ``origin`` is `~astropy.coordinates.ICRS`, this object's
    ``lat`` will be pointed in the direction of Dec, while ``lon`` will point
    in the direction of RA. For more on offset frames, see
    :ref:`astropy-offset-frames`.

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    origin : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the origin of this frame. Note that this
        origin is used purely for on-sky location/rotation.  It can have a
        ``distance`` but it will not be used by this ``OffsetFrame``.
    rotation : `~astropy.coordinates.Angle` or `~astropy.units.Quantity` with angle units
        The final rotation of the frame about the ``origin``. The sign of
        the rotation is the left-hand rule.  That is, an object at a
        particular position angle in the un-rotated system will be sent to
        the positive latitude (z) direction in the final frame.

    Notes
    -----
    ``OffsetFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``OffsetFrame``.  Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``origin``.
    """

    rotation = QuantityAttribute(default=0, unit=u.deg)
    origin = CoordinateAttribute(default=None, frame=None)

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # an offset frame for this class.
        if not (issubclass(cls, OffsetFrame) and cls is not OffsetFrame):
            # We get the origin argument, and handle it here.
            try:
                origin_frame = kwargs["origin"]
            except KeyError:
                raise TypeError(
                    "Can't initialize an OffsetFrame without origin= keyword."
                )
            if hasattr(origin_frame, "frame"):
                origin_frame = origin_frame.frame
            newcls = make_offset_cls(origin_frame.__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.origin is not None and not self.origin.has_data:
            raise ValueError(
                "The origin supplied to OffsetFrame has no " "data."
            )
        if self.has_data:
            self._set_offset_data_lon_wrap_angle(self.data)

    @staticmethod
    def _set_offset_data_lon_wrap_angle(data):
        if hasattr(data, "lon"):
            data.lon.wrap_angle = 180.0 * u.deg
        return data

    def represent_as(self, base, s="base", in_frame_units=False):
        """
        Ensure the wrap angle for any spherical
        representations.
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        self._set_offset_data_lon_wrap_angle(data)
        return data


##############################################################################
# END

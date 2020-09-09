# -*- coding: utf-8 -*-

"""Attributes.

:mod:`~astropy.coordinates.attributes`

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "OrbitRepresentationBase",
    "OrbitSkyOffsetRepresentation",
    "OrbitOffsetCartesianRepresentation",
    "OrbitOffsetCylindricalRepresentation",
]


##############################################################################
# IMPORTS

import abc
import inspect
import typing as T
from collections import OrderedDict

import astropy.units as u
import numpy as np
from astropy.coordinates import BaseDifferential, BaseRepresentation
from astropy.coordinates.distances import Distance
from astropy.utils import classproperty

##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


class OrbitRepresentationBase(BaseRepresentation):
    """Base class for Orbit representations.

    first parameter is the affine parameter, generally time or arc-length

    .. todo:

        - support an inverse function t(d), where d is the distance along the
          orbit, so can support the arc length as the affine parameter.
        - support hidden variables, like _PA

    """

    attr_classes = OrderedDict([("afn", u.Quantity)])
    """Attribute Classes. Should be an OrderedDict."""

    def __new__(cls, *args, **kwargs):
        afn_name = kwargs.pop("afn_name", None)
        # print("afn_name: ", afn_name)

        if afn_name is not None:
            setattr(cls, afn_name, property(lambda self: self.afn))

        self = super().__new__(cls)
        self._afn_name = afn_name

        return self

    def __init__(self, afn, *args, **kwargs):
        """OrbitRepresentationBase, afn along orbit as affine parameter."""
        super().__init__(afn, *args, **kwargs)

    # /def

    @abc.abstractmethod
    def from_cartesian(self, cartesian):  # TODO
        """From Cartesian. Abstract method. Should be classmethod."""
        pass

    # /def

    @abc.abstractmethod
    def to_cartesian(self):  # TODO
        """To Cartesian. Abstract method."""
        pass

    # /def


# /class


# -------------------------------------------------------------------


class OrbitSkyOffsetUnitRepresentation(OrbitRepresentationBase):
    """Define a Pseudo-Spherical projected-Orbit on-sky coordinate system.

    Parameterized by the:

        - afn, the affine parameter along the orbit (instantiation coordinate)
        - sep, the on-sky angular separation from the orbit at position `afn`
        - distance, the radial distance from coordinate center.


    .. todo::

        - Allow differentials
        - OrbitOffsetCartesianRepresentation transformation
        - OrbitOffsetCylindricalRepresentation transformation

    Parameters
    ----------
    afn : Quantity
    sep : Angle or Quantity
    _PA : Angle or Quantity


    """

    attr_classes = OrderedDict(
        [
            ("afn", u.Quantity),  # affine parameter
            ("sep", u.Quantity),  # sky separation
            ("_PA", u.Quantity),  # position-angle, hidden variable
        ]
    )

    @classproperty
    def _dimensional_representation(cls):
        return OrbitSkyOffsetRepresentation

    # @u.quantity_input(sep="angle", _PA="angle")
    def __init__(
        self,
        afn,
        sep,
        _PA=np.NaN * u.deg,
        *,
        differentials=None,
        copy: bool = True,
        afn_name: T.Optional[str] = None
    ):
        if differentials is not None:  # TODO, allow differentials
            raise ValueError()

        super().__init__(
            afn, sep, _PA=_PA, copy=copy, differentials=differentials
        )

    # /def

    def scale_factors(self):
        """Physical scale factor in component directions.

        Returns
        -------
        Returns a dict with a Quantity for each component with the appropriate
        physical scale factor for a unit change in that direction.

        """
        # TODO is self.afn.unit fixed?
        sf_afn = np.broadcast_to(1.0 / self.afn.unit, self.shape, subok=True)
        sf_sep = np.broadcast_to(1.0 / u.radian, self.shape, subok=True)

        return OrderedDict((("afn", sf_afn), ("sep", sf_sep)))

    # /def

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is an AlongOrbit representation

        # TODO: shortcut even if a differential_class is passed in,
        # using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(
                other_class, OrbitSkyOffsetRepresentation
            ):  # TODO differential
                return OrbitSkyOffsetRepresentation(
                    afn=self.afn,
                    sep=self.sep,
                    distance=1 * u.dimensionless_unscaled,
                    _PA=self._PA,
                )
            elif issubclass(other_class, OrbitOffsetCartesianRepresentation):
                raise Exception("IOU transformation")
            elif issubclass(other_class, OrbitOffsetCylindricalRepresentation):
                raise Exception("IOU transformation")

        return super().represent_as(other_class, differential_class)

    # /def

    @classmethod
    def from_cartesian(self, cartesian):  # TODO
        raise Exception("There is no cartesian representation.")

    # /def

    def to_cartesian(self):  # TODO
        raise Exception("There is no cartesian representation.")

    # /def


# /class

# -------------------------------------------------------------------


class OrbitSkyOffsetRepresentation(OrbitRepresentationBase):
    """Define a Pseudo-Spherical projected-Orbit on-sky coordinate system.

    Parameterized by the:

        - afn, the afn along the orbit (from instantiation coordinate)
        - sep, the on-sky angular separation from the orbit at position afn
        - distance, the radial distance from coordinate center.


    .. todo::

        - Allow differentials
        - OrbitOffsetCartesianRepresentation transformation
        - OrbitOffsetCylindricalRepresentation transformation

    """

    attr_classes = OrderedDict(
        [
            ("afn", u.Quantity),
            ("sep", u.Quantity),  # sky separation
            ("distance", u.Quantity),
            ("_PA", u.Quantity),  # position-angle, hidden variable
        ]
    )
    _unit_representation = OrbitSkyOffsetUnitRepresentation

    # @u.quantity_input(sep="angle", _PA="angle")
    def __init__(
        self,
        afn,
        sep,
        distance,
        _PA=np.NaN * u.deg,
        *,
        differentials=None,
        copy: bool = True,
        afn_name: T.Optional[str] = None
    ):
        if differentials is not None:  # TODO, allow differentials
            raise ValueError()

        super().__init__(
            afn,
            sep,
            distance,
            _PA=_PA,
            copy=copy,
            differentials=differentials,
        )

        if self._distance.unit.physical_type == "length":
            try:
                self._distance = Distance(self._distance, copy=False)
            except ValueError as e:
                if e.args[0].startswith("Distance must be >= 0"):
                    raise ValueError(
                        "Distance must be >= 0. To allow negative "
                        "distance values, you must explicitly pass"
                        " in a `Distance` object with the the "
                        "argument 'allow_negative=True'."
                    )
                else:
                    raise

    # /def

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is an AlongOrbit representation

        # TODO: shortcut even if a differential_class is passed in,
        # using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(
                other_class, OrbitSkyOffsetUnitRepresentation
            ):  # TODO differential
                return OrbitSkyOffsetUnitRepresentation(
                    afn=self.afn, sep=self.sep, _PA=self._PA,
                )
            elif issubclass(other_class, OrbitOffsetCartesianRepresentation):
                raise Exception("IOU transformation")
            elif issubclass(other_class, OrbitOffsetCylindricalRepresentation):
                raise Exception("IOU transformation")

        return super().represent_as(other_class, differential_class)

    # /def

    @classmethod
    def from_cartesian(self, cartesian):
        raise Exception("There is no cartesian representation.")

    # /def

    def to_cartesian(self):
        raise Exception("There is no cartesian representation.")

    # /def

    def unit_vectors(self):
        """Unit vectors in component directions.

        Returns
        -------
        Returns a dict with a CartesianRepresentation of unit vectors
        in the direction of each component.

        """
        # return {'comp1': CartesianRepresentation(...),
        #         'comp2': CartesianRepresentation(...),
        #         'comp3': CartesianRepresentation(...)}
        raise Exception("Not yet implemented")

    # /def

    def scale_factors(self):
        """Physical scale factor in component directions.

        Returns
        -------
        Returns a dict with a Quantity for each component with the appropriate
        physical scale factor for a unit change in that direction.

        """
        # TODO is self.afn.unit fixed?
        sf_afn = np.broadcast_to(1.0 / self.afn.unit, self.shape, subok=True)
        sf_sep = np.broadcast_to(1.0 / u.radian, self.shape, subok=True)
        sf_distance = np.broadcast_to(1.0 * u.one, self.shape, subok=True)

        return OrderedDict(
            (("afn", sf_afn), ("sep", sf_sep), ("distance", sf_distance))
        )

    # /def


# /class

# -------------------------------------------------------------------


class OrbitSkyOffsetDifferential(BaseDifferential):  # TODO
    base_representation = OrbitSkyOffsetRepresentation


# /class


# -------------------------------------------------------------------


class OrbitOffsetCartesianRepresentation(OrbitRepresentationBase):
    """Define a pseudo-Cartesian along-Orbit coordinate system.

    A tube around the orbit.

    Parameterized by the:

        - afn, the affine parameter along the orbit (instantiation coordinate)
        - x, the distance from the orbit, in the plane of the orbit
        - y, the distance from the orbit, perpendicular to the plane

    .. todo::

        - Allow differentials
        - get correct phi value in `represent_as`
        - OrbitSkyOffsetRepresentation transformation
        - dynamically defined by the orbit, for Cartesian transformations.
           maybe make this a base-class with abstract methods

    """

    attr_classes = OrderedDict(
        [
            ("afn", u.Quantity),
            ("x", u.Quantity),
            ("y", u.Quantity),
            ("_d_afn", u.Quantity),
        ]
    )

    # @u.quantity_input(x="length", y="length")
    def __init__(
        self,
        afn,
        x,
        y,
        _d_afn=np.NaN * u.pc,
        differentials=None,
        copy: bool = True,
        afn_name: T.Optional[str] = None,
    ):
        if differentials is not None:  # TODO, allow differentials
            raise ValueError()
        super().__init__(
            afn, x, y, _d_afn=_d_afn, copy=copy, differentials=differentials
        )

    # /def

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is an OrbitOffser representation

        # TODO: shortcut even if a differential_class is passed in,
        # using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, OrbitOffsetCylindricalRepresentation):
                return other_class(
                    afn=self.afn,
                    rho=np.sqrt(self.x ** 2 + self.y ** 2),
                    phi=np.arctan2(self.y, self.x),  # TODO correct sign?
                    _d_afn=self._d_afn,
                    copy=False,
                    afn_name=self._afn_name,
                )
            elif issubclass(other_class, OrbitSkyOffsetRepresentation):
                raise Exception("IOU transformation")

        return super().represent_as(other_class, differential_class)

    # /def

    @classmethod
    def from_cartesian(self, cartesian):  # TODO
        raise Exception("There is no cartesian representation, currently.")

    # /def

    def to_cartesian(self):  # TODO
        # Step 1: use afn to get orbit coordinate at that afn value
        # Step 2: use `x` and `y` offset to get real-space position
        # Step 3: representation machinery to return Cartesian
        raise Exception("There is no cartesian representation, currently.")

    # /def


# /class

# -------------------------------------------------------------------


class OrbitOffsetCartesianDifferential(BaseDifferential):  # TODO
    base_representation = OrbitSkyOffsetRepresentation


# /class


# -------------------------------------------------------------------


class OrbitOffsetCylindricalRepresentation(OrbitRepresentationBase):
    """Define a pseudo-Cartesian along-Orbit coordinate system.

    Parameterized by the:

        - afn, the afn along the orbit (from instantiation coordinate)
        - rho, the distance from the orbit
        - phi, the angle around from the orbit
            0 pointing toward observer from orbit start


    .. todo::

        - Allow differentials
        - make sure x, y conversions correct in `represent_as`
        - dynamically defined by the orbit, for Cartesian transformations.
           maybe make this a base-class with abstract methods

    Parameters
    ----------
    afn : Quantity
    rho : Quantity
    phi : Quantity
    _d_afn : Quantity
    differentials
    copy : bool
    afn_name : str, optional

    """

    attr_classes = OrderedDict(
        [
            ("afn", u.Quantity),
            ("rho", u.Quantity),
            ("phi", u.Quantity),
            ("_d_afn", u.Quantity),
        ]
    )

    # @u.quantity_input(rho="length", phi="angle")
    def __init__(
        self,
        afn,
        rho,
        phi,
        _d_afn=np.NaN * u.pc,
        differentials=None,
        copy: bool = True,
        afn_name: T.Optional[str] = None,
    ):
        """Initialize class."""
        if differentials is not None:  # TODO, allow differentials
            raise ValueError()
        super().__init__(
            afn,
            rho,
            phi,
            _d_afn=_d_afn,
            copy=copy,
            differentials=differentials,
        )

    # /def

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is an AlongOrbit representation

        # TODO: shortcut even if a differential_class is passed in,
        # using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, OrbitOffsetCartesianRepresentation):
                return other_class(
                    afn=self.afn,
                    x=self.rho * np.cos(self.phi),  # TODO get correct
                    y=self.rho * np.sin(self.phi),  # TODO get correct
                    _d_afn=self._d_afn,
                    copy=False,
                    afn_name=self._afn_name,
                )
            elif issubclass(other_class, OrbitSkyOffsetRepresentation):
                raise Exception("IOU transformation")

        return super().represent_as(other_class, differential_class)

    # /def

    @classmethod
    def from_cartesian(self, cartesian):
        """From Cartesian.

        .. todo::

            implement

        """
        raise Exception("There is no cartesian representation, currently.")

    # /def

    def to_cartesian(self):
        """To Cartesian.

        .. todo::

            implement

        """
        # Step 1: use afn to get orbit coordinate at that afn value
        # Step 2: use `x` and `y` offset to get real-space position
        # Step 3: representation machinery to return Cartesian
        raise Exception("There is no cartesian representation, currently.")

    # /def


# /class

# -------------------------------------------------------------------


class OrbitOffsetCylindricalDifferential(BaseDifferential):
    """Orbit-Offset Cylindrical Differential.

    .. todo::

        implement

    """

    base_representation = OrbitSkyOffsetRepresentation


# /class


##############################################################################
# END

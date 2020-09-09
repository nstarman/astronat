# -*- coding: utf-8 -*-

"""Galpy-Gala Interface for Initial Conditions.

Both `galpy` and `gala` use a Galactocentric (Gc) reference frame under
the hood. The `galpy` Orbit instantiation is very flexible, but also
requires some care since `galpy`'s Gc frame conversions are different than
Astropy's. The provided solution is to either pass a Gc frame class
(preferred), or if using a SkyCoord in another reference frame, to provide
the conversion information. `galpy` will use the information in leiu of
its own conversions. This is described `here
<https://docs.galpy.org/en/v1.5.0/orbit.html#orbinit>`_ and shown.

    >>> from astropy.coordinates import CartesianDifferential, SkyCoord
    >>> from galpy.orbit import Orbit
    >>> c = SkyCoord(
    ...     ra=20.*u.deg,dec=30.*u.deg,distance=2.*u.kpc,
    ...     pm_ra_cosdec=-10.*u.mas/u.yr,pm_dec=20.*u.mas/u.yr,
    ...     radial_velocity=50.*u.km/u.s,
    ...     galcen_distance=8.*u.kpc, z_sun=15.*u.pc,
    ...     galcen_v_sun=CartesianDifferential([10.0,235.,7.]*u.km/u.s)
    ... )
    >>> o = Orbit(c)

The current solution in *gala* for instantiating a *PhaseSpacePosition*,
the precursor to an orbit, is to demand a GC coordinates as an array or
BaseRepresentation. The Galactocentric properties, like ``galcen_v_sun``,
are determined by Astropy's ``galactocentric_frame_defaults``
ScienceState. Note, this implementation might change as the "frame" key of
PhaseSpacePosition is developed.

In order to maintain consistency between the two packages we write a
custom PhaseSpacePosition that accepts SkyCoords and BaseCoordinateFrame
objects (like GC) in addition to `gala`'s standard inputs. The SkyCoords
and BaseCoordinateFrames are stored and internally converted to the GC
frame. This conversion is done with the
:class:`~astropy.coordinates.Galactocentric`, which is set by
:func:`~astropy.coordinates.galactocentric_frame_defaults`. If a SkyCoord
has any of the GC frame properties (``galcen_coord``, ``roll``,
``galcen_distance``, ``galcen_v_sun``, ``z_sun``), these are used instead
of the defaults, permitting the above :mod:`~galpy` example.

Notes
-----
.. todo::

    - It would be nice to have a more general PhaseSpacePosition class,
      like a PhaseSpacePositionBundle that can contain many difference IC PSPs,
      like a PhotometricSpacePosition class.

    - functions to generate a PhaseSpacePosition, like for a stream

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["Jo Bovy", "Adrian Price-Whelan"]


__all__ = ["PhaseSpacePosition"]


##############################################################################
# IMPORTS

import typing as T

import astropy.coordinates as coords

try:
    import gala
except ImportError:
    _HAS_GALA = False
else:
    gala.__version__
    _HAS_GALA = True

from utilipy.utils.collections import ReferenceBase

##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################

if _HAS_GALA:

    from gala.dynamics import PhaseSpacePosition as PhaseSpacePositionBase

else:

    class PhaseSpacePositionBase:
        """Gala-style Initial Conditions Class.

        The initial condition (IC) construction process in :mod:`~galpy` is
        very flexible: the Orbits class accepts anything from a list of floats
        to `~astropy.units.Quantity` arrays to a full
        `~astropy.coordinates.SkyCoord`. Conversely, the :mod:`~gala` package
        implements a single class for ICs, the
        :mod:`~gala.dynamics.PhaseSpacePosition`. To interface between the
        two packages we implement a dummy PhaseSpacePosition that mimics
        the Gala package.

        """

        def __init__(self, pos, vel=None, frame=None):
            """Phase Space Position Proxy.

            Parameters
            ----------
            pos : SkyCoord or BaseCoordinateFrame

            """
            if not isinstance(pos, (coords.Galactocentric,)):
                raise TypeError("currently only support Galactocentric")
            else:
                self.pos = pos

            if vel is not None:
                raise TypeError("currently only support `pos`")

            if frame is not None:
                raise TypeError("currently only support `pos`")

        # /def

        def to_coord_frame(self, frame, galactocentric_frame=None, **kwargs):
            """To coordinate frame.

            Parameters
            ----------
            frame : :class:`~astropy.coordinates.BaseCoordinateFrame`
                frame to which to transform
            galactocentric_frame : None, optional
            kwargs

            Returns
            -------
            pos

            """
            if galactocentric_frame is not None:
                raise TypeError("currently only support `frame`")
            if galactocentric_frame != {}:
                raise TypeError("currently only support `frame`")

            if isinstance(self.pos, (coords.Galactocentric)):
                return self.pos
            else:
                raise TypeError("currently only support Galactocentric")

        # /def

    # /class


# -------------------------------------------------------------------


# @format_doc(
#     None,
#     # pspdoc=PhaseSpacePositionBase.__init__.__doc__
# )
class PhaseSpacePosition(PhaseSpacePositionBase, ReferenceBase):
    """Phase Space Position, for controlling initial conditions.

    Both `galpy` and `gala` use a Galactocentric (Gc) reference frame under
    the hood. The `galpy` Orbit instantiation is very flexible, but also
    requires some care since `galpy`'s Gc frame conversions are different than
    Astropy's. The provided solution is to either pass a Gc frame class
    (preferred), or if using a SkyCoord in another reference frame, to provide
    the conversion information. `galpy` will use the information in leiu of
    its own conversions. This is described `here
    <https://docs.galpy.org/en/v1.5.0/orbit.html#orbinit>`_ and shown.

        >>> from astropy.coordinates import CartesianDifferential, SkyCoord
        >>> from galpy.orbit import Orbit
        >>> c = SkyCoord(
        ...     ra=20.*u.deg,dec=30.*u.deg,distance=2.*u.kpc,
        ...     pm_ra_cosdec=-10.*u.mas/u.yr,pm_dec=20.*u.mas/u.yr,
        ...     radial_velocity=50.*u.km/u.s,
        ...     galcen_distance=8.*u.kpc, z_sun=15.*u.pc,
        ...     galcen_v_sun=CartesianDifferential([10.0,235.,7.]*u.km/u.s)
        ... )
        >>> o = Orbit(c)

    The current solution in *gala* for instantiating a *PhaseSpacePosition*,
    the precursor to an orbit, is to demand a GC coordinates as an array or
    BaseRepresentation. The Galactocentric properties, like ``galcen_v_sun``,
    are determined by Astropy's ``galactocentric_frame_defaults``
    ScienceState. Note, this implementation might change as the "frame" key of
    PhaseSpacePosition is developed.

    In order to maintain consistencey between the two packages we write a
    custom PhaseSpacePosition that accepts SkyCoords and BaseCoordinateFrame
    objects (like GC) in addition to `gala`'s standard inputs. The SkyCoords
    and BaseCoordinateFrames are stored and internally converted to the GC
    frame. This conversion is done with the
    :class:`~astropy.coordinates.Galactocentric`, which is set by
    :func:`~astropy.coordinates.galactocentric_frame_defaults`. If a SkyCoord
    has any of the GC frame properties (``galcen_coord``, ``roll``,
    ``galcen_distance``, ``galcen_v_sun``, ``z_sun``), these are used instead
    of the defaults, permitting the above :mod:`~galpy` example.

    """

    # ::
    #
    #       {pspdoc}

    def __init__(
        self,
        pos,
        vel=None,
        frame=None,
        name: T.Optional[str] = None,
        *,
        reference: T.Union[str, dict, None] = None,
    ):
        """Represent phase-space positions.

        i.e. positions and conjugate momenta (velocities).

        The class can be instantiated with Astropy SkyCoords, coordinate
        objects (e.g., :class:`~astropy.coordinates.Galactocentric`),
        representation objects (e.g.,
        :class:`~astropy.coordinates.CartesianRepresentation`), Astropy
        :class:`~astropy.units.Quantity` objects, or plain Numpy arrays.

        If passing in SkyCoords or coordinate or representation objects,
        the default representation is taken to be the class that is passed in.

        If passing in Quantity or Numpy array instances for both position and
        velocity, they are assumed to be Cartesian. Array inputs are
        interpreted as dimensionless quantities. The input position and
        velocity objects can have an arbitrary number of (broadcastable)
        dimensions. For Quantity or array inputs, the first axis (0) has
        special meaning:

            - `axis=0` is the coordinate dimension (e.g. x, y, z for Cartesian)

        So if the input position array, `pos`, has shape `pos.shape = (3,
        100)`, this would represent 100 3D positions (`pos[0]` is `x`,
        `pos[1]` is `y`, etc.). The same is true for velocity.

        Parameters
        ----------
        pos : SkyCoord or BaseCoordinateFrame or BaseRepresentation or Quantity
            Positions. If a numpy array (e.g., has no units), this will be
            stored as a dimensionless :class:`~astropy.units.Quantity`. See
            the note above about the assumed meaning of this object's axes.
        vel : BaseDifferential or Quantity or array_like
            Velocities. If a numpy array (e.g., has no units), this will be
            stored as a dimensionless :class:`~astropy.units.Quantity`. See
            the note above about the assumed meaning of this object's axes.
        frame : :class:`~gala.potential.FrameBase` (optional)
            The reference frame of the input phase-space positions.
        name : str, optional
            The name of the phase space position, ex Pal 5

        Other Parameters
        ----------------
        reference : str or dict, optional
            stored as a dictionary.

        """
        self.name = name
        observed = None  # start with None, create later
        gc_frame = None  # start with None, create later

        if isinstance(pos, (coords.SkyCoord, coords.BaseCoordinateFrame)):
            # 1) store the "observed" coordinate without modification
            observed = pos

            # 2) get information
            representation_type = pos.representation_type

            # the Galactocentric frame defaults (gcfd)
            gcfd: dict = coords.galactocentric_frame_defaults.get()
            for n in gcfd.keys():  # iterate through Gc frame keys
                # check if object has key (only SkyCoord)
                # this is needed for Galpy's SkyCoord overload method
                if hasattr(pos, n):
                    v = getattr(pos, n)
                    if v is not None:  # only assign if really exists
                        gcfd[n] = v
            # /for

            # 3) Convert to correct Galactocentric BaseRepresentation-type
            pos = pos.transform_to(coords.Galactocentric(**gcfd))

            # store the Galactocentric form
            gc_frame = pos

            pos = pos.represent_as(representation_type, in_frame_units=True)

        # /if

        # make PhaseSpacePosition from base-class
        super().__init__(pos=pos, vel=vel, frame=frame)
        # and add references
        ReferenceBase.__init__(self, reference=reference)

        # set the "observed" coordinates
        if observed is None:  # if passed Gc arrays, not Frame-type
            observed = self.to_coord_frame(coords.Galactocentric)
        self.observed = observed

        # set the Galactocentric BaseCoordinateFrame
        if gc_frame is None:  # if passed Gc arrays, not Frame-type
            gc_frame = self.to_coord_frame(coords.Galactocentric)
        self._gc_frame = gc_frame  # store the frame

        return

    # /def

    @property
    def vxvv(self):
        """The Galpy-style position.

        The :mod:`~galpy` Orbit instantiation is very flexible, but also
        requires some care since Galpy's Gc frame conversions are different
        than :mod:`~astropy`. The provided solution is to pass a
        :class:`~astropy.coordinates.SkyCoord` in the Galactocentric frame.

        Returns
        -------
        :class:`~astropy.coordinates.SkyCoord` in Galactocentric frame.

        """
        return coords.SkyCoord(self._gc_frame)

    # /def


# /class


##############################################################################
# END

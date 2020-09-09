# -*- coding: utf-8 -*-

"""Parameters."""


__all__ = [
    # functions
    "solar_system_vesc_params",
    "vesc_sun_at_R",
    "twobody_vesc",
    "multibody_vesc",
]


##############################################################################
# IMPORTS

import functools
import itertools
import typing as T

import numpy as np
from astropy import units as u
from astropy.utils.decorators import format_doc
from astropy.utils.state import ScienceState
from utilipy.math import as_quantity
from utilipy.utils.typing import QuantityType

##############################################################################
# PARAMETERS

_KMS = u.km / u.s

_ref_B = (
    "Explanatory Supplement to the Astronomical Almanac. "
    "1992. K. P. Seidelmann, Ed., p.706 (Table 15.8) and p.316 "
    "(Table 5.8.1), University Science Books, Mill Valley, California."
)
_ref_C = (
    "Seidelmann, P.K. et al. 2007. 'Report of the IAU/IAG Working "
    "Group on cartographic coordinates and rotational elements: 2006' "
    "Celestial Mech. Dyn. Astr. 98:155-180."
)
_ref_D = (
    "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group "
    "on cartographic coordinates and rotational elements: 2015' "
    "Celestial Mech. Dyn. Astr. 130:22."
)

_sqrt2 = np.sqrt(2)


##############################################################################
# CODE
##############################################################################


class solar_system_vesc_params(ScienceState):
    """Solar System Parameters."""

    _latest_value: str = "default"
    _references: T.Optional[T.Dict[str, T.Any]] = None
    _value: T.Optional[dict] = None

    _registry: T.Dict[str, T.Dict[str, T.Dict[str, T.Any]]] = {
        "DEFAULT": {
            "params": {
                "Sun": 617.5 * _KMS,
                "Mercury": 4.25 * _KMS,
                "Venus": 10.36 * _KMS,
                "Earth": 11.19 * _KMS,
                "Mars": 5.03 * _KMS,
                "Jupiter": 60.20 * _KMS,
                "Saturn": 36.09 * _KMS,
                "Uranus": 21.38 * _KMS,
                "Neptune": 23.56 * _KMS,
                "Pluto": 1.21 * _KMS,
            },
            "references": {
                "_source": "https://ssd.jpl.nasa.gov/?planet_phys_par",
                "Sun": (_ref_B, _ref_C, _ref_D),
                "Mercury": (_ref_B, _ref_C, _ref_D),
                "Venus": (_ref_B, _ref_C, _ref_D),
                "Earth": (_ref_B, _ref_C, _ref_D),
                "Mars": (_ref_B, _ref_C, _ref_D),
                "Jupiter": (_ref_B, _ref_C, _ref_D),
                "Saturn": (_ref_B, _ref_C, _ref_D),
                "Uranus": (_ref_B, _ref_C, _ref_D),
                "Neptune": (_ref_B, _ref_C, _ref_D),
                "Pluto": (_ref_B, _ref_C, _ref_D),
            },
        }
    }

    @classmethod
    def get_solar_params_from_string(cls, arg: str):
        """Get parameters from registry."""
        # Resolve the meaning of 'latest'
        if arg == "latest":
            arg = cls._latest_value

        if arg.lower() == "default":

            info = cls._registry["DEFAULT"]

        elif arg in cls._registry:

            info = cls._registry[arg]

        else:
            raise ValueError(
                f"Invalid string input to retrieve solar "
                f'parameters for Galactocentric frame: "{arg}"'
            )

        return info["params"], info["references"]

    # /def

    @classmethod
    def validate(cls, value: T.Union[None, str, dict]):
        """Validate `value`, from string or dict."""
        if value is None:
            value = cls._latest_value

        if isinstance(value, str):
            params, refs = cls.get_solar_params_from_string(value)
            cls._references = refs
            return params

        elif isinstance(value, dict):
            return value

        else:
            raise ValueError(
                "Invalid input to retrieve solar parameters."
                "Input must be a string, or dict"
            )

    # /def

    @classmethod
    def register(cls, name: str, params: dict, references: dict):
        """Register a set of parameters.

        Parameters
        ----------
        name : str
        params : dict
        references : dict

        """
        cls._registry[name] = {"params": params, "references": references}

    # /def

    @classmethod
    def set(
        cls,
        value: dict,
        register_as: T.Optional[str] = None,
        references: T.Optional[dict] = None,
    ):
        """Set (see ScienceState) with optional registering.

        Parameters
        ----------
        value : dict
        register_as : str, optional
            the name of the science state to set
        references : dict, optional
            references for `value`. Only used if `register_as` is str.

        """
        super().set(value)

        if isinstance(register_as, str):
            cls._registry[register_as] = {
                "params": value,
                "references": references or {},
            }

    # /def


# /class


# -------------------------------------------------------------------


vesc_sun_at_earth = 42.1 * _KMS
"""Escape velocity from the sun, starting at 1 AU [1]_.

References
----------
.. [1] https://en.wikipedia.org/wiki/Escape_velocity.

"""


# -------------------------------------------------------------------


@u.quantity_input(R="length")
def vesc_sun_at_R(R):
    r"""Escape velocity from the sun, starting at position R.

    The Newtonian escape velocity from a spherical mass distribution,
    starting at a position r_0, is :math:`v_0=\sqrt{2gr_0}`. If this value
    is known at some position r0, than it is known at all R by virtue
    of ratios. We calculate this value for the sun relative to the
    known value at the earth -- 42.1 km / s.

    Parameters
    ----------
    R : Distance
        from the sun.

    Returns
    -------
    vesc : Quantity

    """
    ratio: float = R.to_value(u.AU)  # b/c r_earth = 1 AU
    return np.sqrt(ratio) * vesc_sun_at_earth


# /def

# -------------------------------------------------------------------


_multibody_escape_wikipedia = r"""
    When escaping a compound system, such as a moon orbiting a planet or a
    planet orbiting a sun, a rocket that leaves at escape velocity
    (:math:`ve_1`) for the first (orbiting) body, (e.g. Earth) will not travel
    to an infinite distance because it needs an even higher speed to escape
    gravity of the second body (e.g. the Sun). Near the Earth, the rocket's
    trajectory will appear parabolic, but it will still be gravitationally
    bound to the second body and will enter an elliptical orbit around that
    body, with an orbital speed similar to the first body.

    To escape the gravity of the second body once it has escaped the first
    body the rocket will need to be traveling at the escape velocity for the
    second body (:math:`ve_2`) (at the orbital distance of the first body).
    However, when the rocket escapes the first body it will still have the
    same orbital speed around the second body that the first body has (vo). So
    its excess velocity as it escapes the first body will need to be the
    difference between the orbital velocity and the escape velocity. With a
    circular orbit, escape velocity is sqrt{2} times the orbital speed.
    Thus the total escape velocity vte when leaving one body orbiting a second
    and seeking to escape them both is, under simplified assumptions:

    .. math::

        v_{te}=\sqrt{(v_{e2} - v_o)^2 + v_{e1}^2}
        = \sqrt{\left(k v_{e2}\right)^2 + v_{e1}^2}

    where :math:`k=1−1/\sqrt{2} \sim 0.2929` for circular orbits.
"""  # TODO instead indent by function


@u.quantity_input(ve1="speed", ve2="speed")
@format_doc(None, wikipedia=_multibody_escape_wikipedia)
def twobody_vesc(
    ve1: u.Quantity,
    ve2: u.Quantity,
    vo: T.Union[None, T.Sequence[u.Quantity]] = None,
):
    r"""Two-body escape velocity.

    Parameters
    ----------
    ve1, ve2: :class:`~astropy.units.Quantity`
        Escape velocities.
    vo : Quantity or None, optional
        The orbital velocity of object 1 around object 2.

    Returns
    -------
    vesc : :class:`~astropy.units.Quantity`
        The compound escape velocity.

    Examples
    --------
    For a Galactic macro falling into the potential well of the Earth,
    the minimum observed velocity

        >>> vesc_earth = 11.186 * u.km / u.s
        >>> vesc_sun_at_earth = 42.1 * u.km / u.s
        >>> twobody_vesc(vesc_earth, vesc_sun_at_earth)
        <Quantity 16.6485836 km / s>

    Notes
    -----
    From `Wikipedia <https://en.wikipedia.org/wiki/Escape_velocity>`_:

    {wikipedia}

    """
    vo = vo or ve2 / _sqrt2  # None -> circular

    return np.sqrt(ve1 ** 2 + (ve2 - vo) ** 2)


# /def


def _norm_v1_v2(v1: T.Sequence, v2: T.Sequence) -> T.Sequence:
    return np.sqrt(v1 ** 2.0 + v2 ** 2.0)


# /def


@format_doc(None, wikipedia=_multibody_escape_wikipedia)
class multibody_vesc:
    """Multi-body escape velocity.

    Parameters
    ----------
    *vescs: Quantity
        velocities, ordered from 1st to last body.

    vo : list of Quantity or None, optional
        The orbital velocity of object `vescs`(i+1) around object `vescs`(i).
        if list of quantities, must match `vescs` in length
        if None (default) then orbits are assumed circular.

    accumulate : bool
        whether to return the accumulative escape velocity for each larger
        system (True), or just the total escape velocity (False, default).

    Returns
    -------
    :class:`~astropy.units.Quantity`
        The compound escape velocity
        if `accumulate` False (default) then scalar, else accumulated vector.

    Examples
    --------
    For a Galactic macro falling into the potential well of the Earth,
    the minimum observed velocity

        >>> vesc_earth = 11.186 * u.km / u.s
        >>> vesc_sun_at_earth = 42.1 * u.km / u.s
        >>> vesc_gal_at_sun = 550 * u.km / u.s
        >>> multibody_vesc(vesc_earth, vesc_sun_at_earth, vesc_gal_at_sun)
        <Quantity 161.94929058 km / s>

        >>> multibody_vesc.accumulate(vesc_earth, vesc_sun_at_earth, vesc_gal_at_sun)
        <Quantity [ 11.186     ,  16.6485836 , 161.94929058] km / s>

    Notes
    -----
    From `Wikipedia <https://en.wikipedia.org/wiki/Escape_velocity>`_:

    {wikipedia}

    """

    @staticmethod
    def prepare_vels(
        *vescs: QuantityType, vo: QuantityType = None
    ) -> QuantityType:
        r"""Constructs effective velocities.

        In a two-body escape,

        .. math::

            v_{te}=\sqrt{(v_{e2} - v_o)^2 + v_{e1}^2}

        This can be rewritten as

        .. math::

            v_{te}=\sqrt{v_{e2}^{eff}^2 + v_{e1}^{eff}^2}

        where :math:`v_{e1}^{eff} = v_{e1}` and all higher terms
        have the orbital velocity subtracted.

        Returns
        -------
        vs : :class:`~astropy.units.Quantity`
            same len as vescs

        """
        vs: u.Quantity = as_quantity(vescs)

        if vo is None:
            vs[1:] = vs[1:] * (1 - 1 / np.sqrt(2))
        else:
            vs[1:] = vs[1:] - vo

        return vs

    # /def

    # instead of call
    def __new__(
        cls,
        *vescs: QuantityType,
        vo: T.Union[None, T.Sequence[QuantityType]] = None,
    ) -> QuantityType:
        """Evaluate multi-body escape velocity, finding final velocity."""
        vs = cls.prepare_vels(*vescs, vo=vo)
        return functools.reduce(_norm_v1_v2, vs)

    # /def

    @classmethod
    def accumulate(
        cls,
        *vescs: QuantityType,
        vo: T.Union[None, T.Sequence[QuantityType]] = None,
    ) -> QuantityType:
        """Evaluate multi-body escape velocity, accumulated."""
        vs = cls.prepare_vels(*vescs, vo=vo)
        return as_quantity(tuple(itertools.accumulate(vs, _norm_v1_v2)))

    # /def


# /class


##############################################################################
# END

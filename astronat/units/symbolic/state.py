# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   :
# AUTHOR  :
# PROJECT :
#
# ----------------------------------------------------------------------------

"""**DOCSTRING**.

Description.

"""

__author__ = ""
# __copyright__ = "Copyright 2019, "
# __credits__ = [""]
# __license__ = ""
# __version__ = "0.0.0"
# __maintainer__ = ""
# __email__ = ""
# __status__ = "Production"


# __all__ = [
#     # functions
#     "",
#     # other
#     "",
# ]


##############################################################################
# IMPORTS

# BUILT-IN

import typing as T

from astropy.utils.state import ScienceState
from sympy.physics.units import DimensionSystem, UnitSystem
from sympy.physics.units.systems import SI

from .astro import dimsys_astro, galactic
from .base import dimsys_base  # , unit_base

# THIRD PARTY


# PROJECT-SPECIFIC


##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


class dimension_system(ScienceState):
    """Sympy Dimension System."""

    _latest_value: str = "default"
    _references: T.Optional[T.Dict[str, T.Any]] = None
    _value: T.Optional[DimensionSystem] = None

    _registry: T.Dict[str, T.Dict[str, T.Dict[str, T.Any]]] = {
        "DEFAULT": {"dimsys": dimsys_astro, "references": {}},
        "base": {"dimsys": dimsys_base, "references": {}},
    }

    @classmethod
    def get_dimsys_from_string(cls, arg: str):
        """Get parameters from registry."""
        # Resolve the meaning of 'latest'
        #         if arg == "latest":
        #             arg = cls._latest_value

        if arg.lower() == "default":
            info = cls._registry["DEFAULT"]

        elif arg in cls._registry:
            info = cls._registry[arg]

        else:
            raise ValueError("Invalid string input.")

        return info["dimsys"], info["references"]

    # /def

    @classmethod
    def validate(cls, value: T.Union[None, str, DimensionSystem] = None):
        """Validate `value`, from string or DimensionSystem."""
        if value is None:
            value = "default"

        if isinstance(value, str):
            params, refs = cls.get_dimsys_from_string(value)
            cls._references = refs
            return params

        elif isinstance(value, DimensionSystem):
            return value

        else:
            raise ValueError(
                "Invalid input; must be a string, or DimensionSystem."
            )

    # /def

    @classmethod
    def register(cls, name: str, dimsys: DimensionSystem, references: dict):
        """Register a set of parameters.

        Parameters
        ----------
        name : str
        dimsys : :class:`~sympy.physics.units.DimensionSystem`
        references : dict

        """
        cls._registry[name] = {"dimsys": dimsys, "references": references}

    # /def

    @classmethod
    def set(
        cls,
        value: T.Union[str, DimensionSystem],
        register_as: T.Optional[str] = None,
        references: T.Optional[dict] = None,
    ):
        """Set (see ScienceState) with optional registering.

        Parameters
        ----------
        value : DimensionSystem
        register_as : str, optional
            the name of the science state to set
        references : dict, optional
            references for `value`. Only used if `register_as` is str.

        """
        super().set(value)

        if isinstance(register_as, str):
            cls._registry[register_as] = {
                "dimsys": value,
                "references": references or {},
            }

    # /def


# /class


class unit_system(ScienceState):
    """Sympy Dimension System."""

    _latest_value: str = "SI"
    _references: T.Optional[T.Dict[str, T.Any]] = None
    _value: T.Optional[UnitSystem] = None

    _registry: T.Dict[str, T.Dict[str, T.Dict[str, T.Any]]] = {
        "SI": {"unitsys": SI, "references": {}},
        "galactic": {"unitsys": galactic, "references": {}},
    }

    @classmethod
    def get_unitsys_from_string(cls, arg: str):
        """Get parameters from registry."""
        if arg.lower() == "default":
            info = cls._registry["SI"]

        elif arg in cls._registry:
            info = cls._registry[arg]

        else:
            raise ValueError("Invalid string input.")

        return info["unitsys"], info["references"]

    # /def

    @classmethod
    def validate(cls, value: T.Union[None, str, UnitSystem] = None):
        """Validate `value`, from string or UnitSystem."""
        if value is None:
            value = "default"

        if isinstance(value, str):
            params, refs = cls.get_unitsys_from_string(value)
            cls._references = refs
            return params

        elif isinstance(value, UnitSystem):
            return value

        else:
            raise ValueError("Invalid input; must be a string or UnitSystem.")

    # /def

    @classmethod
    def register(cls, name: str, unitsys: UnitSystem, references: dict):
        """Register a set of parameters.

        Parameters
        ----------
        name : str
        unitsys : :class:`~sympy.physics.units.UnitSystem`
        references : dict

        """
        cls._registry[name] = {"unitsys": unitsys, "references": references}

    # /def

    @classmethod
    def set(
        cls,
        value: UnitSystem,
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
                "unitsys": value,
                "references": references or {},
            }

    # /def


# /class


##############################################################################
# END

# -*- coding: utf-8 -*-

"""Units Definition."""

__all__ = [
    "kyr",
    "kiloyear",
    "Myr",
    "megayear",
    "pc",
    "parsec",
    "kpc",
    "kiloparsec",
    "Mpc",
    "megaparsec",
    "solMass",
    "MSun",
]


##############################################################################
# IMPORTS

# BUILT-IN

# THIRD PARTY

# PROJECT-SPECIFIC

from sympy import pi
from sympy.physics.units.definitions import unit_definitions
from sympy.physics.units.definitions.unit_definitions import *
from sympy.physics.units.definitions.unit_definitions import (
    astronomical_unit,
    year,
)
from sympy.physics.units.prefixes import kilo, mega
from sympy.physics.units.quantities import Quantity

__all__ += (
    unit_definitions.__all__ if hasattr(unit_definitions, "__all__") else []
)

##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################

yr = Quantity("year", abbrev="yr")
yr.set_global_relative_scale_factor(1.0, year)

kyr = kiloyear = Quantity("kiloyear", abbrev="kyr")
kyr.set_global_relative_scale_factor(kilo, year)

Myr = megayear = Quantity("megayear", abbrev="Myr")
Myr.set_global_relative_scale_factor(mega, year)

# -------------------------------------------------------------------

pc = parsec = Quantity("parsec", abbrev="pc")
pc.set_global_relative_scale_factor(180 * 3600 / pi, astronomical_unit)

kpc = kiloparsec = Quantity("kiloparsec", abbrev="kpc")
kpc.set_global_relative_scale_factor(kilo, parsec)

Mpc = megaparsec = Quantity("megaparsec", abbrev="Mpc")
Mpc.set_global_relative_scale_factor(mega, parsec)

# -------------------------------------------------------------------

solMass = MSun = Quantity("solar_mass", abbrev="MSun")

##############################################################################
# END

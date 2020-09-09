# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Symbolic Units, built on :mod:`~sympy.physics.units`."""

# __all__ = [
#     # modules
#     "",
#     # functions
#     "",
#     # other
#     "",
# ]


##############################################################################
# IMPORTS

# BUILT IN

# THIRD PARTY

from sympy import S
from sympy.physics.units import UnitSystem
from sympy.physics.units.prefixes import PREFIXES, prefix_unit

# from .units_definition import kpc, Myr, solMass, rad
from . import unit_definitions as u
from .base import dimsys_base

# PROJECT-SPECIFIC


##############################################################################
# PARAMETERS

One = S.One


# Prefixes of units like g, J, N etc get added using `prefix_unit`
# in the for loop, but the actual units have to be added manually.
units = []
units.extend([u.m, u.g, u.s, u.J, u.N, u.W, u.Pa, u.Hz])  # MKS
units.extend([u.A, u.V, u.ohm, u.S, u.F, u.H, u.C, u.T, u.Wb])  # MKSA
units.extend(
    [
        u.mol,
        u.cd,
        u.K,
        u.lux,
        u.hertz,
        u.newton,
        u.pascal,
        u.joule,
        u.watt,
        u.coulomb,
        u.volt,
        u.farad,
        u.ohm,
        u.siemens,
        u.weber,
        u.tesla,
        u.henry,
        u.candela,
        u.lux,
        u.becquerel,
        u.gray,
        u.katal,
    ]
)  # SI
units.extend([u.pc, u.yr])  # astronomical

all_units = units[:]

for unit in units:
    all_units.extend(prefix_unit(unit, PREFIXES))

all_units.extend([u.G, u.c])  # MKS
all_units.extend([u.Z0])  # MKSA
all_units.extend([u.mol, u.cd, u.K, u.lux])  # SI
all_units.extend([u.Myr])  # astronomical


##############################################################################
# CODE
##############################################################################

dimsys_astro = dimsys_base.extend(
    [],
    new_derived_dims=dict(
        # cosmological
        hubble=dict(time=-1)
    ),
)

# -------------------------------------------------------------------

galactic = UnitSystem(
    base_units=(u.kpc, u.Myr, u.solMass, u.rad),
    units=all_units,
    name="Galactic",
    dimension_system=dimsys_astro,
)

# galactic = unit_base.extend(
#     [u.kpc, u.Myr, u.solMass, u.rad],
#     name="Galactic",
#     dimension_system=dimsys_astro,
# )

# galactic.set_quantity_scale_factor(u.kpc, One)
# galactic.set_quantity_scale_factor(u.Myr, One)
# galactic.set_quantity_scale_factor(u.rad, One)

# galactic.set_quantity_scale_factor(u.MSun, One)
galactic.set_quantity_scale_factor(u.kilogram, 1 / 1.9891e30 * u.MSun)
galactic.set_quantity_scale_factor(u.MSun, 1.9891e30 * u.kilogram)


# -------------------------------------------------------------------

# dims = (
#     velocity,
#     acceleration,
#     momentum,
#     force,
#     energy,
#     power,
#     pressure,
#     frequency,
#     action,
# )

# units = [m, g, s, J, N, W, Pa, Hz]
# all_units = []


# # Prefixes of units like g, J, N etc get added using `prefix_unit`
# # in the for loop, but the actual units have to be added manually.
# all_units.extend([g, J, N, W, Pa, Hz])

# for u in units:
#     all_units.extend(prefix_unit(u, PREFIXES))
# all_units.extend([G, c])


##############################################################################
# END

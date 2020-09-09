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

# THIRD PARTY

from sympy import Rational, S, pi, sqrt
from sympy.physics.units import DimensionSystem, UnitSystem
from sympy.physics.units.definitions import dimension_definitions as d

from . import unit_definitions as u
from .utils import _set_quantity

# from sympy.physics.units.definitions import unit_definitions as u


# PROJECT-SPECIFIC


###############################################################################
# Parameters

One = S.One


###############################################################################
# DIMENSION SYSTEM
###############################################################################

dimsys_base = DimensionSystem(
    # Dimensional dependencies
    [
        d.length,
        d.mass,
        d.time,
        d.current,
        d.temperature,
        d.amount_of_substance,
        d.luminous_intensity,
        d.information,
    ],
    # Dimensional dependencies for derived dimensions
    dimensional_dependencies=dict(
        velocity=dict(length=1, time=-1),
        acceleration=dict(length=1, time=-2),
        momentum=dict(mass=1, length=1, time=-1),
        force=dict(mass=1, length=1, time=-2),
        energy=dict(mass=1, length=2, time=-2),
        power=dict(length=2, mass=1, time=-3),
        pressure=dict(mass=1, length=-1, time=-2),
        frequency=dict(time=-1),
        action=dict(length=2, mass=1, time=-1),
        volume=dict(length=3),
        voltage=dict(mass=1, length=2, current=-1, time=-3),
        impedance=dict(mass=1, length=2, current=-2, time=-3),
        conductance=dict(mass=-1, length=-2, current=2, time=3),
        capacitance=dict(mass=-1, length=-2, current=2, time=4),
        inductance=dict(mass=1, length=2, current=-2, time=-2),
        charge=dict(current=1, time=1),
        magnetic_density=dict(mass=1, current=-1, time=-2),
        magnetic_flux=dict(length=2, mass=1, current=-1, time=-2),
    ),
)


# -------------------------------------------------------------------
# Assign dimensions

# -----------------------------------------------
# Base Units

# length
_set_quantity(u.meter, d.length, One, system=dimsys_base)

# mass
_set_quantity(u.gram, d.mass, One, system=dimsys_base)

# time
_set_quantity(u.second, d.time, One, system=dimsys_base)

# current
_set_quantity(u.ampere, d.current, One, system=dimsys_base)

# dimsys_base.set_quantity_dimension(u.statampere, d.current)
# set_quantity_scale_factor(statampere, u.statcoulomb/u.second)  # TODO


# -----------------------------------------------
# Derived Units

# volume
_set_quantity(u.liter, d.length ** 3, u.meter ** 3 / 1000, system=dimsys_base)

# force
_set_quantity(
    u.newton, d.force, u.kilogram * u.meter / u.second ** 2, system=dimsys_base
)

# energy
_set_quantity(u.joule, d.energy, u.newton * u.meter, system=dimsys_base)

# power
_set_quantity(u.watt, d.power, u.joule / u.second, system=dimsys_base)

# pressure
_set_quantity(
    u.pascal, d.pressure, u.newton / u.meter ** 2, system=dimsys_base
)

# charge
# coulomb done in unit_definitions
_set_quantity(
    u.statcoulomb, d.charge, 5 * u.coulomb / 149896229, system=dimsys_base
)

# voltage
# volt done in unit_definitions
_set_quantity(u.statvolt, d.voltage, u.erg / u.statcoulomb, system=dimsys_base)

# magnetic_density
# dimsys_base.set_quantity_dimension(u.gauss, d.magnetic_density)
# dimsys_base.set_quantity_dimension(u.tesla, d.magnetic_density)

# magnetic_flux
# dimsys_base.set_quantity_dimension(u.maxwell, d.magnetic_flux)

# luminous density
# dimsys_base.set_quantity_dimension(u.lux, d.luminous_intensity / d.length ** 2)

# rate
# dimsys_base.set_quantity_dimension(u.katal, d.amount_of_substance / d.time)

# gray is the SI unit of absorbed dose
# dimsys_base.set_quantity_dimension(u.gray, d.energy / d.mass)

# u.becquerel is the SI unit of radioactivity
# dimsys_base.set_quantity_dimension(u.becquerel, 1 / d.time)


# INVERSE UNITS

_set_quantity(u.dioptre, 1 / d.length, 1 / u.meter, system=dimsys_base)

_set_quantity(u.hertz, d.frequency, One, system=dimsys_base)


#####################################################################
# Constants

# Newton constant
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.gravitational_constant,
    d.length ** 3 * d.mass ** -1 * d.time ** -2,
    6.67430e-11 * u.meter ** 3 / (u.kilogram * u.second ** 2),
    system=dimsys_base,
)

# speed of light

_set_quantity(
    u.speed_of_light,
    d.velocity,
    299792458 * u.meter / u.second,
    system=dimsys_base,
)

# Planck constant
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.planck, d.action, 6.62607015e-34 * u.joule * u.second, system=dimsys_base
)

# Reduced Planck constant
# REF: NIST SP 959 (June 2019)

_set_quantity(u.hbar, d.action, u.planck / (2 * pi), system=dimsys_base)

# elementary charge
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.elementary_charge,
    d.charge,
    1.602176634e-19 * u.coulomb,
    system=dimsys_base,
)

# Electronvolt
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.electronvolt, d.energy, 1.602176634e-19 * u.joule, system=dimsys_base
)

# Avogadro number
# REF: NIST SP 959 (June 2019)

_set_quantity(u.avogadro_number, One, 6.02214076e23, system=dimsys_base)

# Avogadro constant

_set_quantity(
    u.avogadro_constant,
    d.amount_of_substance ** -1,
    u.avogadro_number / u.mol,
    system=dimsys_base,
)

# Boltzmann constant
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.boltzmann_constant,
    d.energy / d.temperature,
    1.380649e-23 * u.joule / u.kelvin,
    system=dimsys_base,
)

# Stefan-Boltzmann constant
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.stefan_boltzmann_constant,
    d.energy * d.time ** -1 * d.length ** -2 * d.temperature ** -4,
    pi ** 2
    * u.boltzmann_constant ** 4
    / (60 * u.hbar ** 3 * u.speed_of_light ** 2),
    system=dimsys_base,
)

# Atomic mass
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.atomic_mass_constant,
    d.mass,
    1.66053906660e-24 * u.gram,
    system=dimsys_base,
)

# Molar gas constant
# REF: NIST SP 959 (June 2019)

_set_quantity(
    u.molar_gas_constant,
    d.energy / (d.temperature * d.amount_of_substance),
    u.boltzmann_constant * u.avogadro_constant,
    system=dimsys_base,
)

# Faraday constant

_set_quantity(
    u.faraday_constant,
    d.charge / d.amount_of_substance,
    u.elementary_charge * u.avogadro_constant,
    system=dimsys_base,
)

# Josephson constant

_set_quantity(
    u.josephson_constant,
    d.frequency / d.voltage,
    0.5 * u.planck / u.elementary_charge,
    system=dimsys_base,
)

# Von Klitzing constant

_set_quantity(
    u.von_klitzing_constant,
    d.voltage / d.current,
    u.hbar / u.elementary_charge ** 2,
    system=dimsys_base,
)

# Acceleration due to gravity (on the Earth surface)

_set_quantity(
    u.acceleration_due_to_gravity,
    d.acceleration,
    9.80665 * u.meter / u.second ** 2,
    system=dimsys_base,
)

# magnetic constant:

_set_quantity(
    u.magnetic_constant,
    d.force / d.current ** 2,
    4 * pi / 10 ** 7 * u.newton / u.ampere ** 2,
    system=dimsys_base,
)

# electric constant:

_set_quantity(
    u.vacuum_permittivity,
    u.capacitance / d.length,
    1 / (u.u0 * u.c ** 2),
    system=dimsys_base,
)

# vacuum impedance:

_set_quantity(u.vacuum_impedance, u.impedance, u.u0 * u.c, system=dimsys_base)

# Coulomb's constant:
_set_quantity(
    u.coulomb_constant,
    d.force * d.length ** 2 / d.charge ** 2,
    1 / (4 * pi * u.vacuum_permittivity),
    system=dimsys_base,
)

# psi
_set_quantity(
    u.psi, d.pressure, u.pound * u.gee / u.inch ** 2, system=dimsys_base
)

_set_quantity(
    u.mmHg,
    d.pressure,
    u.dHg0 * u.acceleration_due_to_gravity * u.kilogram / u.meter ** 2,
    system=dimsys_base,
)

_set_quantity(
    u.milli_mass_unit, d.mass, u.atomic_mass_unit / 1000, system=dimsys_base
)

_set_quantity(
    u.quart, d.length ** 3, Rational(231, 4) * u.inch ** 3, system=dimsys_base
)

# Other convenient units and magnitudes

_set_quantity(
    u.lightyear, d.length, u.speed_of_light * u.julian_year, system=dimsys_base
)

_set_quantity(
    u.astronomical_unit, d.length, 149597870691 * u.meter, system=dimsys_base
)

# Fundamental Planck units:

_set_quantity(
    u.planck_mass,
    d.mass,
    sqrt(u.hbar * u.speed_of_light / u.G),
    system=dimsys_base,
)

_set_quantity(
    u.planck_time,
    d.time,
    sqrt(u.hbar * u.G / u.speed_of_light ** 5),
    system=dimsys_base,
)

_set_quantity(
    u.planck_temperature,
    d.temperature,
    sqrt(u.hbar * u.speed_of_light ** 5 / u.G / u.boltzmann ** 2),
    system=dimsys_base,
)

_set_quantity(
    u.planck_length,
    d.length,
    sqrt(u.hbar * u.G / u.speed_of_light ** 3),
    system=dimsys_base,
)

_set_quantity(
    u.planck_charge,
    d.charge,
    sqrt(4 * pi * u.electric_constant * u.hbar * u.speed_of_light),
    system=dimsys_base,
)

# Derived Planck units:

_set_quantity(
    u.planck_area, d.length ** 2, u.planck_length ** 2, system=dimsys_base
)

_set_quantity(
    u.planck_volume, d.length ** 3, u.planck_length ** 3, system=dimsys_base
)

_set_quantity(
    u.planck_momentum,
    d.mass * d.velocity,
    u.planck_mass * u.speed_of_light,
    system=dimsys_base,
)

_set_quantity(
    u.planck_energy,
    d.energy,
    u.planck_mass * u.speed_of_light ** 2,
    system=dimsys_base,
)

_set_quantity(
    u.planck_force,
    d.force,
    u.planck_energy / u.planck_length,
    system=dimsys_base,
)

_set_quantity(
    u.planck_power,
    d.power,
    u.planck_energy / u.planck_time,
    system=dimsys_base,
)

_set_quantity(
    u.planck_density,
    d.mass / d.length ** 3,
    u.planck_mass / u.planck_length ** 3,
    system=dimsys_base,
)

_set_quantity(
    u.planck_energy_density,
    d.energy / d.length ** 3,
    u.planck_energy / u.planck_length ** 3,
    system=dimsys_base,
)

_set_quantity(
    u.planck_intensity,
    d.mass * d.time ** (-3),
    u.planck_energy_density * u.speed_of_light,
    system=dimsys_base,
)

_set_quantity(
    u.planck_angular_frequency,
    1 / d.time,
    1 / u.planck_time,
    system=dimsys_base,
)

_set_quantity(
    u.planck_pressure,
    d.pressure,
    u.planck_force / u.planck_length ** 2,
    system=dimsys_base,
)

_set_quantity(
    u.planck_current,
    d.current,
    u.planck_charge / u.planck_time,
    system=dimsys_base,
)

_set_quantity(
    u.planck_voltage,
    d.voltage,
    u.planck_energy / u.planck_charge,
    system=dimsys_base,
)

_set_quantity(
    u.planck_impedance,
    u.impedance,
    u.planck_voltage / u.planck_current,
    system=dimsys_base,
)

_set_quantity(
    u.planck_acceleration,
    d.acceleration,
    u.speed_of_light / u.planck_time,
    system=dimsys_base,
)

# Older units for radioactivity

_set_quantity(
    u.curie, 1 / d.time, 37000000000 * u.becquerel, system=dimsys_base
)

_set_quantity(
    u.rutherford, 1 / d.time, 1000000 * u.becquerel, system=dimsys_base
)


#####################################################################
# Unit System
#####################################################################

# # unit system
# unit_base = UnitSystem(
#     base_units=(),
#     units=[],
#     name="empty_unit_base",
#     dimension_system=dimsys_base,
# )

##############################################################################
# END

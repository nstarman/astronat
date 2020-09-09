# -*- coding: utf-8 -*-

"""Test :mod:`~astronat.dynamics.vesc`."""


__all__ = [
    "test_solar_system_vesc_params",
    "test_vesc_sun_at_R",
]


##############################################################################
# IMPORTS

import astropy.units as u
import pytest

from .. import vesc

##############################################################################
# PARAMETERS

_KMS = u.km / u.s


##############################################################################
# CODE
##############################################################################


def test_solar_system_vesc_params():
    """Test :class:`~macro_lightning.parameters.solar_system_vesc_params`."""
    # test set
    vesc.solar_system_vesc_params.set("latest")

    # test setting set references
    refs = vesc.solar_system_vesc_params._references
    expected = vesc.solar_system_vesc_params._registry["DEFAULT"]["references"]
    assert refs == expected

    # test get
    ss = vesc.solar_system_vesc_params.get()
    expected = vesc.solar_system_vesc_params._registry["DEFAULT"]["params"]
    assert ss == expected

    # test register
    new_params = {k: v + 2 * _KMS for k, v in ss.items()}
    vesc.solar_system_vesc_params.register(
        "new", params=new_params, references={"_source": None},
    )

    # test can set
    vesc.solar_system_vesc_params.set("new")

    # test getting works
    new_ss = vesc.solar_system_vesc_params.get()
    assert new_ss == new_params

    pass


# /def


# -------------------------------------------------------------------


def test_vesc_sun_at_R():
    """Test :func:`~macro_lightning.parameters.vesc_sun_at_R`."""
    # Unity test
    assert vesc.vesc_sun_at_R(1 * u.AU) == vesc.vesc_sun_at_earth

    # quadruple distance = half the escape velocity
    assert vesc.vesc_sun_at_R(4 * u.AU) == 2.0 * vesc.vesc_sun_at_earth

    # Error test
    with pytest.raises(Exception):
        vesc.vesc_sun_at_R(1)

    with pytest.raises(Exception):
        vesc.vesc_sun_at_R(1 * u.deg)


# /def


# -------------------------------------------------------------------

# def test_twobody_vesc


##############################################################################
# END

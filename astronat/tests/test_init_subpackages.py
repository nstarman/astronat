# -*- coding: utf-8 -*-

"""Tests for :mod:`~utilipy.utils`."""


__all__ = [
    "test_init_constants",
    "test_constants_top_level_imports",
    "test_init_units",
    "test_init_utils",
    "test_utils_top_level_imports",
]


##############################################################################
# IMPORTS

# BUILT-IN

import os


# PROJECT-SPECIFIC

from .. import (
    constants,
    units,
    utils,
)


##############################################################################
# TESTS
##############################################################################


def test_init_constants():
    """Test :mod:`~utilipy.constants` initialization."""
    # Expectations
    local = [
        "conf",
        "frozen",
        "FrozenConstants",
        "default_values",
        "ConstantsValues",
        # top-level
        "values",
        "_frozen",
    ]

    # test __all__ conforms to module
    for name in constants.__all__:
        assert hasattr(constants, name)

    # test __all__ matches expectations
    for name in constants.__all__:
        assert name in local

    return


# /def


def test_constants_top_level_imports():
    """Test top-level imports of :mod:`~utilipy.constants`."""
    # First test they exist
    subpkg: str
    for subpkg in constants.__all_top_imports__:
        assert hasattr(constants, subpkg)

    # Next test that top-levels are all the possible top-levels
    drct: str = os.path.split(constants.__file__)[0]  # directory
    donottest = ("__pycache__", "data", "tests")  # stuff not to test

    for file in os.listdir(drct):  # iterate through directory
        # test?
        if os.path.isdir(drct + "/" + file) and file not in donottest:
            assert file in constants.__all_top_imports__
        else:  # nope, chuck testa.
            pass

    return


# /def


# --------------------------------------------------------------------------


def test_init_units():
    """Test :mod:`~utilipy.units` initialization."""
    # Expectations
    local = [
        # top-level
        "amuse",
        "composite",
        "core",
        "decorators",
        "full_amuse",
        # core
        "quantity_return_",
        "ExpandedUnitType",
        # decorators
        "quantity_output",
        "quantity_io",
        "QuantityInputOutput",
    ]

    if hasattr(units, "__all__"):

        # test __all__ conforms to module
        for name in units.__all__:
            assert hasattr(units, name)

        # test __all__ matches expectations
        for name in units.__all__:
            assert name in local

    elif local:  # has locals, but no __all__. must fail.
        raise ValueError

    return


# /def


# --------------------------------------------------------------------------


def test_init_utils():
    """Test :mod:`~utilipy.utils` initialization."""
    # Expectations
    local = [
        "typing",
        # top-level imports
        # "logging",
    ]

    # test __all__ conforms to module
    for name in utils.__all__:
        assert hasattr(utils, name)

    # test __all__ matches expectations
    for name in utils.__all__:
        assert name in local

    return


# /def


def test_utils_top_level_imports():
    """Test Top-Level Imports."""
    # First test they exist
    subpkg: str
    for subpkg in utils.__all_top_imports__:
        assert hasattr(utils, subpkg)

    # Next test that top-levels are all the possible top-levels
    drct: str = os.path.split(utils.__file__)[0]  # directory
    donottest = ("tests", "__pycache__")  # stuff not to test

    for file in os.listdir(drct):  # iterate through directory
        # test?
        if os.path.isdir(drct + "/" + file) and file not in donottest:
            assert file in utils.__all_top_imports__
        else:  # nope, chuck testa.
            pass

    return


# /def


# --------------------------------------------------------------------------


##############################################################################
# END

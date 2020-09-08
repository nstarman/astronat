# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Test Unit Module."""


__all__ = [
    # core and decorators
    "test_core",
    "test_decorators",
    # added units
    "test_amuse",
    "test_composite",
    "test_full_amuse",
]


##############################################################################
# IMPORTS

from . import test_composite  # core and decorators; added units
from . import test_amuse, test_core, test_decorators, test_full_amuse

##############################################################################
# END

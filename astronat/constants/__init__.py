# -*- coding: utf-8 -*-

"""Astropy Constants. Extended.

Astropy constants, with a frozen version for reproducibility.

float versions of the constants accessible through `values` module
this includes frozen version for reproducibility
to access frozen version, set ``frozen_constants=True`` in utilipy config


References
----------
References [#]_.

.. [#] Astropy Collaboration et al., 2018, AJ, 156, 123.


Notes
-----
.. todo::

    - use UnitSystem to choose how to represent the units, then values.
    - worry about when do stuff like ``import imperial``.

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "frozen",
    "FrozenConstants",
    "default_values",
    "ConstantsValues",
    "conf",
]


##############################################################################
# IMPORTS

from astropy import constants
from astropy.constants import *  # noqa

from ._frozen import FrozenConstants, frozen
from .setup_package import conf
from .values import ConstantsValues
from .values import values as default_values

# -------------------------------------------------------------------
# __ALL__

__all_top_imports__ = ("values", "_frozen")

__all__ += list(__all_top_imports__)
# __all__ += constants.__all__


#############################################################################
# END

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

# THIRD PARTY


# PROJECT-SPECIFIC


##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


def _set_quantity(q, dim, scale, system=None):

    system.set_quantity_dimension(q, dim)
    system.set_quantity_scale_factor(q, scale)

    return q


# /def


def set_quantity(q, dim, scale, unit_system=None):
    from .state import unit_system as unit_system_state

    if unit_system is None:
        unit_system = unit_system_state.get()

    return _set_quantity(q, dim, scale, system=unit_system)


# /def


# -------------------------------------------------------------------


##############################################################################
# END

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Setup for :mod:`~utilipy.constants`."""


###############################################################################
# IMPORTS

from __future__ import absolute_import

from astropy import config as _config

###############################################################################
# CODE
###############################################################################


class Conf(_config.ConfigNamespace):
    """Configuration parameters for `utilipy.constants`."""

    frozen_constants = _config.ConfigItem(
        True,
        description=(
            "constants set by data file (True), "
            "or by .value on astropy constants at runtime (False)"
        ),
        cfgtype="boolean(default=True)",
    )


conf = Conf()
# /class


###############################################################################
# INFO


__all__ = [
    "conf",  # configuration instance
]


###############################################################################
# END

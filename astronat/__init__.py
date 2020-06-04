# -*- coding: utf-8 -*-
# see LICENSE.rst

"""The `astronat` Package.

Welcome to astronat, for astrophysics  concomitant codes.
There are modules for astronomical functions, improving astropy units
and quantity-enabled functions, and much more.

Routine Listings
----------------
`help`
    *astronat* help function. Online search or offline overview.

`online_help`
    Search the online *astronat* documentation for the given query.

References
----------
Utilipy [1]_, Astropy [2]_

.. [1] nstarman. (2020, March 23). nstarman/utilipy: astropy_template (Version
    astropy_template). Zenodo. http://doi.org/10.5281/zenodo.3724822
.. [2] Astropy Collaboration et al., 2018, AJ, 156, 123.


Examples
--------

help
^^^^

To do a specific search of `astronat`'s docs, use the ``online_help`` function.
If you don't have a specific query, that's fine too,
`astronat` will open the general search page.

As an example, here we query RTD for the documentation on `LogFile`.

    >>> import astronat
    >>> astronat.online_help(query="distance_modulus") # doctest: +SKIP

The same can be accomplished with the general `help` function.

    >>> import astronat
    >>> astronat.help(query="distance_modulus", online=True) # doctest: +SKIP

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    # modules
    "constants",
    "dynamics",
    "units",
    "utils",
    # functions
    "help",
    "online_help",
    "reload_config",
]

__all_top_imports__ = (
    "constants",
    "dynamics",
    # "extern",
    # "phot",
    # "sc",
    "units",
    "utils",
)

##############################################################################
# IMPORTS

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# (sets the __version__)
from ._astropy_init import *  # noqa
from ._astropy_init import __version__  # noqa

# CUSTOM

import utilipy
from utilipy.utils import typing as T


# THIRD PARTY

import astropy.config as config


# PROJECT-SPECIFIC

# import packages into top-level namespace
from . import (  # noqa
    utils,
    constants,
    dynamics,
    # phot,
    # sc,
    units,
)


#############################################################################
# CONFIG FUNCTIONS


def reload_config():
    """Reload astronat configuration.

    See Also
    --------
    :mod:`~astropy.config`

    """
    config.reload_config("astronat")


# /def


#############################################################################
# HELP FUNCTIONS


def online_help(query: T.Optional[str] = None):
    """Search the online documentation for the given query.

    Opens the results in the default web browser.
    Requires an active internet connection.

    Parameters
    ----------
    query : str, optional
        The search query for `RTD <https://astronat.readthedocs.io>`_.
        None (default) or "" is an empty search.

    """
    from urllib.parse import urlencode
    import webbrowser

    # process the query
    if query is None:  # empty query, empty search
        query = ""
    else:  # encode the query
        query: str = urlencode({"q": query})

    # first get version to search
    version = __version__
    if "dev" in version:
        version = "latest"
    else:
        version = "v" + version

    url = f"https://astronat.readthedocs.io/en/{version}/search.html?{query}"

    webbrowser.open(url)

    return


# /def


@utilipy.decorators.code_dev.indev
def help(query: T.Optional[str] = None, online: bool = False):
    """astronat help function.

    Parameters
    ----------
    query : str, optional
        The search query.
    online : bool, optional
        Whether to open the online help or just print some help
        documentation (default)

    """
    if online:
        return online_help(query=query)
    # else:

    print("This function is a work in progress")

    print("".join(["-"] * 79))

    return


# /def


#############################################################################
# END

# -*- coding: utf-8 -*-

"""Attributes.

:mod:`~astropy.coordinates.attributes`

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "PotentialAttribute",
]


##############################################################################
# IMPORTS

import typing as T

from astropy.coordinates import Attribute

from ..common import PotentialType

# try:  # TODO do this in a central location
#     from galpy.potential import Potential
# except ImportError:
#     _HAS_GALPY = False
# else:
#     galpy.__version__
#     _HAS_GALPY = True

# try:  # TODO do this in a central location
#     import gala
# except ImportError:
#     _HAS_GALA = False
# else:
#     gala.__version__
#     _HAS_GALA = True


##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


class PotentialAttribute(Attribute):
    """A frame attribute that is a Potential.

    Can be `None`, which should be used for special cases in associated
    frame transformations like "this quantity should be ignored" or similar.

    Parameters
    ----------
    default : value or Quantity or None
        Default value for the attribute if the user does not supply one. If a
        Quantity, it must be consistent with ``unit``, or if a value, ``unit``
        cannot be None.
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.

    """

    def convert_input(self, value: T.Union[None, PotentialType]):
        """Convert input value to a Potential object and validate.

        .. todo::

            validate value

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.

        """
        if value is None:  # nothing to convert
            return None, False
        elif isinstance(value, PotentialAttribute):  # already good
            return None, False

        # TODO validate
        out = value
        converted = False

        return out, converted

    # /def


# /class


# -------------------------------------------------------------------


##############################################################################
# END

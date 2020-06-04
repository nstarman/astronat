# -*- coding: utf-8 -*-

"""Coordinate Transformations.

https://docs.astropy.org/en/stable/api/astropy.coordinates.CoordinateTransform.html

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["Astropy"]

__all__ = [
    "FunctionWithKwargTransform",
    "FunctionWithKwargTransformWithFiniteDifference",
]


##############################################################################
# IMPORTS

# BUILT-IN

import functools
import typing as T
from warnings import warn


# THIRD PARTY

from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates.representation import BaseRepresentation
from astropy.coordinates.transformations import (
    FunctionTransform,
    FunctionTransformWithFiniteDifference,
)
import astropy.units as u

from astropy.utils.decorators import format_doc

# PROJECT-SPECIFIC


##############################################################################
# PARAMETERS

# RepresentationType = T.TypeVar[BaseRepresentation]

##############################################################################
# CODE
##############################################################################


@format_doc(None, original_docstring=FunctionTransform.__doc__)
class FunctionWithKwargTransform(FunctionTransform):
    """FunctionTransform generalized to accept keyword arguments.

    {original_docstring}

    """

    @format_doc(
        FunctionTransform.__init__.__doc__
    )  # TODO add func_kwargs parameters
    def __init__(
        self,
        func: T.Callable,
        fromsys,
        tosys,
        priority: int = 1,
        register_graph=None,
        func_kwargs: T.Dict[str, T.Any] = {},
    ):
        """Initalize.

        Other Parameters
        ----------------
        func_kwargs : dict
            arguments into `func`

        """
        # store the func_kwargs
        self.func_kwargs = func_kwargs
        # construct a partial version of func with these kwargs
        # this can stand in for func anywhere.
        func_partial = functools.partial(func, **self.func_kwargs)

        super().__init__(
            func_partial,
            fromsys,
            tosys,
            priority=priority,
            register_graph=register_graph,
        )

    # /def


# /class


# -------------------------------------------------------------------


@format_doc(
    None, original_docstring=FunctionTransformWithFiniteDifference.__doc__
)
class FunctionWithKwargTransformWithFiniteDifference(
    FunctionTransformWithFiniteDifference
):
    r"""FunctionTransformWithFiniteDifference with keyword arguments.

    {original_docstring}

    """

    @format_doc(
        FunctionTransformWithFiniteDifference.__init__.__doc__
    )  # TODO add func_kwargs parameters
    def __init__(
        self,
        func: T.Callable,
        fromsys,
        tosys,
        priority: int = 1,
        register_graph=None,
        finite_difference_frameattr_name: str = "obstime",
        finite_difference_dt=1 * u.second,
        symmetric_finite_difference=True,
        func_kwargs: T.Dict[str, T.Any] = {},
    ):
        """Initalize.

        Other Parameters
        ----------------
        func_kwargs : dict
            arguments into `func`

        """
        # store the func_kwargs
        self.func_kwargs = func_kwargs
        # construct a partial version of func with these kwargs
        # this can stand in for func anywhere.
        func_partial = functools.partial(func, **self.func_kwargs)

        super().__init__(
            func_partial,
            fromsys,
            tosys,
            priority=priority,
            register_graph=register_graph,
            finite_difference_frameattr_name=finite_difference_frameattr_name,
            finite_difference_dt=finite_difference_dt,
            symmetric_finite_difference=symmetric_finite_difference,
        )

    # /def


# /class

##############################################################################
# END

# -*- coding: utf-8 -*-

"""Units Utilities."""

__author__ = "Nathaniel Starkman"


__all__ = [
    "quantity_return_",
]


##############################################################################
# IMPORTS

# BUILT-IN

import textwrap

from typing import Union, Sequence, Any, TypeVar


# THIRD PARTY

from astropy.units.core import *

import astropy.units as u
from astropy.units.core import Unit, IrreducibleUnit
from astropy.utils.decorators import format_doc

from utilipy.utils.typing import UnitableType


##############################################################################
# PARAMETERS

_doc_base_params = """unit: :class:`~astropy.units.Unit`, optional
    sets the unit for the returned value.
    if None, returns value unchanged, unless `to_value` is used
    if blank string, decomposes
to_value: bool, optional
    whether to return ``.to_value(unit)``
    see ``astropy.units.Quantity.to_value``
equivalencies: list, optional
    equivalencies for ``.to()`` and ``.to_value()``
    only used if `unit` to `to_value` are not None/False
decompose: bool or list, optional
    unit decomposition
    default, False

    * bool: True, False for decomposing.
    * list: bases for ``.decompose(bases=[])``.
      Will first decompose, then apply `unit`, `to_value`, `equivalencies`.

    Decomposing then converting wastes time, since
    ``.to(unit, equivalencies)`` internally does conversions.
    The only use for combining `decompose` with other ``quantity_return_``
    parameters is with

    .. code-block:: python

        unit=None, to_value=True, equivalencies=[]

    since this will decompose to desired bases then return the value in those bases

    .. note::

        **experimental feature:**
        for things which are not (:class:`~astropy.units.Unit`),
        tries wrapping in ``Unit()``. This would normally return an error, but now
        allows for conversions such as:

        >>> x = 10 * u.km * u.s
        >>> bases = [u.Unit(2 * u.km), u.s]
        >>> x.decompose(bases=bases) # doctest: +SKIP
        <Quantity 5.0 2 km s>"""

_doc_base_raises = """
:class:`~ValueError`
    if unit not astropy compatible
:class:`~astropy.units.UnitConversionError`
    if conversion not legit
"""


###############################################################################
# CODE
###############################################################################


@format_doc(
    None,
    parameters=textwrap.indent(_doc_base_params, "    "),
    raises=textwrap.indent(_doc_base_raises, "    "),
)
def quantity_return_(
    res: Any,
    unit: Union[None, UnitableType] = None,
    to_value: bool = False,
    equivalencies: list = [],
    decompose: Union[bool, Sequence] = False,
):
    """Control function return of :class:`~astropy.units.Quantity`.

    Parameters
    ----------
    res: :class:`~astropy.units.Quantity`, optional
        the result
    {parameters}

    Returns
    -------
    res:
        function output, converted / decomposed / evaluated to desired units

    Raises
    ------
    {raises}

    Examples
    --------
    How to apply in a function directly

    >>> def example_function(x, **kw):
    ...     return quantity_return_(x, unit=kw.get('unit', None),
    ...                             to_value=kw.get('to_value', False),
    ...                             equivalencies=kw.get('equivalencies', []),
    ...                             decompose=kw.get('decompose', []))
    >>> example_function(10*u.km, unit=u.m, to_value=True)
    10000.0

    """
    # fast checks to do nothing
    # nothing required
    if not hasattr(res, "to"):
        return res
    # nothing asked
    elif (
        (unit is None)
        & (to_value is False)
        & (equivalencies == [])
        & (decompose is False)
    ):
        return res

    # --------------------
    # Decomposing

    if decompose is False:
        pass
    elif decompose is True:
        res = res.decompose()
    elif decompose:  # decompose is NOT empty list
        clss = (Unit, IrreducibleUnit)
        bases = [
            Unit(x) if not issubclass(x.__class__, clss) else x
            for x in decompose
        ]
        res = res.decompose(bases=bases)

    # --------------------
    # Returning

    if (unit is None) and (to_value is False):  # nothing further required
        pass
    elif to_value is True:  # return value
        unit = u.Unit(unit) if unit is not None else unit
        res = res.to_value(unit, equivalencies=equivalencies)
    else:  # return, with unit conversion
        if equivalencies:  # has equivalencies, must use "to"
            res = res.to(u.Unit(unit), equivalencies=equivalencies)
        else:  # no equivalencies, can convert in-place
            res <<= u.Unit(unit)

    return res


# /def


# ------------------------------------------------------------------------


###############################################################################
# END

# -*- coding: utf-8 -*-

"""Table Utilities."""


__all__ = [
    "rename_columns",
    "cast_columns",
]


##############################################################################
# IMPORTS

import typing as T
import warnings
from collections.abc import Sequence

from utilipy.utils.typing import TableType

##############################################################################
# CODE
##############################################################################


def rename_columns(
    table: TableType, rename: T.Dict[str, str], warn_mismatch: bool = True
) -> TableType:
    """Rename columns in Table.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
    rename : dict
        Rename columns. Keys are current names, values are new names.
        Store old names in column meta and rename table in `table` meta.
    warn_mismatch : bool
        Warn if keys in `rename` not columns in `table`.

    Returns
    -------
    :class:`~astropy.table.Table`
        `table`, columns renamed in-place.

    """
    renamed = rename.copy()  # copy of rename dict, for putting in meta.

    # iterate through dictionary, renaming columns
    name: str
    new_name: str
    for name, new_name in rename.items():

        if name in table.colnames:
            table.rename_column(name=name, new_name=new_name)

            # store original name in Table column meta
            # fails if QTable, since columns are Quantities, which have no meta
            if hasattr(table[new_name], "meta"):
                table[new_name].meta["rename"] = name

        # if name is not in `table` (and `warn_mismatch`) raise warning.
        else:
            if warn_mismatch:
                warnings.warn(f"{name} not in table.")

            del renamed[name]  # del so don't add to meta

        # /if

    # /for

    table.meta["rename"] = renamed  # store rename dict in meta

    return table


# /def


# ------------------------------------------------------------------------


def cast_columns(
    table: TableType,
    recast: T.Dict[str, T.Union[T.Any, T.Tuple[T.Any, bool]]],
    warn_mismatch: bool = True,
    elementwise: bool = False,
) -> TableType:
    """Cast Table column types.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
    rename : dict
        Rename columns. Keys are current names, values are new names.
        Store old names in column meta and rename table in Table meta.
    warn_mismatch : bool
        Warn if keys in `rename` not columns in `table`.

    Returns
    -------
    :class:`~astropy.table.Table`
        `table`, columns re-typed in-place.

    """
    recast = recast.copy()  # copy of rename dict, for putting in meta.

    # iterate through dictionary, renaming columns
    names: T.Tuple[str] = tuple(recast.keys())
    name: str

    for name in names:  # need to preload so don't change size
        new_type = recast[name]

        if name not in table.colnames:
            if warn_mismatch:
                warnings.warn(f"{name} not in table.")
            del recast[name]

            continue
        # /def

        if isinstance(new_type, Sequence):
            new_type, _elementwise = new_type
        else:
            _elementwise = elementwise

        if not _elementwise:
            # first try applying to whole column
            # this has the advantage of not creating ``object``s as
            # the value in each row of the column
            try:
                col = table[name].astype(new_type)
            # if that fails, apply function element-wise
            except Exception:
                try:
                    col = new_type(table[name])
                except Exception:
                    _elementwise = True
                else:
                    table[name] = col
            else:
                table[name] = col

        if _elementwise:
            table[name] = [new_type(x) for x in table[name]]

    # /for

    table.meta["recast"]: str = str(recast)  # store recast in meta
    # TODO change function

    return table


# /def


# ------------------------------------------------------------------------


##############################################################################
# END

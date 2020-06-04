# -*- coding: utf-8 -*-

"""Table Utilities."""


__all__ = [
    "rename_columns",
    "cast_columns",
]


##############################################################################
# IMPORTS

# BUILT-IN

import typing as T
import warnings


# THIRD PARTY


# PROJECT-SPECIFIC


##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


def rename_columns(
    table, rename: T.Dict[str, str], warn_mismatch: bool = True
):
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

    # /for

    table.meta["rename"] = renamed  # store rename dict in meta

    return table


# /def


# ------------------------------------------------------------------------


def cast_columns(
    table,
    recast: T.Dict[str, str],
    warn_mismatch: bool = True,
    elementwise=False,
):
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
    for name, new_type in recast.items():

        if isinstance(new_type, (list, tuple)):
            new_type, _elementwise = new_type
        else:
            _elementwise = elementwise

        if name in table.colnames:
            if not _elementwise:
                # first try applying to whole column
                # this has the advantage of not creating ``object``s as
                # the value in each row of the column
                try:
                    table[name] = new_type(table[name])
                # if that fails, apply function element-wise
                except Exception:
                    table[name] = [new_type(x) for x in table[name]]
            else:
                table[name] = [new_type(x) for x in table[name]]

        else:
            if warn_mismatch:
                warnings.warn(f"{name} not in table.")

            del recast[name]

    # /for

    table.meta["recast"] = str(recast)  # store rename in meta

    return table


# /def


# ------------------------------------------------------------------------


##############################################################################
# END

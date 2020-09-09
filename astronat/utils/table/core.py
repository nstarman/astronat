# -*- coding: utf-8 -*-

"""Table Core Functions.

.. todo::

    - astropy table-like-objects
    for method of reducing a TableList into a single table.

    - change reader and writer to astropy asdf
    which can serialize a series of tables.

    - have an export to astropy NDData class
    need to understand uncertainties type. metadata?, NDDataArray

"""

__all__ = [
    "TableList",
    "QTableList",
    "TablesList",
]


##############################################################################
# IMPORTS

import collections.abc as cabc
import pathlib
import typing as T
import warnings
from collections import OrderedDict

import asdf
from astropy.table import QTable, Table
from astropy.utils.collections import HomogeneousList
from astropy.utils.introspection import resolve_name
from astropy.utils.metadata import MetaData
from astroquery.exceptions import InputWarning
from utilipy.utils.typing import OrderedDictType

##############################################################################
# CODE
##############################################################################


class TablesList(HomogeneousList):
    """Grouped Tables.

    A subclass of list that contains only elements of a given type or
    types.  If an item that is not of the specified type is added to
    the list, a `TypeError` is raised. Also includes some pretty printing
    methods for an OrderedDict of :class:`~astropy.table.Table` objects.

    """

    _types = None

    meta = MetaData(copy=False)

    def __init__(
        self,
        inp: OrderedDictType = [],
        *,
        name: T.Optional[str] = None,
        reference: T.Optional[T.Any] = None,
        **metadata,
    ):
        """Astroquery-style table list.

        Parameters
        ----------
        inp : sequence, optional
            An initial set of tables.
        name : str, optional
            name of the list of tables.
        reference : citation, optional
            citation.
        **metadata : Any
            arguments into meta

        """
        # meta
        self.meta["name"] = name
        self.meta["reference"] = reference
        for k, v in metadata.items():
            self.meta[k] = v

        inp = self._validate(inp)  # ODict, & ensure can assign values

        # Convert input to correct to type
        # If None, can be anything
        if self._types is not None:
            for k, val in inp.items():
                inp[k] = self._types(val)  # TODO handle multiple "_types"

        # finally add the input
        # _dict store the indices for the keys
        self._dict = {k: i for i, k in enumerate(inp.keys())}

        # need to bypass HomogeneousList init, which uses ``extend``
        list.__init__(self, inp.values())

    # /def

    # -----------------

    def _assert(self, x):
        """Check `x` is correct type (set by _type)."""
        if self._types is None:  # allow any type
            return
        super()._assert(x)

    # /def

    def _validate(self, value):
        """Validate `value` compatible with table."""
        if isinstance(value, TablesList):  # tablelist or subclass
            pass
        elif not isinstance(value, OrderedDict):
            try:
                value = OrderedDict(value)
            except (TypeError, ValueError):
                raise ValueError(
                    "Input to TableList must be an OrderedDict "
                    "or list of (k,v) pairs"
                )

        return value

    # /def

    # -----------------
    # Properties

    @property
    def name(self) -> str:
        """Name."""
        return self.meta["name"]

    # /def

    @property
    def __reference__(self):
        """Get reference from metadata, if exists."""
        return self.meta.get("reference", None)

    # /def

    # -----------------
    # Dictionary methods

    def keys(self):
        """Set-like object giving table names."""
        return self._dict.keys()

    # /def

    def sortedkeys(self):
        """Set-like object giving table names.

        Ordered by value. Does not update with table.

        """
        sorted_dict = dict(
            sorted(self._dict.items(), key=lambda item: item[0])
        )
        return cabc.KeysView(sorted_dict.keys())

    # /def

    def values(self):
        """Tuple object providing a view on tables.

        Note that the tables can be edited

        """
        return tuple(self)

    # /def

    def items(self):
        """Generator providing iterator over name, table."""
        return cabc.ItemsView(zip(self.keys(), self.values()))

    # /def

    # -----------------
    # Get / Set

    def index(self, key: str) -> int:
        """Index of `key`.

        Parameters
        ----------
        key : str

        Returns
        -------
        int

        """
        return self._dict[key]

    # /def

    def __getitem__(self, key: T.Union[int, slice, str]):
        """Get item or slice.

        Parameters
        ----------
        key : str or int or slice
            if str, the dictionary key.
            if int, the dictionary index
            if slice, slices dictionary.values
            supports string as slice start or stop

        Returns
        -------
        Table

        Raises
        ------
        TypeError
            if key is not int or key

        """
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, slice):
            start, stop = key.start, key.stop
            # string replacement for start, stop values
            # replace by int
            if isinstance(start, str):
                start = self.index(start)
            if isinstance(stop, str):
                stop = self.index(stop)
            key = slice(start, stop, key.step)
            return super().__getitem__(key)
        else:
            return super().__getitem__(self.index(key))

    # /def

    def __setitem__(self, key: str, value):
        """Set item, but only if right type (managed by super)."""
        if not isinstance(key, str):
            raise TypeError

        # first try if exists
        if key in self._dict:
            ind = self._dict[key]
            return super().__setitem__(ind, value)  # (super _assert)

        # else append to end
        else:
            ind = len(self)
            self._dict[key] = ind
            return super().append(value)  # (super _assert)

    # /def

    def __delitem__(self, key):
        """Delete Item. Forbidden."""
        raise NotImplementedError("Forbidden.")

    # /def

    def update(self, other):
        """Update TableList using OrderedDict update method."""
        values = self._validate(other)  # first make sure adding a key-val pair
        for k, v in values.items():  # TODO better
            self[k] = v  # setitem manages _dict

    # /def

    def extend(self, other):
        """Extend TableList. Unlike update, cannot have duplicate keys."""
        values = self._validate(other)  # first make sure adding a key-val pair

        if any((k in self.keys() for k in values.keys())):
            raise ValueError("cannot have duplicate keys")

        self.update(values)

    # /def

    def __iadd__(self, other):
        """Add in-place."""
        return super().__iadd__(other)

    # /def

    def append(self, key: str, value):
        """Append, if unique key and right type (managed by super)."""
        if key in self._dict:
            raise ValueError("cannot append duplicate key.")

        self._dict[key] = len(self)
        return super().append(value)

    # /def

    def pop(self):
        """Pop. Forbidden."""
        raise NotImplementedError("Forbidden.")

    # /def

    def insert(self, value):
        """Insert. Forbidden."""
        raise NotImplementedError("Forbidden.")

    # /def

    # -----------------
    # string representation

    def __repr__(self):
        """String representation.

        Overrides the `OrderedDict.__repr__` method to return a simple summary
        of the `TableList` object.

        Returns
        -------
        str

        """
        return self.format_table_list()

    # /def

    def format_table_list(self) -> str:
        """String Representation of list of Tables.

        Prints the names of all :class:`~astropy.table.Table` objects, with
        their respective number of row and columns, contained in the
        `TableList` instance.

        Returns
        -------
        str

        """
        ntables = len(list(self.keys()))
        if ntables == 0:
            return "Empty {cls}".format(cls=self.__class__.__name__)

        header_str = "{cls} with {keylen} tables:".format(
            cls=self.__class__.__name__, keylen=ntables
        )
        body_str = "\n".join(
            [
                "\t'{t_number}:{t_name}' with {ncol} column(s) "
                "and {nrow} row(s) ".format(
                    t_number=t_number,
                    t_name=t_name,
                    nrow=len(self[t_number]),
                    ncol=len(self[t_number].colnames),
                )
                for t_number, t_name in enumerate(self.keys())
            ]
        )

        return "\n".join([header_str, body_str])

    # /def

    def print_table_list(self):
        """Print Table List.

        calls ``format_table_list``

        """
        print(self.format_table_list())

    # /def

    def pprint(self, **kwargs):
        """Helper function to make API more similar to astropy.Tables.

        .. todo::

            uses "kwargs"

        """
        if kwargs != {}:
            warnings.warn(
                "TableList is a container of astropy.Tables.", InputWarning
            )

        self.print_table_list()

    # /def

    # -----------------
    # I/O

    def _save_table_iter(self, format, **table_kw):
        for i, name in enumerate(self.keys()):  # order-preserving

            # get kwargs for table writer
            # first get all general keys (by filtering)
            # then update with table-specific dictionary (if present)
            kw = {
                k: v for k, v in table_kw.items() if not k.startswith("table_")
            }
            kw.update(table_kw.get("table_" + name, {}))

            if isinstance(format, str):
                fmt = format
            else:
                fmt = format[i]

            yield name, fmt, kw

    # /def

    def write(
        self,
        drct: str,
        format="asdf",
        split=True,
        serialize_method=None,
        **table_kw,
    ):
        """Write to ASDF.

        Parameters
        ----------
        drct : str
            The main directory path.
        format : str or list, optional
            save format. default "asdf"
            can be list of same length as TableList
        split : bool, optional
            *Applies to asdf `format` only*

            Whether to save the tables as individual file
            with `file` coordinating by reference.

        serialize_method : str, dict, optional
            Serialization method specifier for columns.

        **table_kw
            kwargs into each table.
            1. dictionary with table name as key
            2. General keys

        """

        # -----------
        # Path checks

        path = pathlib.Path(drct)

        if path.suffix == "":  # no suffix
            path = path.with_suffix(".asdf")

        if path.suffix != ".asdf":  # ensure only asdf
            raise ValueError("file type must be `.asdf`.")

        drct = path.parent  # directory in which to save

        # -----------
        # TableType

        if self._types is None:
            table_type = [
                tp.__class__.__module__ + "." + tp.__class__.__name__
                for tp in self.values()
            ]
        else:
            table_type = self._types.__module__ + "." + self._types.__name__

        # -----------
        # Saving

        TL = asdf.AsdfFile()
        TL.tree["meta"] = tuple(self.meta.items())
        TL.tree["table_names"] = tuple(self.keys())  # in order
        TL.tree["save_format"] = format
        TL.tree["table_type"] = table_type

        if format == "asdf" and not split:  # save as single file
            for name in self.keys():  # add to tree
                TL.tree[name] = self[name]

        else:  # save as individual files
            for name, fmt, kw in self._save_table_iter(format, **table_kw):

                # name of table
                table_path = drct.joinpath(name)
                if table_path.suffix == "":  # TODO currently always. CLEANUP
                    table_path = table_path.with_suffix(
                        "." + fmt.split(".")[-1]
                    )

                # internal save
                if format == "asdf":
                    kw["data_key"] = kw.get("data_key", name)
                self[name].write(
                    table_path,
                    format=fmt,
                    serialize_method=serialize_method,
                    **kw,
                )

                # save by relative reference
                if format == "asdf":
                    with asdf.open(table_path) as f:
                        TL.tree[name] = f.make_reference(path=[name])
                else:
                    TL.tree[name] = str(table_path.relative_to(drct))

        # Need to add a "data" key to not break asdf
        if not format == "asdf":
            TL.tree["data"] = TL.tree["table_names"]

        # /if
        TL.write_to(str(path))  # save directory

    # /def

    @classmethod
    def _read_table_iter(cls, f, format, **table_kw):
        names = f.tree["table_names"]
        # table type, for casting
        # so that QTableList can open a saved TableList correctly
        # TablesList specifies no type, so must rely on saved info
        if cls._types is None:
            table_type = f.tree["table_type"]
            if not isinstance(table_type, cabc.Sequence):
                table_type = [table_type] * len(names)
            ttypes = [resolve_name(t) for t in table_type]
        else:
            ttypes = [cls._types] * len(names)

        for i, name in enumerate(names):  # order-preserving

            if isinstance(format, str):
                fmt = format
            else:
                fmt = format[i]

            # get kwargs for table writer
            # first get all general keys (by filtering)
            # then update with table-specific dictionary (if present)
            kw = {
                k: v for k, v in table_kw.items() if not k.startswith("table_")
            }
            kw.update(table_kw.get("table_" + name, {}))

            yield name, ttypes[i], fmt, kw

    # /def

    @classmethod
    def read(
        cls,
        drct: str,
        format: T.Union[str, T.Sequence] = None,
        suffix: T.Optional[str] = None,
        **table_kw,
    ):
        """Write to ASDF.

        Parameters
        ----------
        drct : str
            The main directory path.
        format : str or list, optional
            read format. default "asdf"
            can be list of same length as TableList
        suffix : str, optional
            suffix to apply to table names.
            will be superceded by an "fnames" argument, when added

        **table_kw
            kwargs into each table.
            1. dictionary with table name as key
            2. General keys

        """
        # -----------
        # Path checks

        path = pathlib.Path(drct)

        if path.suffix == "":  # no suffix
            path = path.with_suffix(".asdf")

        if path.suffix != ".asdf":  # ensure only asdf
            raise ValueError("file type must be `.asdf`.")

        drct = path.parent  # directory

        # -----------
        # Read

        TL = cls()
        with asdf.open(path) as f:
            f.resolve_references()

            # load in the metadata
            TL.meta = OrderedDict(f.tree["meta"])

            if format is None:
                format = f.tree["save_format"]

            # iterate through tables
            for name, ttype, fmt, kw in cls._read_table_iter(
                f, format, **table_kw
            ):
                tl = f.tree[name]

                # TODO what if tuple of str as path to name?
                if not isinstance(tl, str):  # only for asdf
                    TL[name] = ttype(f.tree[name], **kw)  # TODO need kw?
                else:
                    table_path = drct.joinpath(tl)  # .with_suffix(suffix)
                    TL[name] = ttype.read(str(table_path), format=fmt, **kw)

        return TL

    # /def

    def copy(self):
        """Shallow copy."""
        out = self.__class__(self)
        out.meta = self.meta

        return out

    # /def


# /class

# -------------------------------------------------------------------


class TableList(TablesList):
    """Homogeneous TablesList."""

    _types = Table


# /class


# -------------------------------------------------------------------


class QTableList(TablesList):
    """Astroquery-style QTable list.

    Attributes
    ----------
    meta : :class:`~astropy.utils.metadata.MetaData`

    """

    _types = QTable


# /class


##############################################################################
# END

# -*- coding: utf-8 -*-

"""Table Core Functions.

.. todo::

    - astropy table-like-objects
    for method of reducing a TableList into a single table.

    - change reader and writer to astropy asdf
    which can serialize a series of tables.

    - Inherit from ReferenceClass to get automatic meta and __reference__

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

# BUILT-IN

import json
import pathlib
import typing as T
import warnings

from collections import OrderedDict


# THIRD PARTY

from astropy.table import Table, QTable
from astropy.utils.introspection import resolve_name
from astropy.io.misc import fnpickle, fnunpickle

from astroquery.exceptions import InputWarning

from utilipy.utils.collections import ReferenceBase


##############################################################################
# PARAMETERS


def _type_pass_thru(x):
    return x


# TODO only necessary because `typing` doesn't have OrderedDict for py3.6
OrderedDictType = T.TypeVar(
    "OrderedDictType", OrderedDict, T.Sequence[T.Tuple[str, T.Any]]
)


##############################################################################
# CODE
##############################################################################


class TableList(OrderedDict, ReferenceBase):
    """docstring for TableList."""

    _types = Table
    # meta = MetaData()  # in ReferenceBase

    def __init__(
        self,
        inp: OrderedDictType = [],
        *,
        name: T.Optional[str] = None,
        reference=None,
        **kwargs,
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
        **kwargs : Any
            arguments into meta

        """
        # meta
        self.meta["name"] = name
        self.meta["reference"] = reference
        for k, v in kwargs.items():
            self.meta[k] = v

        super().__init__()
        inp = self._validate(inp)  # ODict, & ensure can assign values

        # convert input to correct to type
        if self._types is None:
            converter = _type_pass_thru
        elif not callable(self._types):
            converter = self._types[0]
        else:
            converter = self._types
        for k, val in inp.items():
            inp[k] = converter(val)

        # finally add the input.
        self.update(inp)

        return

    # /def

    # -----------------

    def _assert(self, x, name=None):
        """Check `x` is correct type (set by _types)."""
        if self._types is None:  # allow any type
            pass
        elif not isinstance(x, self._types):  # catch wrong types
            raise TypeError(
                "Object {} not of type '{}'".format(name or "", self._types)
            )

    # /def

    def _validate(self, value):
        """Validate `value` compatible with table and assert correct type."""
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

        # for key, val in value.items():
        #     print(type(val))
        #     self._assert(val, name=key)

        return value

    # /def

    # -----------------
    # Properties

    @property
    def name(self):
        """Name."""
        return self.meta["name"]

    # /def

    # -----------------
    # Get / Set

    def index(self, key: str):
        """Index of `key`.

        Parameters
        ----------
        key : str

        Returns
        -------
        int

        """
        return list(self.keys()).index(key)

    # /def

    def __getslice__(self, slicer):
        """Slice dictionary.

        slice(star, stop, step)
        supports string as slice start or stop

        Examples
        --------
        First creating slices
            >>> vizier = Table([[1], [2]], names=["a", "b"])
            >>> simbad = Table([[3], [5]], names=["d", "e"])
            >>> od = TableList([("vizier", vizier), ("simbad", simbad)],
            ...                name="test")

        Example of standard slicing

            >>> od[1:]  # doctest: +SKIP
            Table([[3], [5]], names=["d", "e"])

        Now slicing with string

            >>> od["vizier":]  # doctest: +SKIP
            Table([[3], [5]], names=["d", "e"])


        """
        # check if slice, if it is, then check for string inputs
        # it it isn't a slice, then still try to apply to list of values
        if isinstance(slicer, slice):
            keys = list(self.keys())  # list of keys

            start, stop = slicer.start, slicer.stop
            # string replacement for start, stop values
            # replace by int
            if isinstance(start, str):
                start = keys.index(start)
            if isinstance(stop, str):
                stop = keys.index(stop)
            slicer = slice(start, stop, slicer.step)

        return list(self.values())[slicer]

    # /def

    def __getitem__(self, key):
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
            return super().__getitem__(tuple(self.keys())[key])
        elif isinstance(key, slice):
            return self.__getslice__(key)
        else:
            return super().__getitem__(key)

    # /def

    def __setitem__(self, key, value):
        """Set item, but only if right type."""
        # self._assert(value)
        super().__setitem__(key, value)

    # /def

    def update(self, other):
        """Update TableList using OrderedDict update method."""
        values = self._validate(other)  # first make sure adding a key-val pair
        # print(values)
        OrderedDict.update(self, values)  # TODO use super

    # /def

    def extend(self, other):
        """Extend TableList."""
        self.update(other)

    # /def

    def __iadd__(self, other):
        """Add in-place."""
        self.extend(other)
        return self

    # /def

    def append(self, other):
        """Forbidden to append."""
        raise Exception("Forbidden to append.")

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

    def format_table_list(self):
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
            return "Empty TableList"

        header_str = "TableList with {keylen} tables:".format(keylen=ntables)
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

        calls format_table_list

        """
        print(self.format_table_list())

    # /def

    def pprint(self, **kwargs):
        """Helper function to make API more similar to astropy.Tables."""
        if kwargs != {}:
            warnings.warn(
                "TableList is a container of astropy.Tables.", InputWarning
            )
        self.print_table_list()

    # /def

    # -----------------
    # I/O

    def _renamedkeycopy(self, name=None):
        """Shallow copy with keys renamed."""
        if name is None:
            name = self.name
        out = self.__class__(((name + "_" + k, v) for k, v in self.items()))
        out.meta = self.meta
        return out

    # /def

    def copy(self):
        """Shallow copy."""
        out = self.__class__(self)
        out.meta = self.meta
        return out

    # /def

    def write(self, drct: str, *args, **kwargs):
        """Write TableList.

        .. todo::
            write table meta better than pickle
            correct bad characters is figured out names

        Parameters
        ----------
        drct : str
            the output directory.
        *args : tuple, optional
            Positional arguments passed through to data writer.
            If supplied the first argument is the output filename.

                - None will make own file names, trying `tables` in meta
                - List must match TableList length and order
        format : str
            File format specifier.
        serialize_method : str, dict, optional
            Serialization method specifier for columns.

        """
        # arguments are into each table. normally they accept the filename
        # as the first argument. now the first argument can be a list of
        # filenames. The rest of the arguments are passed to each table.
        # if no arguments are given, then "fnames" in metadata is used to
        # try and figure out the names. If that doesn't exist it will name
        # the tables by index.
        if not args:  # args are not None, need to get name list.
            if "fnames" in self.meta.keys():
                fnames = self.meta["fnames"]  # name list
                subargs = list(args[1:])  # rest of args

            else:  # need to figure out names
                raise ValueError("TODO, figure out names")
                # fnames = self.keys()  # TODO correct bad characters
                # subargs = ()
        else:
            fnames = self.meta["fnames"]
            subargs = ()
        # /if

        # Step 0: information
        tli = {
            "name": self.name,
            "format": kwargs.get("format", None),
            "tables": [(k, v) for k, v in zip(self.keys(), fnames)],
            "meta": "meta.pkl",
        }
        # table_type
        if not callable(self._types):  # it's a list, choose the first
            table_type = self._types[0]
        else:
            table_type = self._types
        table_type = (table_type.__module__ + "." + table_type.__name__,)
        tli["table_type"] = table_type

        # Step 1: directory
        # make or determine if should overwrite
        drct = pathlib.Path(drct)
        try:
            drct.mkdir(parents=True)

        except FileExistsError:  # directory exists
            # now need to determine whether to allow overwriting
            if not kwargs.get("overwrite", False):
                raise FileExistsError(f"director {drct} already exists.")

        # Step 2: Save JSON TableList descriptor
        with open(
            drct.joinpath("tablelist.json"), "w", encoding="utf-8"
        ) as file:
            json.dump(tli, file, ensure_ascii=False, indent=4)

        # Step 3: write meta
        with open(drct.joinpath(tli["meta"]), "wb") as file:
            fnpickle(self.meta, file)

        # Step 4: write individual tables
        # print(tli["tables"])
        for key, fname in tli["tables"]:
            self[key].write(drct.joinpath(fname), *subargs, **kwargs)
        # /for

        return

    # /def

    @staticmethod
    def validate_tablelist_json(tli: T.Union[str, dict]):
        """Checks if JSON is a valid tablelistinfo."""
        if isinstance(tli, str):  # open if file path
            with open(tli, "r") as file:
                tli = json.load(file)

        keys = tli.keys()

        # first check has valid keys
        if not all(
            [k in keys for k in ("name", "format", "tables", "table_type")]
        ):
            raise ValueError("does not have write keys")
        # now check values are legit
        if not isinstance(tli["name"], str):
            raise TypeError("name is not a str")
        if not isinstance(tli["format"], (str, type(None))):
            raise TypeError("format is not str or None")
        if not isinstance(tli["tables"], (list, tuple)):
            raise TypeError("tables must be list or tuple")

        for t in tli["tables"]:
            if not isinstance(t[0], str) or not isinstance(t[1], str):
                raise TypeError("table element must be strings")

        if not isinstance(tli["table_type"], (str, list, tuple)):
            raise TypeError("table_type must be str, list, or tuple")

    # /def

    @classmethod
    def read(cls, tli, *args, **kwargs):
        """Read and parse a data directory and return as a TableList.

        .. todo::

            support providing in all the info contained in tablist json

        Parameters
        ----------
        tli : str or dict
            the tablelist info
                - name
                - format: table format
                - tables: (key, save location) tuples
                - meta: metadata pickle location
        format : str
            File format specifier.
        *args : tuple, optional
            Positional arguments passed through to data reader
            Do NOT include input filename.
        **kwargs : dict, optional
            Keyword arguments passed through to data reader.

        Returns
        -------
        TableList
            with name, meta, and contents set by `tli` info

        """
        # --------------
        # loader (it's a json file)

        if isinstance(tli, str):  # open if file path
            # get directory from file name
            drct = pathlib.Path(tli).resolve()
            if not drct.is_dir():
                drct = drct.parent
            # json loader
            with open(tli, "r") as file:
                tli = json.load(file)
        else:
            drct = pathlib.Path(".")  # paths relative

        cls.validate_tablelist_json(tli)

        # --------
        # construct tablelist

        tablelist = OrderedDict()

        # read tables
        # pop format, with tli value as default
        fmt = kwargs.pop("format", tli.get("format", None))

        if cls._types is not None:
            table_type = [cls._types] * len(tli["tables"])
        else:
            tt = kwargs.pop("table_type", tli.get("table_type", "Table"))
            if isinstance(tt, str):
                _types = resolve_name(tt)
                table_type = [_types] * len(tli["tables"])
            else:
                if not len(tt) == len(tli["tables"]):
                    raise ValueError

        for reader, (key, fname) in zip(table_type, tli["tables"]):
            tablelist[key] = reader.read(
                drct.joinpath(fname), *args, format=fmt, **kwargs
            )

        # read metadata
        with open(drct.joinpath(tli["meta"]), "rb") as file:
            meta = fnunpickle(file)

        # make TableList
        TL = cls(tablelist)
        TL.meta = meta

        return TL

    # /def


# /class


# -------------------------------------------------------------------


class QTableList(TableList):
    """Astroquery-style QTable list.

    Attributes
    ----------
    meta : :class:`~astropy.utils.metadata.MetaData`
    _types : ClassType
        the element type
    _dict : dictionary
        where the keys are stored.

    """

    _types = QTable


# /class


# -------------------------------------------------------------------


class TablesList(TableList):
    """Non-homogeneous TableList."""

    _types = None

    # def _assert(self, x):  # allow any type to be added
    #     return True

    # /def


# /class


##############################################################################
# END

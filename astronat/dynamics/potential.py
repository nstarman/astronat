# -*- coding: utf-8 -*-

"""Potentials."""

__author__ = "Nathaniel Starkman"

__all__ = [
    "GalpyCompositePotential",
    # Types
    "GalpyBasePotentialType",
    "GalpyPotentialType",
]


##############################################################################
# IMPORTS

import collections.abc as abc
import typing as T
import warnings
from collections import OrderedDict

from galpy.potential import Potential

##############################################################################
# Parameters

GalpyBasePotentialType = T.TypeVar(
    "GalpyPotentialType", Potential, T.Sequence[Potential],
)

##############################################################################
# CODE
##############################################################################


class GalpyCompositePotential(list):
    """A composite Galpy potential composed of distinct components.

    Galpy supports composition through lists, and this class is fully
    compatible with all list methods, as it subclasses list. Additionally,
    many dictionary methods are supported, such as item access by key, view
    methods like ``.keys()`` and ``.values()``, and some other convenience
    methods.

    Attributes
    ----------
    _components : OrderedDict
        keys are component names, values are index into self

    Methods
    -------
    pack_pots
        For packing a potential, list of potentials (in a few
        different formats), or dictionary into an OrderedDict for Potential
        composition.
    __getitem__
        In addition to all standard list inputs, allows component names (str,
        in `_components`) as keys or start, stop value in a slice object.
    keys
        same as in dict
    values
        same as in dict
    valuesdict
        OrderedDict of component name, Potential

    """

    def __init__(
        self,
        pots: T.Optional[GalpyBasePotentialType] = None,  # TODO also self type
        names: T.Union[str, T.Sequence, None] = None,
        **components,
    ):
        """Create Composite Galpy Potential.

        Parameters
        ----------
        pots : Potential or Mapping or Sequence, optional
            Accepts a few different input options, which determine `names`'s
            input. Note this is compatible with galpy's list of potentials
            syntax.

        names : str or Sequence, optional
            component names for `pots`. Accepted argument depends on `pots`.

            - single Potential, then `names` must be a string.
            - list of potentials, `names` must be equal-length sequence or None
              If `names` is None, then will autogenerate names
            - list of components and potentials ie [("name", Potential), ...]
              means that `names` is ignored
            - dictionary likewise ignores `names`

        **components : Potential or Sequence
            Potential (or composite potential), with keys determining the
            component to which the potential corresponds.

        """
        self._components = OrderedDict()  # where the indices are stored

        if pots is not None:
            comps = self.pack_pots(pots, names=names)
        else:  # pots is None
            comps = OrderedDict()

        intersect = set(comps).intersection(components)
        if any(intersect):
            raise KeyError(
                "Component names must be unique. "
                f"Passed duplicate {intersect.keys()}"
            )
        else:
            comps.update(**components)
        # to prevent double storing the potentials,
        # have to make an indexdict from comps
        ID: T.OrderedDict[str, int]
        ID = OrderedDict([(k, i) for i, k in enumerate(comps.keys())])

        self._components = ID
        super().__init__(comps.values())

    # /def

    def __getitem__(self, key):
        """Get component by key.

        In addition to all standard list inputs, allows component names (str,
        in `_components`) as keys or start, stop value in a slice object.

        """
        if isinstance(key, str):
            key = self._components[key]

        elif isinstance(key, slice):
            if isinstance(key.start, str):
                start = self._components[key.start]
            else:
                start = key.start

            if isinstance(key.stop, str):
                stop = self._components[key.stop]
            else:
                stop = key.stop

            key = slice(start, stop, key.step)

        return super().__getitem__(key)

    # /def

    # def append(self, pots, names=None):
    #     comps = self.pack_pots(pots, names)

    #     intersect = set(comps).intersection(self._components.keys())
    #     if any(intersect):
    #         raise KeyError(
    #             "Component names must be unique. "
    #             f"Passed duplicate {intersect.keys()}"
    #         )
    #     ID = OrderedDict(
    #         [(k, len(self) + i) for i, k in enumerate(comps.keys())]
    #     )

    #     super().append(comps.values())
    #     self._components.update(ID)

    # # /def

    def extend(
        self,
        pots: T.Union[GalpyBasePotentialType, T.Sequence, T.Mapping],
        names: T.Union[str, T.Sequence, None] = None,
    ):
        """Extend Composite Potential.

        Parameters
        ----------
        pots : dict

        """
        comps = self.pack_pots(pots, names)

        intersect = set(comps).intersection(self._components.keys())
        if any(intersect):
            raise KeyError(
                "Component names must be unique. "
                f"Passed duplicate {intersect.keys()}"
            )
        ID = OrderedDict(
            [(k, len(self) + i) for i, k in enumerate(comps.keys())]
        )

        super().extend(comps.values())
        self._components.update(ID)

    # /def

    def insert(self, index: int, values):
        # need to also update the _components
        raise Exception("NOT YET IMPLEMENTED")

    # /def

    def __delitem__(self):
        # need to also update the _components
        raise Exception("NOT YET IMPLEMENTED")

    # /

    def remove(self):
        # need to also update the _components
        raise Exception("NOT YET IMPLEMENTED")

    def pop(self):
        # need to also update the _components
        raise Exception("NOT YET IMPLEMENTED")

    def sort(self):
        # need to also update the _components
        raise Exception("NOT YET IMPLEMENTED")

    def clear(self):
        super().clear()
        self._components.clear()

    # /def

    def __contains__(self, key):
        if isinstance(key, str):
            return key in self._components.keys()

        return super().__contains__(key)

    # /def

    def index(self, key):
        if isinstance(key, str):
            return self._components[key]

        return super().index(key)

    # /def

    def pack_pots(
        self,
        pots: T.Union[
            GalpyBasePotentialType, T.Any
        ],  # TODO change in GalpyCompositePotential
        names: T.Union[str, T.Sequence, None] = None,
    ):
        """Pack Galpy Potentials into an OrderedDict.

        Parameters
        ----------
        pots : Potential or Mapping or Sequence, optional
        Accepts a few different input options, which determine `names`'s input.
        Note this is compatible with galpy's list of potentials syntax.

        names : str or Sequence, optional
            component names for `pots`. Accepted argument depends on `pots`.

            - single Potential, then `names` must be a string.
            - list of potentials, `names` must be equal-length sequence or None
              If `names` is None, then will autogenerate names
            - list of components and potentials ie [("name", Potential), ...]
              means that `names` is ignored
            - dictionary likewise ignores `names`

        Returns
        -------
        OrderedDict

        """
        comps = OrderedDict()

        if isinstance(pots, Potential):  # actually just one potential
            if names is None:
                names = str(hash(pots) + len(self))
            elif not isinstance(names, str):
                raise ValueError("Potential name must be a str")

            if names in self._components.keys():
                raise KeyError(
                    "Component names must be unique. "
                    f"Passed duplicate {names}"
                )
            comps[names] = pots

            return comps

        elif isinstance(pots, abc.Mapping):
            if names is not None:
                warnings.warn("Ignoring names in favor of `potentials`' keys")
            comps.update(OrderedDict(pots))

            return comps

        # list of potentials of list of ("name", Potential)
        elif isinstance(pots, abc.Sequence):
            # determine which one by testing first element
            # TODO better validation
            if isinstance(pots[0], abc.Sequence):  # [("name", Potential), ...]
                if not isinstance(pots[0][0], str):
                    raise Exception(
                        "Nested lists of potentials. "
                        "Pass as [('name', Potential), ]"
                    )
                if names is not None:
                    warnings.warn("Ignoring `names`")

                toadd = OrderedDict(pots)

            elif not isinstance(names, T.Sequence) and names is not None:
                raise ValueError(
                    "Must pass a list of names (or None) for potentials"
                )
            elif names is None:
                toadd = OrderedDict(
                    [(str(hash(pot) + i), pot) for i, pot in enumerate(pots)]
                )

            elif len(pots) == len(names):
                toadd = OrderedDict(zip(names, pots))
            else:
                raise ValueError("`pots` & `names` list must be same length")

            comps.update(toadd)

        else:
            raise ValueError

        return comps

    # /def

    def keys(self):
        """Components names in a view."""
        return self._components.keys()

    # /def

    def values(self):
        """Components in a view."""
        return abc.ValuesView(self)

    # /def

    def valuesdict(self):
        """Ordered dict of components."""
        return OrderedDict([(k, self[i]) for k, i in self._components.items()])

    # /def

    def __repr__(self):
        # TODO detect if short enough to collapse
        # TODO properly implement nesting
        s = str(self.valuesdict())
        s = s.replace("OrderedDict([", "GalpyCompositePotential: ([\n\t")
        s = s.replace(")])", ")\n])")
        s = s.replace("),", "),\n\t")
        return s

    # /def


# /class


# -------------------------------------------------------------------


GalpyPotentialType = T.TypeVar(
    "GalpyPotentialType", GalpyBasePotentialType, GalpyCompositePotential,
)


##############################################################################
# END

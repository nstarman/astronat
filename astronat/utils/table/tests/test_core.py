# -*- coding: utf-8 -*-

"""Tests for `~astronat.utils.table.core`."""

# __all__ = [
#     # functions
#     "",
#     # other
#     "",
# ]


##############################################################################
# IMPORTS

from collections import OrderedDict

import astropy.units as u
import numpy as np
import pytest
from astropy.table import QTable, Table

from astronat.utils.table import TablesList, TableList, QTableList

import tempfile

##############################################################################
# PARAMETERS

# Quantity Tables
qt1 = QTable(data={"a": [1, 2, 3] * u.m, "b": [4, 5, 6] * u.s})
"""QTable with units."""
qt2 = Table(data={"c": [7, 8, 9], "d": [10, 11, 12]})
"""QTable with no units."""

##############################################################################
# TESTS
##############################################################################


class Test_TablesList(object):
    """Tests for `~asatronat.utils.table.core.TablesList`."""

    _cls = TablesList

    @classmethod
    def setup_class(cls):
        """Setup any state specific to the execution."""
        # tables
        cls.qt1 = qt1
        cls.qt2 = qt2

        # tables
        cls.inp = {"table1": cls.qt1, "table2": cls.qt2}

        # ordered tables
        cls.ordered_inp = OrderedDict(cls.inp)

        # the list of tables
        cls.QT = cls._cls(
            inp=cls.ordered_inp,
            name="Name",
            reference="Reference",
            extra="Extra",
        )

        cls.tempdir = tempfile.TemporaryDirectory(dir="./")

    # /def

    @classmethod
    def teardown_class(cls):
        """Teardown any state specific to the execution."""
        cls.tempdir.cleanup()

    # /def

    def test__assert_and__validate(self):
        """Test the hidden methods of TablesList."""
        # Assert
        self.QT._assert(self.qt1)
        self.QT._assert(self.qt2)

        # Validate
        # a TablesList is returned as-is
        assert self.QT._validate(self.QT) == self.QT
        # so is an ordered dictionary
        assert self.QT._validate(self.ordered_inp) == self.ordered_inp
        # a dictionary is ordered
        assert isinstance(self.QT._validate(self.inp), OrderedDict)
        # so is another compatible mapping
        assert (
            self.QT._validate(list(self.ordered_inp.items()))
            == self.ordered_inp
        )
        # things that cannot be turned into an ordered dict error
        with pytest.raises(ValueError):
            self.QT._validate(ValueError)

    # /def

    def test_attributes(self):
        """Test attributes of TablesList."""
        assert self.QT.name == "Name"
        assert self.QT.__reference__ == "Reference"

    # /def

    def test_dictionary_methods(self):
        """Test dictionary methods of TablesList."""
        assert set(self.QT.keys()) == set(self.inp.keys())

        # check equality
        for v1, v2 in zip(self.QT.values(), self.inp.values()):
            assert np.all(v1 == v2)

        for (k1, v1), (k2, v2) in zip(self.QT.items(), self.inp.items()):
            assert k1 == k2
            assert np.all(v1 == v2)

    # /def

    def test_index(self):
        """Test index methods of TablesList."""
        assert self.QT.index("table1") == 0
        assert self.QT.index("table2") == 1

    # /def

    def test_getsetdel_item(self):
        """Test getitem methods of TablesList."""
        # Get
        assert np.all(self.QT["table1"] == self.qt1)
        assert np.all(self.QT[1] == self.qt2)

        for q1, q2 in zip(self.QT[:], [self.qt1, self.qt2]):
            assert np.all(q1 == q2)

        # Set
        with pytest.raises(TypeError):
            self.QT[2] = self.qt1

        # add something that already exists
        self.QT["table1"] = self.qt1

        # add something new
        self.QT["table3"] = self.qt1

        # Del
        del self.QT["table3"]

    # /def

    def test_update(self):
        """Test update method of TablesList."""
        self.QT.update({"table3": self.qt1})

        # cleanup
        del self.QT["table3"]

    # /def

    def test_extend(self):
        """Test update method of TablesList."""
        self.QT.extend({"table3": self.qt1})

        with pytest.raises(ValueError):
            self.QT.extend({"table3": self.qt1})

        # cleanup
        del self.QT["table3"]

    # /def

    def test_iadd(self):
        """Test in-place add method of TablesList."""
        self.QT += {"table3": self.qt1}

        # cleanup
        del self.QT["table3"]

    # /def

    def test_append(self):
        """Test append method of TablesList."""
        self.QT.append("table3", self.qt1)

        with pytest.raises(ValueError):
            self.QT.append("table3", self.qt1)

        # cleanup
        del self.QT["table3"]

    # /def

    def test_pop(self):
        """Test pop method of TablesList."""
        with pytest.raises(NotImplementedError):
            self.QT.pop()

    # /def

    def test_insert(self):
        """Test insert method of TablesList."""
        with pytest.raises(NotImplementedError):
            self.QT.insert(None)

    # /def

    # ----------------------

    def test_IO(self):
        """Test read/write.

        .. todo::

            A lot more tests

        """
        self.QT.write(drct=self.tempdir.name + "/saved")

        QT = self._cls.read(self.tempdir.name + "/saved")

        for q1, q2 in zip(self.QT, QT):
            assert np.all(q1 == q2)

    # /def


# /class


#####################################################################


class Test_TableList(Test_TablesList):
    """Tests for `~asatronat.utils.table.core.TableList`."""

    _cls = TableList


# /class


#####################################################################


class Test_QTableList(Test_TablesList):
    """Tests for `~asatronat.utils.table.core.QTableList`."""

    _cls = QTableList

    @classmethod
    def setup_class(cls):
        """Setup any state specific to the execution.

        Need to change qt2 to a QTable

        """
        super().setup_class()

        cls.qt2 = QTable(cls.qt2, copy=False)
        cls.inp["table2"] = cls.qt2
        cls.ordered_inp["table2"] = cls.qt2
        cls.QT["table2"] = cls.qt2

    # /def


# /classs


##############################################################################
# END

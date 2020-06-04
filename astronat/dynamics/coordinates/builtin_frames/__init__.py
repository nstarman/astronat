# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Coordinate frames inter-operable with Astropy frames.

Users shouldn't use this module directly, but rather import from the
`astropy.coordinates` module.  While it is likely to exist for the long-term,
the existence of this package and details of its organization should be
considered an implementation detail, and is not guaranteed to hold for future
versions.

Notes
-----
The builtin frame classes are all imported automatically into this package's
namespace, so there's no need to access the sub-modules directly.
To implement a new frame, a developer should add the frame as a new
module in this package.  Any "self" transformations (i.e., those that transform
from one frame to another frame of the same class) should be included in that
module.  Transformation functions connecting the new frame to other frames
should be in a separate module, which should be imported in this package's
``__init__.py`` to ensure the transformations are hooked up when this package
is imported.  Placing the transformation functions in separate modules avoids
circular dependencies, because they need references to the frame classes.

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["Astropy"]

__all__ = [
    # modules
    "orbitskyoffset",
    "orbitoffset",
    # functions
    "OrbitSkyOffsetFrame",
    "OrbitOffsetFrame",
]


##############################################################################
# IMPORTS

# THIRD PARTY

from astropy.coordinates.baseframe import frame_transform_graph


# PROJECT-SPECIFIC

from .orbitskyoffset import OrbitSkyOffsetFrame
from .orbitoffset import OrbitOffsetFrame

from . import orbitskyoffset, orbitoffset


# TODO: Astropy < 3.1.2 does not have make_transform_graph_docs
try:
    from astropy.coordinates.builtin_frames import make_transform_graph_docs
except ImportError:
    pass
else:
    _transform_graph_docs = make_transform_graph_docs(frame_transform_graph)

    # Here, we override the module docstring so that sphinx renders the
    # transform graph without the developer documentation in the main
    # docstring above.
    __doc__ = _transform_graph_docs

# /try


##############################################################################
# END

# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.
"""
Root-level database for an ORANGE input.
"""

from omnutils.db.entry import Database

from .geometry import geometry
from .universe import universes
from .universe.utils import CheckUniverseNames

root_db = Database(
    'root', "Top-level options for ORANGE input", [
        geometry,
        # arrays, # TODO
        universes,
        CheckUniverseNames('/geometry/global'),
    ])

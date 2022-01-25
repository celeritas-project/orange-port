# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.

from omnutils.db.command import UnzipMap
from omnutils.db.postprocessor import CheckListLengths
from omnutils.db.entry import (Parameter, Database, Deleted, )
from omnutils.db.validator import to_clean_string, to_bool, to_tol, to_pos
from omnutils.db.listvalidator import to_list_string, to_list_nonneg_int

geometry = Database(
    'geometry', "Global geometry options", [
        Parameter(
            'global', to_clean_string,
            "Name of the global universe",),
        Deleted(
            'deduplication_warning',
            "Surface deduplication no longer warns",
            short='dedupe_warn'),
        Parameter(
            'write_kdtree', to_bool,
            "Output the k-D tree representation for each universe",
            default=False,
            hidden=True),
        Deleted(
            'check_overlapping_volumes',
            "Overlapping volumes are now automatically checked"),
        Parameter(
            'tolerance', to_tol,
            "Global tolerance for geometry construction and particle bumping",
            short='tol',
            default=1.e-8),
        Parameter(
            'length_scale', to_pos,
            "Characteristic length scale of the problem",
            default=1.0),
        # Comp name -> matid declaration
        UnzipMap('comps', 'comp', 'matid'),
        Parameter(
            'composition', to_list_string,
            "Provide an optional explicit mapping for composition names",
            short='comp',
            default=[]),
        Parameter(
            'matid', to_list_nonneg_int,
            "Matids corresponding to the given commposition names",
            default=[]),
        CheckListLengths(('composition', 'matid')),
    ])

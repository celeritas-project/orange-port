# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.
"""
Universe components and universes.
"""

from omnutils.exceptions import OmniValueError, ParameterError
from omnutils.db.applicability import (Context,)
from omnutils.db.entry import (
    Database, Polymorphic, Sublist, Parameter, PostProcessor, make_default,
    make_postprocessor, empty_dict, empty_list, optional_entry, Deleted,
    Deprecated
)
from omnutils.db.validator import (
    Option, PathValidator, to_clean_string, to_pos_int, to_float,
    to_tol, to_nonneg_int, to_bool
)
from omnutils.db.listvalidator import (
    ListValidator, to_space_vector, not_empty, ZERO_VECTOR,
)

from .shapes import shape_list
from .utils import (
    to_comp_string,
    to_universe, to_list_sshapes, to_list_universe, to_dim_vector,
    delete_identity, boundary_command, euler_command, check_interior_shapes,
    check_cell_shapes, CheckUniverseNames,
    rotate_param, translate_param
)


def sum_less_than_unity(values):
    summed = sum(values)
    if summed >= 1.0:
        raise ParameterError("sum of volume fractions must be less than one, "
                             "but it is {summed:g}".format(**locals()))


to_list_vf = ListValidator(to_tol, (not_empty, sum_less_than_unity),
                           "fractions that sum to less than one")


class CheckArraySize(PostProcessor):
    __slots__ = ('params',)

    def __init__(self, params):
        self.params = params
        PostProcessor.__init__(self, description="Check array sizes")

    def __call__(self, stack):
        db = stack.current.out

        expected = 1
        for p in self.params:
            expected *= db[p]
        actual = len(db['fill'])

        if expected != actual:
            raise OmniValueError(
                "Given {}, expected {} fill  but got {}".format(
                    ",".join(self.params), expected, actual),
                db['fill']
            )


@make_default("'unit' if 'place' parameter is present")
def default_origin(stack):
    if 'place' in stack.current.inp:
        return 'unit'
    #elif 'origin' in stack.current.inp:
    #    return 'array'
    return 'array'


cell_list = Sublist(Database(
    'cell', "Cell definition", [
        Parameter(
            'name', to_clean_string,
            "Name of the cell",),
        Deleted('matid',
                "'matid' has been replaced with 'comp', the name of the "
                "composition"),
        Parameter(
            'composition', to_comp_string,
            "Name of the composition that fills this cell",
            short='comp'),
        Parameter(
            'volume', to_float,
            "Add a pre-calculated volume for this cell",
            default=optional_entry),
        Parameter(
            'shapes', to_list_sshapes,
            "Senses and shape names defining this cell"),
        check_cell_shapes,
    ]), default=empty_list)


hole_list = Sublist(Database(
    'hole', "Inserted sub-universes", [
        Parameter(
            'name', to_clean_string,
            "Name of the shape that this hole creates"),
        Parameter(
            'fill', to_universe,
            "Name of the universe to fill this hole with"),
        CheckUniverseNames('fill'),
        translate_param,
        euler_command,
        rotate_param,
        delete_identity,
    ]), default=optional_entry)


# TODO: implement placed arrays
# array_list = Sublist(Database(
#     'array', "Arrays placed into local shapes", [
#         Parameter(
#             'name', to_clean_string,
#             "Name of this array placement"),
#         Parameter(
#             'shapes', to_list_sshapes,
#             "Senses and shape names enclosing and truncating the array"),
#         Parameter(
#             'fill', to_universe,
#             "Name of the array universe to place"),
#         CheckUniverseNames('fill'),
#     ]), default=optional_entry)


unit_universe = Database(
    'universe', "Geometry unit (general universe)", [
        Parameter(
            'name', to_clean_string, "Name of the universe",),
        Parameter(
            'implicit_boundary', to_bool,
            "Boundary shapes automatically truncate internal shapes",
            default=False
        ),
        Parameter(
            'implicit_holes', to_bool,
            "Holes automatically truncate media and arrays",
            default=False
        ),
        Parameter(
            'otf_error_checking', to_bool,
            "Force on-the-fly checking of overlapping regions",
            default=False
        ),
        shape_list,
        hole_list,
        # TODO: placed_array_list,
        cell_list,
        boundary_command,
        Parameter(
            'interior', to_list_sshapes,
            "Senses and shape names defining the interior of this universe"),
        check_interior_shapes,
    ])


options_db = Database(
    'options', "Random universe construction options", [
        Parameter(
            'seed', to_nonneg_int,
            "Seed value for instantiating this universe",
            default=0),
        Parameter(
            'failure_batch_size', to_pos_int,
            "Number of samples per batch to test",
            default=1e5),
        Parameter(
            'failure_tolerance', to_tol,
            "Fraction of samples per batch that must be rejected to fail",
            default=(1.0 - 1e-5)),
        Parameter(
            'insert_method', Option(("naive",)),
            "Method for attempting to insert particles into a pebble",
            default="naive"),
    ])


random_universe = Database(
    'universe', "Randomly constructed universe", [
        Parameter(
            'name', to_clean_string, "Name of the universe",),
        shape_list.db,  # A single shape bounding the interior
        Parameter(
            'composition', to_comp_string,
            "Composition outside of the particles (the matrix)",
            short='comp'),
        Parameter(
            'fill', to_list_universe,
            "Names of particle universes to emplace"),
        CheckUniverseNames('fill'),
        Parameter(
            'volume_fraction', to_list_vf,
            "Names of particle universes to emplace",
            short='vf'),
        options_db,
    ])


_array_params = [
    Parameter(
        'name', to_clean_string,
        "Name of the universe",),
    Parameter(
        'fill', to_list_universe,
        "Names of array universe fills, indexed as ZUV"),
    CheckUniverseNames('fill'),

    # GG parameters (TODO: add hidden parameter so that adding these makes the
    # array "enclosed" (i.e. a separate general unit rather than an array to be
    # embedded)
    Parameter(
        'origin_is', Option(('array', 'unit')),
        "Whether ``origin`` is the array lower-left or a unit origin",
        default=default_origin),
    Parameter(
        'origin', to_space_vector,
        "Location of the lower-left corner of the array bounds",
        default=ZERO_VECTOR),
    Parameter(
        'place', to_dim_vector,
        "Center of this unit daughter",
        applicability=Context({'origin_is': ('unit',)})),
    Deleted('clip_composition',
        "Automatic material fill is no longer supported"),
    shape_list,
    hole_list,
    boundary_command,
    Parameter(
        'interior', to_list_sshapes,
        "Sense/shapes defining the interior of the array",
        default=optional_entry),  # Array can auto-set boundary
    check_interior_shapes,
]


rect_array = Database(
    'universe', "Rectangular array", _array_params + [
        Parameter(
            'nx', to_pos_int,
            "Number of units in the X direction",),
        Parameter(
            'ny', to_pos_int,
            "Number of units in the Y direction",),
        Parameter(
            'nz', to_pos_int,
            "Number of units in the Z direction",
            default=1),
        CheckArraySize(('nx', 'ny', 'nz')),
    ])


hex_array = Database(
    'universe', "Hexagonal array", _array_params + [
        Parameter(
            'layout', Option(('rhomb', 'rect')),
            "Layout of the hex array elements",
            default='rhomb'),
        Parameter(
            'nu', to_pos_int,
            "Number of units in the U direction",),
        Parameter(
            'nv', to_pos_int,
            "Number of units in the V direction",),
        Parameter(
            'nz', to_pos_int,
            "Number of units in the Z direction",
            default=1),
        CheckArraySize(('nu', 'nv', 'nz')),
    ])


dod_array = Database(
    'universe', "Dodecahedral array", _array_params + [
        Parameter(
            'nx', to_pos_int,
            "Number of units in the X direction",),
        Parameter(
            'ny', to_pos_int,
            "Number of units in the Y direction",),
        Parameter(
            'nz', to_pos_int,
            "Number of units in the Z direction",
            default=1),
        CheckArraySize(('nx', 'ny', 'nz')),
    ])


rtk_universe = Database(
    'universe', "Insert an RTK geometry from an external file", [
        Parameter(
            'name', to_clean_string,
            "Name of the universe",),
        Parameter(
            'input', PathValidator(mode="read", ext=".xml"),
            "Path to the RTK geometry XML file"),
        shape_list,
        hole_list,
        boundary_command,
        Parameter(
            'interior', to_list_sshapes,
            "Sense/shapes defining the interior of the array",),
        check_interior_shapes,
    ])


core_universe = Database(
    'universe',
    "Insert a VERA-defined reactor core", [
        Parameter('name', to_clean_string, "Name of universe",),
    ])


# Note: using a list here automatically converts to an OrderedDict
universes = Sublist(Polymorphic([
    ('unit', unit_universe),
    ('general', unit_universe), # deprecated
    # ('keno6', unit_universe), # removed: set implicit boundary parameters
    ('random', random_universe),
    ('array', rect_array),
    ('hexarray', hex_array),
    ('dodarray', dod_array),
    ('rtk', rtk_universe),
    ('core', core_universe),
], "Universes"))

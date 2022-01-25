# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.
"""
Utilities for universes and shapes.
"""
from six import string_types

from math import cos, sin, pi

from omnutils.location import str_node, get_location
from omnutils.exceptions import InvalidMultiOptionError

from omnutils.db.entry import (
    make_postprocessor, make_command, described, PostProcessor,
    Parameter, Deprecated)
from omnutils.db.validator import (
    to_str, to_clean_string, to_pos, to_float, to_nonneg_int, to_fraction,
    make_validator
)
from omnutils.db.listvalidator import (
    ListValidator, NumericListValidator, not_empty, monotonic_increasing,
    HasLength, to_space_vector, to_space_matrix, ZERO_VECTOR, IDENTITY_MATRIX)

@make_postprocessor("Squelch identity rotations and null translations",
                    hidden=True)
def delete_identity(stack):
    for params in (stack.current.out, stack.current.exp):
        if params['translate'] == ZERO_VECTOR:
            del params['translate']
        if params.get('rotate', None) == IDENTITY_MATRIX:
            del params['rotate']


rotate_param = Parameter(
    'rotate', to_space_matrix,
    "Rotation matrix",
    default=IDENTITY_MATRIX
)
translate_param = Parameter(
    'translate', to_space_vector,
    "Local-to-global translation vector",
    default=ZERO_VECTOR,
    short='origin',
)

@make_validator("Fraction, or ``---`` (in python, ``None``) for no rotation")
def to_optional_frac(value):
    if value is None:
        return 0.0
    return to_fraction(value)


to_euler = NumericListValidator(
    to_optional_frac, (HasLength(3),),
    description="3 values, each a fraction or ``---`` (i.e. ``None``)"
)


@make_command("euler",
              "Perform an Euler rotation about the :math:`z,x',z''` axes",
              ('rotate'),
              units="turns")
def euler_command(db, value):
    # Input value should be three rotations (alpha/beta/gamma)
    (a, b, c) = [v * 2 * pi for v in to_euler(value)]
    (s1, c1) = (sin(a), cos(a))
    (s2, c2) = (sin(b), cos(b))
    (s3, c3) = (sin(c), cos(c))

    # Result is rotation matrix
    db['rotate'] = [
        c1 * c3 - c2 * s1 * s3, -c1 * s3 - c2 * c3 * s1,  s1 * s2,
        c3 * s1 + c1 * c2 * s3,  c1 * c2 * c3 - s1 * s3, -c1 * s2,
        s2 * s3,                 c3 * s2,       c2]


@make_validator("name of another existing universe")
def to_universe(value):
   return to_clean_string(value)


@make_validator("name of a shape in the local universe")
def to_shape(value):
   return to_clean_string(value)


@described("non-empty string")
def to_comp_string(value):
    # Check that it's a string
    value = to_str(value)

    if not value:
        raise ParameterError("got empty string")

    return value


@described("shape with optional sense (``-``/``+``) or negation ``~``")
def to_sense_shape(value):
    """Ensure that the passed value is a valid shape name with optional sense.

    The output will have a single character that is its sense.
    """
    # Make sure value is a string
    shape = to_str(value)

    # Extract manually specified sense
    sense = shape[0]
    if sense in "-+~!":
        # Convert pos/negation/outside sense to '~'
        sense = '~' if sense in '+!~' else ''
        shape = shape[1:]
    else:
        # Default to 'inside'
        sense = ''

    # Validate shape name
    shape = to_shape(shape)

    return str_node(shape + sense, location=get_location(value))


to_list_sshapes = ListValidator(to_sense_shape, (not_empty,))

to_list_universe = ListValidator(to_universe, (not_empty,))

to_dim_vector = NumericListValidator(
    to_nonneg_int, (HasLength(3),),
    description="length-3 logical position vector"
)

to_extents = NumericListValidator(
    to_float, (HasLength(2), monotonic_increasing),
    description="(min, max) extent values",
)

to_pos_xy = NumericListValidator(to_pos, (HasLength(2),),
                                 description="length-2 (x,y) position"
                                 )
to_pos_xyz = NumericListValidator(to_pos, (HasLength(3),),
                                  description="3 positive floats")

_to_surf_list = ListValidator(to_clean_string)

@make_validator("list of shape surfaces (e.g. mx, co) or ``*`` for all")
def to_surf_list(value):
    if isinstance(value, string_types):
        value = [value]
    if len(value) == 1 and value[0] == "*":
        return value
    return _to_surf_list(value)

# Inside no longer requires '-' sense by default
boundary_command = Deprecated('boundary', 'interior')

def _shape_from_ss(s):
    if s[-1] == '~':
        return s[:-1]
    return s

def _check_shapes(stack, parameter, valid_paths):
    # Get provided shape/sense list
    given_sense_shapes = stack[parameter]
    # Convert to just shapes
    given_shapes = set(_shape_from_ss(s) for s in given_sense_shapes)

    # Create a set of all valid shapes
    valid_shapes = set()
    for path in valid_paths:
        try:
            these_shapes = stack.get(path)
        except KeyError:
            pass
        else:
            valid_shapes.update(these_shapes)

    invalid_shapes = given_shapes - valid_shapes
    if invalid_shapes:
        raise InvalidMultiOptionError(
            "shape", invalid_shapes, valid_shapes,
            location=get_location(given_sense_shapes[0]))


@make_postprocessor("Check that ``interior`` shapes have been defined")
def check_interior_shapes(stack):
    if 'interior' not in stack.current.out:
        # Array universe can be implicitly defined
        return

    _check_shapes(stack, 'interior', ['shape', 'hole'])


@make_postprocessor("Check that ``shapes`` have been defined")
def check_cell_shapes(stack):
    _check_shapes(stack, 'shapes', ['../../shape', '../../hole'])


class CheckShapeNames(PostProcessor):
    """Validate the shape names at 'param' against items in 'shape_lists'

    'shape_lists' can be a relative path such as ``../shape`` or a list of
    paths such as ``['shape','hole']``.
    """

    __slots__ = ('parameter', 'shape_list')

    def __init__(self, parameter, shape_list):
        self.parameter = parameter
        self.shape_list = shape_list

        desc = ("Shape names in ``{parameter}`` must be in the shape list at "
                "``{shape_list}``".format(**locals()))
        PostProcessor.__init__(self, desc)

    def __call__(self, stack):
        # Get provided shape/sense list
        given_sense_shapes = stack[self.parameter]
        # Convert to just shapes
        given_shapes = set(_shape_from_ss(s) for s in given_sense_shapes)

        # Get the dict of valid shapes
        valid_shapes = stack[self.shape_list]
        # Extract just the shape names
        valid_shapes = set(valid_shapes)

        invalid_shapes = given_shapes - valid_shapes

        if invalid_shapes:
            raise InvalidMultiOptionError(
                "shape", invalid_shapes, valid_shapes,
                location=get_location(given_sense_shapes[0]))


class CheckUniverseNames(PostProcessor):
    __slots__ = ('parameter',)

    def __init__(self, parameter):
        self.parameter = parameter

        desc = ("Universe names in ``{parameter}`` must already have been "
                "defined".format(**locals()))
        PostProcessor.__init__(self, desc)

    def __call__(self, stack):
        # Get provided universe
        given_universes = stack[self.parameter]
        if isinstance(given_universes, string_types):
            # Single universe
            given_universes = [given_universes]
        given = set(given_universes)

        valid_universes = stack['/universe']
        valid_universes = set(valid_universes)

        invalid = given - valid_universes

        if invalid:
            raise InvalidMultiOptionError(
                "universe", invalid, valid_universes,
                location=get_location(given_universes[0]))

# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.
"""
Create list of shapes for use in general universe.
"""

from collections import OrderedDict
from copy import deepcopy

from omnutils.db.entry import (
    Database, Polymorphic, Sublist, Parameter, PolymorphicBuilder,
    empty_list, described, make_command
)
from omnutils.db.command import ShorthandParameters
from omnutils.db.postprocessor import CheckLessThan
from omnutils.db.validator import (
    ApplyAndTestValidator, to_clean_string, to_float, to_axis, to_int,
    to_fraction, to_pos, to_bool
)
from omnutils.db.listvalidator import (
    NumericListValidator, HasLength, to_space_vector, ZERO_VECTOR,
)

from .utils import (
    to_dim_vector, to_extents, to_pos_xy, to_pos_xyz, to_surf_list,
    delete_identity, euler_command, translate_param, rotate_param
)


to_ppiped_angle = ApplyAndTestValidator(float, (lambda v: 0 <= v < 0.25),
                                        "real number in [0,.25)", "value outside range")


# Hack to allow shape to get the right extended documentation for euler_command
euler_command = deepcopy(euler_command)

shape_builder = PolymorphicBuilder(
    'shape', "Geometry shapes",
    prepend=[
        Parameter(
            'name', to_clean_string,
            "Name of the shape",),
        translate_param,
        euler_command,
        rotate_param,
        delete_identity,
        Parameter(
            'reflect', to_surf_list,
            "Reflecting boundary faces",
            default=[]),
    ],
    append=[])

shape_builder.add("cuboid", "Box shape", [
    ShorthandParameters('faces',
                            "xmin xmax ymin ymax zmin zmax".split()),
    Parameter(
        'xmin', to_float,
        "Minimum x-coordinate of box source"),
    Parameter(
        'xmax', to_float,
        "Maximum x-coordinate of box source"),
    CheckLessThan('xmin', 'xmax'),
    Parameter(
        'ymin', to_float,
        "Minimum y-coordinate of box source"),
    Parameter(
        'ymax', to_float,
        "Maximum y-coordinate of box source"),
    CheckLessThan('ymin', 'ymax'),
    Parameter(
        'zmin', to_float,
        "Minimum z-coordinate of box source"),
    Parameter(
        'zmax', to_float,
        "Maximum z-coordinate of box source"),
    CheckLessThan('zmin', 'zmax'),
])
shape_builder.alias("cuboid", "box")

shape_builder.add("sphere", "Sphere shape", [
    Parameter(
        'radius', to_float,
        "Radius of sphere",
        short="r"),
])

shape_builder.add("cyl", "Cylinder shape", [
    Parameter(
        'axis', to_axis,
        "Axis along the cylinder"),
    Parameter(
        'extents', to_extents,
        "Negative and positive position along the given axis"),
    Parameter(
        'radius', to_float,
        "Radius of cylinder",
        short="r"),
])

to_arc_angle = ApplyAndTestValidator(float, (lambda v: v >= 0.0 and v <= 0.5),
                                     "real number in [0.0, 0.5]", "value outside range")

shape_builder.add("cylsegment", "Cylinder segment shape", [
    Parameter(
        'extents', to_extents,
        "Negative and positive position along the z-axis"),
    Parameter(
        'inner_radius', to_pos,
        "Inner radius",
        short="ri"),
    Parameter(
        'outer_radius', to_pos,
        "Outer radius",
        short="ro"),
    CheckLessThan('inner_radius', 'outer_radius'),
    Parameter(
        'angle', to_fraction,
        "Beginning angle CCW from :math:`x = 0`",
        short="a",
        units="turns"),
    Parameter(
        'arc', to_arc_angle,
        "Angle subtended by cylindrical segment",
        short="da",
        units="turns"),
])
shape_builder.alias("cylsegment", "pad")

shape_builder.add("ring", "Cylindrical shell shape", [
    Parameter(
        'inner_radius', to_float,
        "Inner radius",
        short="ri"),
    Parameter(
        'outer_radius', to_float,
        "Outer radius",
        short="ro"),
    CheckLessThan('inner_radius', 'outer_radius'),
    Parameter(
        'extents', to_extents,
        "Negative and positive position along the Z axis"),
])
shape_builder.alias("ring", "cylshell")

shape_builder.add("prism", "Regular prism shape", [
    Parameter(
        'num_sides',
        ApplyAndTestValidator(to_int, (lambda v: v >= 3),
                              "integer >= 3", "integer < 3"),
        "Number of sides on the prism"),
    Parameter(
        'apothem', to_float,
        "Inner radius of the prism",
        short="r"),
    Parameter(
        'extents', to_extents,
        "Negative and positive position along the Z axis"),
    Parameter(
        'rotfrac', to_fraction,
        "Angle (fraction of the angle spanned by one face) to rotate",
        default=0.0),
])

shape_builder.add("slab", "Infinite slab shape", [
    Parameter(
        'axis', to_axis,
        "Axis along the slab"),
    Parameter(
        'extents', to_extents,
        "Negative and positive position along the slab axis"),
])

shape_builder.add("plane", "Infinite half-space", [
    Parameter(
        'normal', to_space_vector,
        "Vector (possibly unnormalized) point outward from the plane"),
    Parameter(
        'point', to_space_vector,
        "A point somewhere on the plane"),
])

shape_builder.add("wedge", "Wedge shape", [
    Parameter(
        'corner_pt', to_pos_xy,
        "XY location of the right-angle corner of the wedge",
        short='xy'),
    Parameter(
        'height', to_pos,
        "Height of the wedge (along Z axis)",
        short='zlng'),
    Parameter(
        'width', to_pos,
        "Length of the hypotenuse of the wedge (along the X axis)",
        short='xbase'),
])

shape_builder.add("cone", "Cone shape", [
    Parameter(
        'axis', to_axis,
        "Axis along the centerline of the cone"),
    Parameter(
        'radii',
        NumericListValidator(to_pos, (HasLength(2),),
                             description="2 positive floats"),
        "Radii at the top and bottom",
        short='r'),
    Parameter(
        'extents', to_extents,
        "Negative and positive base positions along the given axis"),
])

shape_builder.add("ellipsoid", "Ellipsoid shape", [
    Parameter(
        'radii', to_pos_xyz,
        "Radii in the x, y, and z directions",
        short='r'),
])

shape_builder.add("hopper", "Hopper shape", [
    Parameter(
        'extents', to_extents,
        "Negative and positive base positions along the Z axis"),
    Parameter(
        'lower_pt', to_pos_xy,
        "X and Y half-lengths on the low side of the hopper",
        short='lo'),
    Parameter(
        'upper_pt', to_pos_xy,
        "X and Y half-lengths on the high side of the hopper",
        short='hi'),
])

shape_builder.add("righttet", "Right tetrahedron shape", [
    Parameter(
        'lengths', to_pos_xyz,
        "Length of the tetrahedron edges along the x, y, and z axes"),
])

shape_builder.add("triprism", "Triangular prism shape", [
    Parameter(
        'lengths', to_pos_xyz,
        "Length of the triangular prism's bounding box"),
])

shape_builder.add("ecylinder", "Elliptical cylinder shape", [
    Parameter(
        'radii', to_pos_xy,
        "Radii in the x and y directions",
        short='r'),
    Parameter(
        'extents', to_extents,
        "Negative and positive position along the given axis"),
])

shape_builder.add("rhombdod", "Rhombic dodecahedron shape", [
    Parameter(
        'apothem', to_pos,
        "Inner radius of the dodecahedron",
        short='r')
])

shape_builder.add("ppiped", "Parallelepiped", [
    Parameter(
        'lengths', to_pos_xyz,
        "Length of the x, xy, and xyz edges"),
    Parameter(
        'y_angle', to_ppiped_angle,
        "Angle between xy edge and y axis (psi)",
        short='psi', rename='psi',
        units="turns"),
    Parameter(
        'xy_angle', to_ppiped_angle,
        "Angle between xy component of xyz edge and x axis (phi)",
        short='phi', rename='phi',
        units="turns"),
    Parameter(
        'xyz_angle', to_ppiped_angle,
        "Angle between xyz edge and z axis (theta)",
        short='theta', rename='theta',
        units="turns"),
])

shape_builder.add("quadric", "General quadric", [
    Parameter(
        'second', to_space_vector,
        "Second-order coefficients :math:`aX^2 + bY^2 + cZ^2`",
        default=ZERO_VECTOR,
        short='abc'),
    Parameter(
        'cross', to_space_vector,
        "Second-order cross coefficients :math:`dXY + eYZ + fXZ`",
        default=ZERO_VECTOR,
        short='def'),
    Parameter(
        'first', to_space_vector,
        "First-order coefficients :math:`gX + hY + iZ`",
        default=ZERO_VECTOR,
        short='ghi'),
    Parameter(
        'scalar', to_float,
        "Scalar coefficient :math:`j`",
        default=0,
        short='j'),
])

shape_list = Sublist(shape_builder.build(), default=empty_list)

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ParallelepipedShape.hh
 * \brief ParallelepipedShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A cuboid skewed by 3 angles
 *
 * is a shape with six faces composed of parallelograms, whose opposing faces
 * are parallel.. It is defined by specifying the length of the faces in the X,
 * Y, and Z axes, and the angles between the X-face and Y-axis, x' plane, the
 * angle between the Y-face and the Z-axis, y' plane, and the angle between x'
 * plane and the X-axis. The bottom face is an XY plane at Z = 0. The top face
 * is an XY plane intersecting Z at Z <= z axis length
 *
 * The x edge is along the x axis, the xy edge is inside the XY plane, and the
 * xyz edge is inside the entire positive octant. When the three angles are
 * zero, the edges are along the x, y, and z axes respectively.
 */
/*!
 */
class ParallelepipedShape final : public Shape
{
    using Base = Shape;

  private:
    //// DATA ////
    //! Scaled basis vectors of the ppiped
    SpaceMatrix basis_;

  public:
    // Construct with different lengths and internal angles
    ParallelepipedShape(Real3     edge_length,
                        real_type psi,
                        real_type theta,
                        real_type phi);

    //// DERIVED INTERFACE ////

    // Shape class name
    std::string type() const;

    // Whether the shape is convex (no internal surface crossings)
    bool is_convex() const;

    // Shape volume
    real_type volume() const;

    // Largest sphere radius that fits in this shape
    real_type inradius() const;

    // Build the shape's surfaces
    void build(ShapeBuilder& builder) const;

  private:
    Real3 calc_upper_vertex() const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

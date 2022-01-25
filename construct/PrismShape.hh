//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/PrismShape.hh
 * \brief PrismShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Regular polygon extruded along the Z axis.
 */
class PrismShape final : public Shape
{
    using Base = Shape;

  public:
    // Constructor
    PrismShape(unsigned int num_sides,
               real_type    apothem,
               real_type    rotate,
               real_type    lo,
               real_type    hi);

    //// ACCESSORS ////

    //! Number of sides
    unsigned int num_sides() const { return num_sides_; }

    //! Apothem
    real_type apothem() const { return apothem_; }

    //! Rotation offset (0 for bottom face at -Y, 1 for rotated 1 side)
    real_type rotation() const { return rotate_offset_; }

    //! Lower point
    real_type lower() const { return lo_; }

    //! Upper point
    real_type upper() const { return hi_; }

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
    //// DATA ////

    // Number of sides
    unsigned int num_sides_;

    // Distance from center to midpoint of its side
    real_type apothem_;

    // Rotational offset (0 has bottom face at -Y, 1 is congruent)
    real_type rotate_offset_;

    // Extents along the z axis
    real_type lo_;
    real_type hi_;

    //// PRIVATE ACCESSORS ////

    // Calculate the circumradius
    real_type calc_circumradius() const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/CylinderSegmentShape.hh
 * \brief CylinderSegmentShape class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Axis-aligned cylinder with radial segment boundaries and finite
 * extents.
 *
 * Cylindrical segments are z-axis aligned segments as shown below:
 *
 * \image html cylinder-segment.png "Cylindrical segment."
 *
 * The vectors bounding the segment, noting that \f$\theta = \pi/2\f$
 * (\f$\sin\theta = 1\f$ and \f$\cos\theta = 0\f$), are
 * \f[
 * \begin{array}{ll}
 * \mathbf{v}  &= \cos\phi_\mbox{begin}\mathrm{i} +
 *                \sin\phi_\mbox{begin}\mathrm{j}\\
 * \mathbf{v}' &= \cos\phi_\mbox{end}\mathrm{i} +
 *                \sin\phi_\mbox{end}\mathrm{j}
 * \end{array}
 * \f]
 * Then, the \b outward normals are defined:
 * \f[
 * \begin{array}{ll}
 *  \mathrm{n}_\mbox{begin} &= v_y\mathrm{i} - v_x\mathrm{j}\\
 *  \mathrm{n}_\mbox{end} &= -v'_y\mathrm{i} + v'_x\mathrm{j}
 * \end{array}
 * \f]
 * This yields the correct result:
 * \f[
 * \begin{array}
 *  \mathrm{v}\cdot\mathrm{n}_\mbox{begin} &= v_xv_y - v_xv_y = 0 \\
 *  \mathrm{v}'\cdot\mathrm{n}_\mbox{end} &= -v'_xv'_y + v'_xv'_y = 0
 * \end{array}
 * \f]
 *
 * We define the normals to be \b outward facing because the inside \c SENSE
 * for each plane is defined \c neg with respect to the normal.
 */
class CylinderSegmentShape final : public Shape
{
    using Base = Shape;

  public:
    // Constructors
    CylinderSegmentShape(real_type inner,
                         real_type outer,
                         real_type beg_rad,
                         real_type angle,
                         real_type lo,
                         real_type hi);

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
    // Cylinder radii
    real_type inner_;
    real_type outer_;

    // Extents along the given axis
    real_type lo_;
    real_type hi_;

    // Beginning and ending azimuthal angle (radians).
    real_type beg_rad_;
    real_type end_rad_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

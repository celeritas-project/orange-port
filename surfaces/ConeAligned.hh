//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/ConeAligned.hh
 * \brief ConeAligned class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class SimpleQuadric;
class GeneralQuadric;

//---------------------------------------------------------------------------//
/*!
 * Axis-aligned real_type-sheeted cone.
 *
 * For a cone parallel to the x axis:
 * \f[
    (y - y_0)^2 + (z - z_0)^2 - t^2 (x - x_0)^2 = 0
   \f]
 */
template<Axis T>
class ConeAligned
{
  public:
    //// ATTRIBUTES ////

    // Get the surface type of this surface.
    inline static constexpr SurfaceType surface_type();

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 4; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 2; }

    // Whether the given quadric can be simplified to this shape
    static bool can_simplify(const SimpleQuadric& sq);

    //@{
    //! Perpendicular axes
    static constexpr Axis U = (T == Axis::X ? Axis::Y : Axis::X);
    static constexpr Axis V = (T == Axis::Z ? Axis::Y : Axis::Z);
    //@}

  public:
    //// CONSTRUCTION ////

    // Construct with origin and tangent angle of opening
    ConeAligned(const Real3& origin, real_type tangent);

    // Construct from a degenerate SQ
    explicit ConeAligned(const SimpleQuadric& sq);

    // Inline construction from flattened coefficients (copying data)
    explicit inline ConeAligned(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return origin_.data(); }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new cone translated by some vector
    ConeAligned<T> translated(const Transform& t) const;

    // Return a cone transformed by a rotate/translate
    GeneralQuadric transformed(const Transform& t) const;

    // Clip a bounding box to a shape with this sense
    void clip(Sense sense, BoundingBox& bbox) const;

    //// CALCULATION ////

    // Determine the sense of the position relative to this surface (init)
    inline SignedSense calc_sense(const Real3& pos) const;

    // Determine distance to intersection
    inline void calc_intersections(const Real3& pos,
                                   const Real3& dir,
                                   bool         on_surface,
                                   real_type*   dist_iter) const;

    // Calculate outward normal at a position
    inline Real3 calc_normal(const Real3& pos) const;

    //// ACCESSORS ////

    //! Get the origin position along the normal axis
    Real3 origin() const { return origin_; }

    //! Get the square of the tangent
    real_type tangent_sq() const { return tsq_; }

  private:
    static constexpr int t_index() { return static_cast<int>(T); }
    static constexpr int u_index() { return static_cast<int>(U); }
    static constexpr int v_index() { return static_cast<int>(V); }

    //// DATA ////

    // Location of the vanishing point
    Real3 origin_;

    // Quadric value
    real_type tsq_;
};

// Print to a stream
template<Axis T>
std::ostream& operator<<(std::ostream&, const ConeAligned<T>&);

//---------------------------------------------------------------------------//
// TYPEDEFS
//---------------------------------------------------------------------------//

using ConeX = ConeAligned<Axis::X>;
using ConeY = ConeAligned<Axis::Y>;
using ConeZ = ConeAligned<Axis::Z>;

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "ConeAligned.i.hh"
//---------------------------------------------------------------------------//

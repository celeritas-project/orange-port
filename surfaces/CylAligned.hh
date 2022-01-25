//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylAligned.hh
 * \brief CylAligned class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class SimpleQuadric;
class GeneralQuadric;
template<Axis T>
class CylCentered;

//---------------------------------------------------------------------------//
/*!
 * Axis-aligned cylinder.
 *
 * The cylinder is centered about the template parameter Axis.
 *
 * For a cylinder along the x axis:
 * \f[
    (y - y_0)^2 + (z - z_0)^2 - R^2 = 0
   \f]
 */
template<Axis T>
class CylAligned
{
  public:
    //// ATTRIBUTES ////

    // Get the surface type of this surface.
    inline static constexpr SurfaceType surface_type();

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 3; }

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

    // Construct with origin and radius
    explicit CylAligned(Real3 origin, real_type radius);

    // Construct from a degenerate quadric
    explicit CylAligned(const SimpleQuadric& sq);

    // Construct from a centered cylinder
    explicit CylAligned(const CylCentered<T>& coc);

    // Inline construction from flattened coefficients (copying data)
    explicit inline CylAligned(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return &origin_u_; }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new cylinder translated by some vector
    CylAligned<T> translated(const Transform& t) const;

    // Return a cylinder transformed by a rotate/translate
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

    //! Get the origin vector along the 'u' axis
    real_type origin_u() const { return origin_u_; }

    //! Get the origin vector along the 'v' axis
    real_type origin_v() const { return origin_v_; }

    //! Get the square of the radius
    real_type radius_sq() const { return radius_sq_; }

  private:
    static constexpr int t_index() { return static_cast<int>(T); }
    static constexpr int u_index() { return static_cast<int>(U); }
    static constexpr int v_index() { return static_cast<int>(V); }

    //// DATA ////

    // Off-axis location
    real_type origin_u_;
    real_type origin_v_;

    // Square of the radius
    real_type radius_sq_;
};

// Print to a stream
template<Axis T>
std::ostream& operator<<(std::ostream&, const CylAligned<T>&);

//---------------------------------------------------------------------------//
// TYPEDEFS
//---------------------------------------------------------------------------//

using CylX = CylAligned<Axis::X>;
using CylY = CylAligned<Axis::Y>;
using CylZ = CylAligned<Axis::Z>;

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CylAligned.i.hh"
//---------------------------------------------------------------------------//

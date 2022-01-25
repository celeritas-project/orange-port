//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/Plane.hh
 * \brief Plane class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class SimpleQuadric;

//---------------------------------------------------------------------------//
/*!
 * Arbitrary defined plane
 *
 * A plane is a first-order quadric that satisfies \f[
    ax + bx + cz - d = 0
    \f]
 */
class Plane
{
  public:
    //// ATTRIBUTES ////

    //! Get the surface type of this surface.
    static constexpr SurfaceType surface_type() { return SurfaceType::p; }

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 4; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 1; }

    // Whether the given quadric can be simplified to this shape
    static bool can_simplify(const SimpleQuadric& sq);

  public:
    //// CONSTRUCTION ////

    // Construct with unnormalized normal and point on the plane
    Plane(const Real3& n, const Real3& p);

    // Construct with normalized normal and displacement
    Plane(const Real3& n, real_type d);

    // Construct from a degenerate SQ
    explicit Plane(const SimpleQuadric& sq);

    // Inline construction from flattened coefficients (copying data)
    explicit inline Plane(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return normal_.data(); }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new plane translated by some vector
    Plane translated(const Transform& t) const;

    // Return a plane transformed by a rotate/translate
    Plane transformed(const Transform& t) const;

    //! Clip a bounding box to a shape with this sense (null-op)
    void clip(Sense, BoundingBox&) const {}

    //// CALCULATION ////

    // Determine the sense of the position relative to this surface (init)
    inline SignedSense calc_sense(const Real3& x) const;

    // Determine distance to intersection
    inline void calc_intersections(const Real3& pos,
                                   const Real3& dir,
                                   bool         on_surface,
                                   real_type*   dist_iter) const;

    //! Calculate outward normal at a position
    Real3 calc_normal(const Real3&) const { return normal_; }

    //// ACCESSORS ////

    //! Normal to the plane
    Real3 normal() const { return normal_; }

    //! Distance from the origin along the normal to the plane
    real_type displacement() const { return d_; }

  private:
    //// DATA ////

    // Normal to plane (a,b,c)
    Real3 normal_;

    // n \dot P (d)
    real_type d_;
};

// Print to a stream
std::ostream& operator<<(std::ostream&, const Plane&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Plane.i.hh"
//---------------------------------------------------------------------------//

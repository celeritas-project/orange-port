//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SimpleQuadric.hh
 * \brief SimpleQuadric class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class GeneralQuadric;

//---------------------------------------------------------------------------//
/*!
 * Axis-aligned quadric
 *
 * Input:
 * \f[
   a(x - x_0)^2 + b(y - y_0)^2 + c(z - z_0)^2
   + 2d(x - x_0) + 2e(y - y_0) + 2f(z - z_0)
   + g = 0
  \f]
 *
 * Stored:
 * \f[
   ax^2 + by^2 + cz^2 + dx + ey + fz + g = 0
  \f]
 *
 * This can represent hyperboloids, ellipsoids, elliptical cylinders, etc.
 */
class SimpleQuadric
{
  public:
    //// ATTRIBUTES ////

    //! Get the surface type of this surface.
    static constexpr SurfaceType surface_type() { return SurfaceType::sq; }

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 7; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 2; }

    // Whether the given quadric can be simplified to this shape
    static bool can_simplify(const GeneralQuadric& gq);

  public:
    //// CONSTRUCTION ////

    // Construct at the origin
    SimpleQuadric(const Real3& abc, const Real3& def, real_type g);

    // Construct with all quadric coordinates
    SimpleQuadric(const Real3& abc,
                  const Real3& def,
                  real_type    g,
                  const Real3& origin);

    // Construct from a degenerate GQ
    explicit SimpleQuadric(const GeneralQuadric& gq);

    // Inline construction from flattened coefficients (copying data)
    explicit inline SimpleQuadric(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return &a_; }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new quadric translated by some vector
    SimpleQuadric translated(const Transform& t) const;

    // Return a quadric transformed by a rotate/translate
    GeneralQuadric transformed(const Transform& t) const;

    //! Clip a bounding box to a shape with this sense (null-op)
    void clip(Sense, BoundingBox&) const {}

    //// CALCULATION ////

    // Determine the sense of the position relative to this surface (init)
    inline SignedSense calc_sense(const Real3& pos) const;

    // Determine distance to intersection
    inline void calc_intersections(const Real3& pos,
                                   const Real3& dir,
                                   bool         on_surface,
                                   real_type*   dist_iter) const;

    //! Calculate outward normal at a position
    inline Real3 calc_normal(const Real3& pos) const;

    //// ACCESSORS ////

    //! Second-order terms
    const Real3& second() const { return {&a_}; }

    //! First-order terms
    const Real3& first() const { return {&d_}; }

    //! Zeroth-order term
    real_type zeroth() const { return g_; }

  private:
    //// DATA ////

    // Second-order terms (a, b, c)
    real_type a_, b_, c_;
    // First-order terms (d, e, f)
    real_type d_, e_, f_;
    // Constant term (g)
    real_type g_;
};

// Print to a stream
std::ostream& operator<<(std::ostream&, const SimpleQuadric&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "SimpleQuadric.i.hh"
//---------------------------------------------------------------------------//

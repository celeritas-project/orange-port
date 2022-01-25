//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/GeneralQuadric.hh
 * \brief GeneralQuadric class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * General quadric surface.
 *
 * \f[
    ax^2 + by^2 + cz^2 + dxy + eyz + fzx + gx + hy + iz + j = 0
   \f]
 */
class GeneralQuadric
{
  public:
    //// ATTRIBUTES ////

    //! Get the surface type of this surface.
    static constexpr SurfaceType surface_type() { return SurfaceType::gq; }

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 10; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 2; }

  public:
    //// CONSTRUCTION ////

    // Construct with all quadric coordinates (second, cross, first, zeroth)
    GeneralQuadric(const Real3& abc,
                   const Real3& def,
                   const Real3& ghi,
                   real_type    j);

    // Inline construction from flattened coefficients (copying data)
    explicit inline GeneralQuadric(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return &a_; }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new quadric translated by some vector
    GeneralQuadric translated(const Transform& t) const;

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

    // Calculate outward normal at a position
    inline Real3 calc_normal(const Real3& pos) const;

    //// ACCESSORS ////

    //! Second-order terms
    const Real3& second() const { return const Real3 & (&a_); }

    //! Cross terms (xy, yz, zx)
    const Real3& cross() const { return const Real3 & (&d_); }

    //! First-order terms
    const Real3& first() const { return const Real3 & (&g_); }

    //! Zeroth-order term
    real_type zeroth() const { return j_; }

  private:
    //// DATA ////

    // Second-order terms (a, b, c)
    real_type a_, b_, c_;
    // Second-order cross terms (d, e, f)
    real_type d_, e_, f_;
    // First-order terms (g, h, i)
    real_type g_, h_, i_;
    // Constant term
    real_type j_;
};

// Print to a stream
std::ostream& operator<<(std::ostream&, const GeneralQuadric&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "GeneralQuadric.i.hh"
//---------------------------------------------------------------------------//

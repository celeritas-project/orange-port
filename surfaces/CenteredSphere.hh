//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CenteredSphere.hh
 * \brief CenteredCenteredSphere class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class Sphere;

//---------------------------------------------------------------------------//
/*!
 * CenteredSphere at an arbitrary spatial point.
 */
class CenteredSphere
{
  public:
    //// ATTRIBUTES ////

    //! Get the surface type of this surface.
    static constexpr SurfaceType surface_type() { return SurfaceType::so; }

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 1; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 2; }

    // Whether the given sphere is at the origin
    static bool can_simplify(const Sphere& sq);

  public:
    //// CONSTRUCTION ////

    // Construct with radius
    explicit CenteredSphere(real_type radius);

    // Inline construction from flattened coefficients (copying data)
    explicit inline CenteredSphere(const real_type* coeff);

    // Construct from a degenerate SQ
    explicit CenteredSphere(const Sphere& sq);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return &radius_sq_; }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new sphere translated by some vector
    Sphere translated(const Transform& t) const;

    //! Return a sphere transformed by a rotate/translate
    Sphere transformed(const Transform& t) const;

    // Clip a bounding box to a shape with this sense
    void clip(Sense sense, BoundingBox& bbox) const;

    //// CALCULATION ////

    // Determine the sense of the position relative to this surface (init)
    inline SignedSense calc_sense(const Real3& x) const;

    // Determine distance to intersection
    inline void calc_intersections(const Real3& pos,
                                   const Real3& dir,
                                   bool         on_surface,
                                   real_type*   dist_iter) const;

    // Calculate outward normal at a position
    inline Real3 calc_normal(const Real3& pos) const;

    //// ACCESSORS ////

    //! Get the square of the radius
    real_type radius_sq() const { return radius_sq_; }

  private:
    //// DATA ////

    // Square of the radius
    real_type radius_sq_;
};

// Print to a stream
std::ostream& operator<<(std::ostream&, const CenteredSphere&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CenteredSphere.i.hh"
//---------------------------------------------------------------------------//

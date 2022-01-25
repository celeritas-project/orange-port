//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/PlaneAligned.hh
 * \brief PlaneAligned class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class Plane;

//---------------------------------------------------------------------------//
/*!
 * Axis-aligned plane
 */
template<Axis T>
class PlaneAligned
{
  public:
    //// ATTRIBUTES ////

    // Get the surface type of this surface.
    inline static constexpr SurfaceType surface_type();

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 1; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 1; }

    // Whether the given plane can be simplified to this surface
    static bool can_simplify(const Plane& p);

  public:
    //// CONSTRUCTION ////

    // Construct with location along the axis
    explicit PlaneAligned(real_type position);

    // Construct from a degenerate plane
    explicit PlaneAligned(const Plane& p);

    // Inline construction from flattened coefficients (copying data)
    explicit inline PlaneAligned(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return &position_; }

    //! Access the type-deleted surface data
    GenericSurfaceRef view() const { return GenericSurfaceRef::from_surface(*this); }

    //// TRANSFORMATION ////

    // Return a new plane translated by some vector
    PlaneAligned<T> translated(const Transform& t) const;

    // Return a plane transformed by a rotate/translate
    Plane transformed(const Transform& t) const;

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

    // Position along the axis
    real_type position() const { return position_; }

  private:
    static constexpr int t_index() { return static_cast<int>(T); }

    //// DATA ////

    // Position along the axis
    real_type position_;
};

// Print to a stream
template<Axis T>
std::ostream& operator<<(std::ostream&, const PlaneAligned<T>&);

//---------------------------------------------------------------------------//
// TYPEDEFS
//---------------------------------------------------------------------------//

using PlaneX = PlaneAligned<Axis::X>;
using PlaneY = PlaneAligned<Axis::Y>;
using PlaneZ = PlaneAligned<Axis::Z>;

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "PlaneAligned.i.hh"
//---------------------------------------------------------------------------//

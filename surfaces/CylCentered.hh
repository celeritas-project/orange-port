//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/CylCentered.hh
 * \brief CylCentered class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "Definitions.hh"

namespace celeritas
{
class GeneralQuadric;
template<Axis T>
class CylAligned;

//---------------------------------------------------------------------------//
/*!
 * Axis-aligned cylinder centered about the origin.
 *
 * The cylinder is centered about the template parameter Axis.
 *
 * For a cylinder along the x axis:
 * \f[
    y^2 + z^2 - R^2 = 0
   \f]
 *
 * This is an optimization of the CylAligned. The motivations are:
 * - Many geometries have units with concentric cylinders centered about the
 *   origin, so having this as a special case reduces the memory usage of those
 *   units (improved memory localization).
 * - The cylindrical mesh geometry has lots of these cylinders, so efficient
 *   tracking through its cells should make this optimization worthwhile.
 */
template<Axis T>
class CylCentered
{
  public:
    //@{
    //! SurfaceType aliases
    using size_type = size_type;
    //@}

  public:
    //// ATTRIBUTES ////

    //! Get the surface type of this surface.
    static inline constexpr SurfaceType surface_type();

    //! Number of values associated with data(), and pointer construction
    static constexpr size_type size() { return 1; }

    //! Number of possible intersections
    static constexpr size_type num_intersections() { return 2; }

    // Whether the given ortho cyl can be simplified to this shape
    static bool can_simplify(const CylAligned<T>& oc);

    //@{
    //! Perpendicular axes
    static constexpr Axis U = (T == Axis::X ? Axis::Y : Axis::X);
    static constexpr Axis V = (T == Axis::Z ? Axis::Y : Axis::Z);
    //@}

  public:
    //// CONSTRUCTION ////

    // Construct with origin and radius
    explicit CylCentered(real_type radius);

    // Construct from a degenerate ortho cyl
    explicit CylCentered(const CylAligned<T>& ortho_cyl);

    // Inline construction from flattened coefficients (copying data)
    explicit inline CylCentered(const real_type* coeff);

    //! Access the data for this object, for inlining into a surface
    const real_type* data() const { return &radius_sq_; }

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

    //! Get the square of the radius
    real_type radius_sq() const { return radius_sq_; }

  private:
    static constexpr int t_index() { return static_cast<int>(T); }
    static constexpr int u_index() { return static_cast<int>(U); }
    static constexpr int v_index() { return static_cast<int>(V); }

    //// DATA ////

    // Square of the radius
    real_type radius_sq_;
};

// Print to a stream
template<Axis T>
std::ostream& operator<<(std::ostream&, const CylCentered<T>&);

//---------------------------------------------------------------------------//
// TYPEDEFS
//---------------------------------------------------------------------------//

using CCylX = CylCentered<Axis::X>;
using CCylY = CylCentered<Axis::Y>;
using CCylZ = CylCentered<Axis::Z>;

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CylCentered.i.hh"
//---------------------------------------------------------------------------//

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/Definitions.hh
 * \brief Class definitions used by surfacesers
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"
#include "Nemesis/containers/Span.hh"
#include "../Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// ENUMS
//---------------------------------------------------------------------------//
/*!
 * The logical state of a point with respect to a surface.
 *
 * For a plane, "outside" is equivalent to
 * \f[
   \vec x \cdot \vec n > 0
 * \f]
 * and "inside" is to the left of the plane's normal (a negative dot product).
 * The exact equality to zero is literally an "edge case" but it can happen
 * with inter-universe coincident surfaces as well as carefully placed
 * particle sources and ray tracing.
 *
 * As an implementataion detail, the "on" case is currently *exact*, but future
 * changes might increase the width of "on" to a finite but small range
 * ("fuzziness").
 */
enum class SignedSense
{
    inside  = -1,
    on      = 0,
    outside = 1
};

//---------------------------------------------------------------------------//
/*!
 * Enumeration for mapping surface classes to integers.
 */
enum class SurfaceType : unsigned char
{
    px = 0, //!< Plane normal to X axis
    py,     //!< Plane normal to Y axis
    pz,     //!< Plane normal to Z axis
    so,     //!< Sphere at the origin
    cxo,    //!< Cylinder along X axis
    cyo,    //!< Cylinder along Y axis
    czo,    //!< Cylinder along Z axis
    p,      //!< General plane
    s,      //!< Sphere
    cx,     //!< Cylinder parallel to X axis
    cy,     //!< Cylinder parallel to Y axis
    cz,     //!< Cylinder parallel to Z axis
    kx,     //!< Cone parallel to X axis
    ky,     //!< Cone parallel to Y axis
    kz,     //!< Cone parallel to Z axis
    sq,     //!< Simple quadric
    gq,     //!< General quadric
    size_
};

//---------------------------------------------------------------------------//
/*!
 * Type-deleted struct for viewing a surface's data.
 */
struct GenericSurfaceRef
{
    span<const real_type> data;
    SurfaceType           type;

    template<class S>
    static GenericSurfaceRef from_surface(const S& surf)
    {
        static_assert(std::is_standard_layout<S>::value,
                      "Surface needs standard layout to get a generic view");
        return {{surf.data(), S::size()}, S::surface_type()};
    }
};

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Evaluate the sense based on the LHS of the quadric equation.
 *
 * This is an optimized jump-free version of:
 * \code
    return quadric == 0 ? SignedSense::on
        : quadric < 0 ? SignedSense::inside
        : SignedSense::outside;
 * \endcode
 * as
 * \code
    int gz = !(quadric <= 0) ? 1 : 0;
    int lz = quadric < 0 ? 1 : 0;
    return static_cast<SignedSense>(gz - lz);
 * \endcode
 * and compressed into a single line.
 *
 * NaN values are treated as "outside".
 */
CELER_FORCEINLINE_FUNCTION SignedSense real_to_sense(real_type quadric)
{
    return static_cast<SignedSense>(!(quadric <= 0) - (quadric < 0));
}

//---------------------------------------------------------------------------//
/*!
 * Convert a signed sense to a Sense enum.
 */
CELER_FORCEINLINE_FUNCTION Sense to_sense(SignedSense s)
{
    return Sense(static_cast<int>(s) >= 0);
}

//---------------------------------------------------------------------------//
// Convert sense to C string version of enum
const char* to_cstring(SignedSense);

//---------------------------------------------------------------------------//
// Convert surface type to C string version of enum
const char* to_cstring(SurfaceType);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

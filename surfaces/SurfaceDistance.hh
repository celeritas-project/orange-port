//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SurfaceDistance.hh
 * \brief SurfaceDistance class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include <cmath>
#include "base/VectorFunctions.hh"
#include "../Fuzziness.hh"
#include "Definitions.hh"
#include "PlaneAligned.hh"
#include "Plane.hh"
#include "CylCentered.hh"
#include "CenteredSphere.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * \struct SurfaceDistance
 * Calculate distances between surfaces.
 *
 * The result should always be positive. For unknown distances we return
 * "no_intersection".
 */
template<class Surface>
struct SurfaceDistance
{
    using Surface_t = Surface;

    //! Return the distance between two surfaces if it's the same at all points
    static real_type global(const Surface_t&, const Surface_t&)
    {
        // Default: we don't know
        return no_intersection();
    }
};

//---------------------------------------------------------------------------//
//! Distance for T-aligned plane
template<Axis T>
struct SurfaceDistance<PlaneAligned<T>>
{
    using Surface_t = PlaneAligned<T>;

    //! Distance between the two planes
    static real_type global(const Surface_t& lhs, const Surface_t& rhs)
    {
        return std::fabs(rhs.position() - lhs.position());
    }
};

//---------------------------------------------------------------------------//
//! Distance for centered spheres
template<>
struct SurfaceDistance<CenteredSphere>
{
    using Surface_t = CenteredSphere;

    //! Distance between the two planes
    static real_type global(const Surface_t& lhs, const Surface_t& rhs)
    {
        return std::fabs(std::sqrt(rhs.radius_sq())
                         - std::sqrt(lhs.radius_sq()));
    }
};

//---------------------------------------------------------------------------//
//! Distance for general plane
template<>
struct SurfaceDistance<Plane>
{
    using Surface_t = Plane;

    //! Distance between two planes
    static real_type global(const Surface_t& lhs, const Surface_t& rhs)
    {
        real_type dot_product = dot_product(lhs.normal(), rhs.normal());
        if (std::fabs(dot_product - 1)
            > fuzziness().surface_simplification_abs())
        {
            // SurfaceContainer aren't parallel, so they intersect
            return 0;
        }

        // Multiply by the dot product in case surfaces point in opposite
        // directions
        return std::fabs(rhs.displacement() - dot_product * lhs.displacement());
    }
};

//---------------------------------------------------------------------------//
//! Distance for on-T centered cylinder
template<Axis T>
struct SurfaceDistance<CylCentered<T>>
{
    using Surface_t = CylCentered<T>;

    //! Distance between the two planes
    static real_type global(const Surface_t& lhs, const Surface_t& rhs)
    {
        return std::fabs(std::sqrt(rhs.radius_sq())
                         - std::sqrt(lhs.radius_sq()));
    }
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

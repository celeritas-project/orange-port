//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file Definitions.hh
 * \brief Definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include "base/Definitions.hh"
#include "base/OpaqueId.hh"
#include "base/Array.hh"
#include "orange/Definitions.hh"

namespace geometria
{
class Transform;
class BoundingBox;
} // namespace geometria

namespace celeritas
{
//---------------------------------------------------------------------------//
// TYPE ALIASES
//---------------------------------------------------------------------------//

//@{
//! Geometry types
using Transform      = geometria::Transform;
using BoundingBox    = geometria::BoundingBox;
using SpanReal3      = geometria::SpanReal3;
using SpanConstReal3 = geometria::SpanConstReal3;
using Real3          = geometria::Real3;
using SpaceMatrix    = geometria::SpaceMatrix;
using Axis           = def::XYZ;
using size_type      = def::size_type;
//@}

//---------------------------------------------------------------------------//
// ID
//---------------------------------------------------------------------------//

namespace id
{
struct Face;
struct Surface;
struct Volume;
struct Universe;
struct Shape;
} // namespace id

//@{
//! Unit-local ID types
using SurfaceId  = OpaqueId<id::Surface, std::uint_least32_t>;
using VolumeId   = OpaqueId<id::Volume, std::uint_least32_t>;
using UniverseId = OpaqueId<id::Universe, std::uint_least32_t>;
using ShapeId    = OpaqueId<id::Shape>;
//@}

//---------------------------------------------------------------------------//
// CONSTANTS
//---------------------------------------------------------------------------//
// std::isinf on old intel compilers (known broken on intel 14, possibly
// others) has a bug that always returns false for constexpr values of inf.
#ifdef __INTEL_COMPILER
#    if __INTEL_COMPILER < 1800
#        define ORANGE_INTEL_CONSTEXPR_BUG
#    endif
#endif

#ifndef ORANGE_INTEL_CONSTEXPR_BUG
// Typical case: constexpr is constexpr
#    define ORANGE_CONSTEXPR constexpr
#else
// Old intel versions: constexpr breaks isinf
#    define ORANGE_CONSTEXPR
#endif

//! Value a real_type shall take to imply that no intersection occurs
CELER_FORCEINLINE_FUNCTION ORANGE_CONSTEXPR real_type no_intersection()
{
    return HUGE_VAL;
}

#undef ORANGE_CONSTEXPR

//---------------------------------------------------------------------------//
// ENUMS
//---------------------------------------------------------------------------//
/*!
 * Whether a position is logically "inside" or "outside" a surface.
 *
 * For a plane, "pos" (outside/true) is equivalent to
 * \f[
   \vec x \cdot \vec n >= 0
 * \f]
 * and "inside" (neg) is to the left of the plane's normal. Likewise, for a
 * sphere, "inside" is where the dot product of the position and outward normal
 * is negative. These are *opposite* the signs from KENO, where `-` has the
 * sense of negating a closed region of space (i.e. the interior of a sphere).
 *
 * \deprecated We should consider different classes for "signs with respect to
 * the plane" versus "negating a CSG node".
 */
enum Sense : bool
{
    neg = 0, //!< Quadric equation evaluates to negative
    pos = 1,
    // Aliases for clarity:
    inside  = neg, //!< '+' sense for KENO region definition vector
    outside = pos, //!< '-' sense for KENO RDV
};

//---------------------------------------------------------------------------//
/*!
 * Whether a position is on or is_crossing a surface.
 *
 * \deprecated These are used only for surfaces and should probably be unified
 * with SignedSense.
 */
enum class SurfaceState : bool
{
    off = 0,
    on  = 1
};

//---------------------------------------------------------------------------//
using zorder_int = std::uint_least16_t;

/*!
 * Priority for "masking" one region with another in KENO cells
 *
 * Regions with higher z order "mask" lower z-order cells.
 * Theoretically you could allow multiple levels of holes, since not having
 * overlapping holes is "a **significant flaw**" ;)
 */
enum class ZOrder : zorder_int
{
    invalid = 0,                        //!< Invalid region
    media,                              //!< Material-filled region or array
    hole,                               //!< Another universe masking this one
    implicit_exterior = zorder_int(-2), //!< Exterior in lower universe
    exterior          = zorder_int(-1), //!< The global problem boundary
};

//---------------------------------------------------------------------------//
// ENUM UTILITIES
//---------------------------------------------------------------------------//
//! Get the C string constant for an axis.
inline const char* to_cstring(Axis ax)
{
    return def::xyz_name(ax);
}

//---------------------------------------------------------------------------//
//! Get a printable character corresponding to a sense.
inline char to_char(Sense s)
{
    return s ? '+' : '-';
}

//---------------------------------------------------------------------------//
//! Convert a boolean value to a Sense enum.
CELER_FORCEINLINE_FUNCTION Sense to_sense(bool s)
{
    return static_cast<Sense>(s);
}

//---------------------------------------------------------------------------//
//! Change the sense across a surface.
CELER_FORCEINLINE_FUNCTION Sense flip_sense(Sense orig)
{
    return static_cast<Sense>(!static_cast<bool>(orig));
}

//---------------------------------------------------------------------------//
//! Whether a zorder corresponds to an exterior region.
CELER_FORCEINLINE_FUNCTION bool is_exterior_zorder(zorder_int z)
{
    return z >= static_cast<zorder_int>(ZOrder::implicit_exterior);
}

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// HASH SPECIALIZATIONS
//---------------------------------------------------------------------------//
namespace std
{
template<>
struct hash<celeritas::Sense>
{
    using argument_type = celeritas::Sense;
    using result_type   = size_type;
    result_type operator()(argument_type sense) const noexcept
    {
        return std::hash<bool>()(static_cast<bool>(sense));
    }
};
} // namespace std

//---------------------------------------------------------------------------//

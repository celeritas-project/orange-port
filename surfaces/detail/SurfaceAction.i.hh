//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/detail/SurfaceAction.i.hh
 * \brief SurfaceAction inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "../PlaneAligned.hh"
#include "../CenteredSphere.hh"
#include "../CylCentered.hh"
#include "../Plane.hh"
#include "../Sphere.hh"
#include "../CylAligned.hh"
#include "../ConeAligned.hh"
#include "../SimpleQuadric.hh"
#include "../GeneralQuadric.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * \struct SurfaceTypeTraits
 * Map surface enumeration to surface type.
 */
template<SurfaceType S>
struct SurfaceTypeTraits;

#define ORANGE_SURFACE_TRAITS(ENUM_VALUE, CLS)        \
    template<>                                        \
    struct SurfaceTypeTraits<SurfaceType::ENUM_VALUE> \
    {                                                 \
        using type = CLS;                             \
    }

// clang-format off
ORANGE_SURFACE_TRAITS(px,  PlaneX);
ORANGE_SURFACE_TRAITS(py,  PlaneY);
ORANGE_SURFACE_TRAITS(pz,  PlaneZ);
ORANGE_SURFACE_TRAITS(so,  CenteredSphere);
ORANGE_SURFACE_TRAITS(cxo, CCylX);
ORANGE_SURFACE_TRAITS(cyo, CCylY);
ORANGE_SURFACE_TRAITS(czo, CCylZ);
ORANGE_SURFACE_TRAITS(p,   Plane);
ORANGE_SURFACE_TRAITS(s,   Sphere);
ORANGE_SURFACE_TRAITS(cx,  CylX);
ORANGE_SURFACE_TRAITS(cy,  CylY);
ORANGE_SURFACE_TRAITS(cz,  CylZ);
ORANGE_SURFACE_TRAITS(kx,  ConeX);
ORANGE_SURFACE_TRAITS(ky,  ConeY);
ORANGE_SURFACE_TRAITS(kz,  ConeZ);
ORANGE_SURFACE_TRAITS(sq,  SimpleQuadric);
ORANGE_SURFACE_TRAITS(gq,  GeneralQuadric);
// clang-format on

#undef ORANGE_SURFACE_TRAITS

//---------------------------------------------------------------------------//
// SURFACE ACTION
//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
template<class F>
SurfaceAction<F>::SurfaceAction(const SurfaceContainer& surfaces, F action)
    : surfaces_(surfaces), action_(action)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply to the surface specified by a surface ID
 */
template<class F>
auto SurfaceAction<F>::operator()(SurfaceId id) -> result_type
{
    CELER_EXPECT(id < surfaces_.size());
#define SURF_APPLY_IMPL(TYPE) \
    case (SurfaceType::TYPE): \
        return this->apply_impl<SurfaceType::TYPE>(id)

    switch (surfaces_.get_type(id))
    {
        SURF_APPLY_IMPL(px);
        SURF_APPLY_IMPL(py);
        SURF_APPLY_IMPL(pz);
        SURF_APPLY_IMPL(so);
        SURF_APPLY_IMPL(cxo);
        SURF_APPLY_IMPL(cyo);
        SURF_APPLY_IMPL(czo);
        SURF_APPLY_IMPL(p);
        SURF_APPLY_IMPL(s);
        SURF_APPLY_IMPL(cx);
        SURF_APPLY_IMPL(cy);
        SURF_APPLY_IMPL(cz);
        SURF_APPLY_IMPL(kx);
        SURF_APPLY_IMPL(ky);
        SURF_APPLY_IMPL(kz);
        SURF_APPLY_IMPL(sq);
        SURF_APPLY_IMPL(gq);
        case SurfaceType::size_:
            CELER_ASSERT_UNREACHABLE();
    }
#undef SURF_APPLY_IMPL
    CELER_ASSERT_UNREACHABLE();
}

//---------------------------------------------------------------------------//
// PRIVATE INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Apply to the surface specified by a surface ID
 */
template<class F>
template<SurfaceType ST>
auto SurfaceAction<F>::apply_impl(SurfaceId id) -> result_type
{
    using Surface_t = typename SurfaceTypeTraits<ST>::type;
    return action_(surfaces_.get<Surface_t>(id));
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SurfaceAction.hh
 * \brief SurfaceAction class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "SurfaceContainer.hh"
#include "detail/SurfaceAction.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Helper function for creating a SurfaceAction instance.
 */
template<class F>
inline detail::SurfaceAction<F>
make_surface_action(const SurfaceContainer& surfaces, F&& action)
{
    return detail::SurfaceAction<F>(surfaces, std::forward<F>(action));
}

//---------------------------------------------------------------------------//
} // namespace celeritas

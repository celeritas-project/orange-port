//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/detail/SurfaceAction.hh
 * \brief SurfaceAction class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <utility>
#include "../SurfaceContainer.hh"
#include "../GeneralQuadric.hh"
#include "../Definitions.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Helper class for applying an action functor to a generic surface.
 *
 * The templated operator() of the given functor F must be a Surface class. The
 * `result_type` type alias here uses GeneralQuadric in to represent the "most
 * generic" type the functor accepts.
 */
//---------------------------------------------------------------------------//

template<class F>
class SurfaceAction
{
  public:
    //@{
    //! Public type aliases
    using result_type
        = decltype(std::declval<F>()(std::declval<GeneralQuadric>()));
    //@}

  private:
    //// DATA ////

    const SurfaceContainer& surfaces_;
    F                       action_;

  public:
    // Constructor
    inline SurfaceAction(const SurfaceContainer& surfaces, F action);

    // Apply to the surface specified by a surface ID
    inline result_type operator()(SurfaceId id);

    //! Get the resulting action
    const F& action() const { return action_; }

  private:
    template<SurfaceType ST>
    inline result_type apply_impl(SurfaceId id);
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "SurfaceAction.i.hh"
//---------------------------------------------------------------------------//

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/BuildLogic.hh
 * \brief BuildLogic class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/containers/Span.hh"
#include "orange/construct/UnitRegion.hh"
#include "LogicEvaluator.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * \struct BuildLogic
 * Convert a logic vector of surface IDs to face (within-cell) IDs.
 */
//---------------------------------------------------------------------------//

struct BuildLogic
{
    //@{
    //! Types
    using CSG_logic_int = CSGCell::logic_int;
    using RPN_logic_int = LogicEvaluator::logic_int;
    //@}

    //! Local surface IDs in the cell
    span<const SurfaceId> faces;

    // Convert the logic/surface ID to an RPN logic value
    inline RPN_logic_int operator()(CSG_logic_int orig);
};

//---------------------------------------------------------------------------//
/*!
 * \struct UnbuildLogic
 * Convert face IDs to surface IDs.
 */
//---------------------------------------------------------------------------//

struct UnbuildLogic
{
    //@{
    //! Types
    using VecSurfaceId  = CSGCell::VecSurfaceId;
    using CSG_logic_int = CSGCell::logic_int;
    using RPN_logic_int = LogicEvaluator::logic_int;
    //@}

    //! Local surface IDs in the cell
    span<const SurfaceId> faces;

    // Convert the logic/surface ID to an RPN logic value
    inline CSG_logic_int operator()(RPN_logic_int orig);
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "BuildLogic.i.hh"
//---------------------------------------------------------------------------//

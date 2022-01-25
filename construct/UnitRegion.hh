//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnitRegion.hh
 * \brief UnitRegion class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "orange/BoundingBox.hh"
#include "CSGCell.hh"
#include "orange/Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * \struct UnitRegion
 * Construction definition of a region of space.
 *
 * This is built by shapes and passed as a construction argument to Unit
 * trackers.
 */
//---------------------------------------------------------------------------//

struct UnitRegion
{
    CSGCell     interior;                 //!< Surface logic for "inside"
    ZOrder      zorder = ZOrder::invalid; //!< Masking precedence
    BoundingBox bbox   = geometria::infinite_bbox(); //!< Enclosing volume
    bool        has_internal_surfaces = true; //!< Requires "complex" tracking
};

//! Functor for sorting by *decreasing* zorder, required for
//! MaskedUnitTracker.
struct ZorderPriority
{
    bool operator()(const UnitRegion& lhs, const UnitRegion& rhs) const
    {
        return lhs.zorder > rhs.zorder;
    }
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

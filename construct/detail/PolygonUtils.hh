//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/PolygonUtils.hh
 * \brief Polygon Utility function declarations
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include <vector>
#include "orange/Definitions.hh"
#include "base/Array.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//

// Whether edges associated with the points overlap
bool edges_overlap(const Array<real_type, 2>& point_A,
                   const Array<real_type, 2>& point_B,
                   const Array<real_type, 2>& point_C);

// Whether edges associated with the points intersect
bool edges_intersect(const Array<real_type, 2>& point_A,
                     const Array<real_type, 2>& point_B,
                     const Array<real_type, 2>& point_C,
                     const Array<real_type, 2>& point_D);

// Whether the three points constitute a right turn
bool is_right_turn(const Array<real_type, 2>& point_A,
                   const Array<real_type, 2>& point_B,
                   const Array<real_type, 2>& point_C);

//---------------------------------------------------------------------------//
} // end namespace detail
} // end namespace celeritas

//---------------------------------------------------------------------------//

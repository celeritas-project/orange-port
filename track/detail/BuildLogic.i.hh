//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/BuildLogic.i.hh
 * \brief BuildLogic inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <algorithm>
#include "orange/construct/UnitRegion.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Convert the logic/surface ID to an RPN logic value.
 */
auto BuildLogic::operator()(CSG_logic_int orig) -> RPN_logic_int
{
    if (!CSGCell::is_operator_token(orig))
    {
        // Value is a surface; convert to a surface ID by finding the
        // location in our list of surfaces.
        auto iter
            = std::lower_bound(faces.begin(), faces.end(), SurfaceId(orig));
        CELER_ASSERT(iter != faces.end() && *iter == SurfaceId(orig));

        // Set destination to the index in the local linearized array
        return RPN_logic_int(iter - faces.begin());
    }

    // Convert from surface RPN (using SurfaceId::size_type,
    // max of probably 2**63) to logical RPN (using logic_int, max of
    // probably 2**16)
    switch (orig)
    {
        // clang-format off
        case CSGCell::LOGIC_TRUE: return LogicEvaluator::LOGIC_TRUE;
        case CSGCell::LOGIC_OR:   return LogicEvaluator::LOGIC_OR;
        case CSGCell::LOGIC_AND:  return LogicEvaluator::LOGIC_AND;
        case CSGCell::LOGIC_NOT:  return LogicEvaluator::LOGIC_NOT;
        default: CELER_ASSERT_UNREACHABLE();
            // clang-format on
    }
}

//---------------------------------------------------------------------------//
/*!
 * Convert the RPN/face logic value to CSG region logic value.
 */
auto UnbuildLogic::operator()(RPN_logic_int orig) -> CSG_logic_int
{
    if (!LogicEvaluator::is_operator_token(orig))
    {
        // Value is a face; convert to a surface ID
        CELER_ASSERT(orig < faces.size());
        return CSG_logic_int(faces[orig].get());
    }

    switch (orig)
    {
        // clang-format off
        case LogicEvaluator::LOGIC_TRUE: return CSGCell::LOGIC_TRUE;
        case LogicEvaluator::LOGIC_OR:   return CSGCell::LOGIC_OR;
        case LogicEvaluator::LOGIC_AND:  return CSGCell::LOGIC_AND;
        case LogicEvaluator::LOGIC_NOT:  return CSGCell::LOGIC_NOT;
        default: CELER_ASSERT_UNREACHABLE();
            // clang-format on
    }
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

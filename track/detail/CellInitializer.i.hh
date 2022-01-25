//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/CellInitializer.i.hh
 * \brief CellInitializer inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "LogicEvaluator.hh"
#include "Utils.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Constructor clears temporary sense vector.
 */
CellInitializer::CellInitializer(const SurfaceContainer&   surfaces,
                                 const LocalState& state)
    : surfaces_(surfaces)
    , pos_(state.pos)
    , surface_(state.surface)
    , sense_(state.sense)
    , temp_senses_(state.temp_senses)
{
    CELER_EXPECT(temp_senses_);
    temp_senses_->clear();
}

//---------------------------------------------------------------------------//
/*!
 * Test the given cell.
 *
 * \return Whether our stored point is inside the cell.
 */
auto CellInitializer::operator()(const Cell& cell_def) -> FoundFace
{
    // Find face inside the cell
    FaceId face = find_face(cell_def.faces, surface_);
    // Update senses, and try to find a local face if not already found
    calc_senses(cell_def, surfaces_, pos_, temp_senses_, &face);

    if (surface_ && face)
    {
        // The particle is *on* a face in the cell being tested, so use the
        // state's sense.
        (*temp_senses_)[face.unchecked_get()] = sense_;
    }

    LogicEvaluator is_inside(cell_def.logic);
    bool           found = is_inside(*temp_senses_);
    return FoundFace{found, sense_, face};
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

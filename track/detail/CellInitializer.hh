//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/CellInitializer.hh
 * \brief CellInitializer class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "orange/surfaces/SurfaceContainer.hh"
#include "../Definitions.hh"
#include "../SenseContainer.hh"
#include "CellContainer.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * See if a position is 'inside' a cell.
 *
 * This both *calculates* and *evaluates* senses. It's assumed that the
 * position is fixed but different cells and senses are being tested.
 */
class CellInitializer
{
  public:
    //@{
    //! Public type aliases
    using Cell = CellContainer::mapped_type;
    //@}

  public:
    // Constructor
    inline CellInitializer(const SurfaceContainer&   surfaces,
                           const LocalState& state);

    // Test the given cell on the given surface with the given sense
    inline FoundFace operator()(const Cell& cell_def);

  private:
    //// DATA ////

    //! Compressed vector of surface definitions
    const SurfaceContainer& surfaces_;

    //! Local position
    SpanConstReal3 pos_;

    //! Local surface
    SurfaceId surface_;

    //! Local sense if on surface
    Sense sense_;

    //! Temporary senses
    SenseContainer* temp_senses_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CellInitializer.i.hh"
//---------------------------------------------------------------------------//

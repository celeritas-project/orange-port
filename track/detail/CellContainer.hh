//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/CellContainer.hh
 * \brief CellContainer class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <vector>
#include "Nemesis/containers/Span.hh"
#include "orange/Definitions.hh"
#include "LogicEvaluator.hh"

namespace celeritas
{
struct UnitRegion;

namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Stored region definitions.
 *
 * The regions must be provided in decreasing zorder. Intersections and flags
 * are expected to be filled in *after* construction. The "number of surfaces"
 * accessor is determined by the highest surface ID given in the input region
 * list.
 */
class CellContainer
{
  public:
    //@{
    //! Public type aliases
    using size_type = VolumeId::size_type;
    using key_type  = VolumeId;
    using logic_int = LogicEvaluator::logic_int;
    using flag_int  = logic_int;
    //@}

    //! Stored cell representation
    struct mapped_type
    {
        zorder_int zorder{0};
        logic_int  num_intersections{0};
        flag_int   flags{0};
        //! Sorted list of surface IDs in this cell
        span<const SurfaceId> faces{};
        //! RPN region definition for this cell, using local surface index
        span<const logic_int> logic{};
    };

    //! Flags
    enum Flags : flag_int
    {
        INTERNAL_SURFACES = 0x1
    };

  public:
    // Empty constructor
    explicit CellContainer() = default;

    // Construct with regions
    explicit CellContainer(const std::vector<UnitRegion>& regions);

    // Allow moving and move assignment, but prevent copies because of span
    CellContainer(CellContainer&&) = default;
    CellContainer& operator=(CellContainer&&) = default;
    CellContainer(const CellContainer&)       = delete;
    CellContainer& operator=(const CellContainer&) = delete;

    //// ACCESSORS ////

    //! Number of cells
    size_type size() const { return cells_.size(); }

    //! Number of corresponding surfaces based on region definitions
    size_type num_surfaces() const { return num_surfaces_; }

    // Get a cell
    inline const mapped_type& operator[](VolumeId id) const;

    // Set flags for a cell
    inline void set_flags(VolumeId id, flag_int flags);

    // Set number of intersections for a cell
    inline void set_num_intersections(VolumeId id, logic_int val);

    // Reconstruct a CSG cell for diagnostic output
    UnitRegion get_region(VolumeId id) const;

  private:
    //// DATA ////

    //! Stored cells
    std::vector<mapped_type> cells_;

    //! Flattened table of surfaces inside each cell (sorted)
    std::vector<SurfaceId> surfids_;

    //! Flattened table of cell logic definitions
    std::vector<logic_int> logic_;

    //! Number of corresponding surfaces (max surface ID + 1)
    size_type num_surfaces_ = 0;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "CellContainer.i.hh"
//---------------------------------------------------------------------------//

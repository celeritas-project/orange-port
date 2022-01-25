//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file Materials.hh
 * \brief Materials class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include "orange/Definitions.hh"
#include "detail/UniverseIndexer.hh"
#include "ORANGEState.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Interface for getting the material ID from the particle state.
 *
 * Currently to satisfy the Shift geometry interface this is rolled into
 * ORANGEGeometry, but as discussed in the context of Oberon and other
 * applications, this could be but one possible interface (templated or virtual
 * or whatever) for obtaining the matid given the state. For example, we could
 * have a OberonMaterials that looks at the local cell at each level to find
 * the unique pin material ID even though the geometry (structure) of the pins
 * might be replicated.
 *
 * As we advance this class (to e.g. add "virtual" radial and axial
 * segments inside a pin to allow sub-pin depletion) it will need to have
 * distance-to-boundary functionality.
 */
class Materials
{
  public:
    //@{
    //! Public type aliases
    using cell_type      = geometria::cell_type;
    using matid_type     = geometria::matid_type;
    using VecMatid       = std::vector<matid_type>;
    using SPConstIndexer = std::shared_ptr<const detail::UniverseIndexer>;
    //@}

  public:
    // Default constructor
    Materials() = default;

    // Construct with matids corresponding to the global cell IDs.
    Materials(VecMatid matids, SPConstIndexer indexer);

    //// ACCESSORS ////

    //! Number of cells
    cell_type num_volumes() const { return matids_.size(); }

    // Get matid from a cell ID (needed for current Shift geometry interface)
    inline matid_type matid(cell_type volid) const;

    //// TEMPLATE INTERFACE ////

    // Get the material ID from the current state.
    inline matid_type matid(const ORANGEState& state) const;

  private:
    //// DATA ////

    VecMatid       matids_;
    SPConstIndexer indexer_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Materials.i.hh"
//---------------------------------------------------------------------------//

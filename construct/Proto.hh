//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/Proto.hh
 * \brief Proto class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <utility>
#include <vector>
#include "orange/Definitions.hh"
#include "orange/Transform.hh"
#include "orange/Definitions.hh"
#include "orange/query/UniverseMetadata.hh"
#include "orange/track/Tracker.hh"

namespace celeritas
{
class PlacedShape;

//---------------------------------------------------------------------------//
/*!
 * Base class for ORANGE universe 'builders'.
 *
 * \note Currently the cell definitions use simple shape/sense pairs, but it
 * could be changed to be a CSG tree with nodes that are shape SPs or
 * operators. (The shapes/detail/CSGNode class could be used for this data
 * structure.)
 */
class Proto
{
  public:
    //@{
    //! Public type aliases
    using matid_type     = geometria::matid_type;
    using SPConstProto   = std::shared_ptr<const Proto>;
    using SPConstShape   = std::shared_ptr<const PlacedShape>;
    using SPMetadata     = std::shared_ptr<UniverseMetadata>;
    using RegionVec      = std::vector<std::pair<Sense, SPConstShape>>;
    using VecProto       = std::vector<SPConstProto>;
    using UPConstTracker = std::unique_ptr<const Tracker>;
    //@}

    struct Daughter
    {
        SPConstProto proto;
        Transform    transform; //!< Parent-to-daughter
    };

    struct BuildArgs
    {
        bool implicit_boundary = true;
    };

    struct BuildResult
    {
        UPConstTracker                               tracker;
        SPMetadata                                   md;
        std::vector<std::pair<VolumeId, matid_type>> matids;
        std::vector<std::pair<VolumeId, Daughter>>   daughters;
    };

  public:
    // Virtual destructor
    virtual ~Proto() = 0;

    //! All metadata
    virtual const ObjectMetadata& metadata() const = 0;

    //! Get the boundary definition, used for defining a hole in a higher level
    virtual const RegionVec& interior() const = 0;

    //! Construct a tracker and metadata from this proto
    virtual BuildResult build(BuildArgs args) const = 0;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

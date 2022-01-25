//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnitBuilder.hh
 * \brief UnitBuilder class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "base/FastHashSet.hh"
#include "orange/query/ObjectMetadata.hh"
#include "orange/query/UnitMetadata.hh"
#include "orange/surfaces/SurfaceContainer.hh"
#include "CSGTree.hh"
#include "UnitRegion.hh"
#include "detail/HashPair.hh"
#include "detail/SurfaceInserter.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct the geometry for a unit from shapes.
 *
 * The added regions must be given in decreasing Z order. At least two cells
 * must be present: an exterior (cell zero) and an interior.
 *
 * This is primarily an implementation detail of UnitProto::build(), but it's
 * also used to construct tests in the Trackers.
 */
class UnitBuilder
{
  public:
    //@{
    //! Public type aliases
    using SPConstShape  = std::shared_ptr<const PlacedShape>;
    using Halfspace     = std::pair<Sense, SPConstShape>;
    using VecSenseShape = std::vector<Halfspace>;
    using VecUnitRegion = std::vector<UnitRegion>;
    //@}

    struct result_type
    {
        SurfaceContainer              surfaces;
        std::vector<UnitRegion>       regions;
        std::shared_ptr<UnitMetadata> md;
    };

  public:
    // Constructor
    UnitBuilder();

    // Reserve space for regions
    void reserve(size_type num_regions);

    // Add an exterior cell using its implicit complement
    void
    exterior(const VecSenseShape& interior, ZOrder zorder, ObjectMetadata md);

    // Build a region in the unit
    VolumeId region(const VecSenseShape& interior,
                    ZOrder               zorder,
                    ObjectMetadata       md,
                    real_type            volume = 0.0);

    // Whether all regions have the same Z order
    bool is_simple() const;

    // Build and invalidate this builder
    result_type operator()(ObjectMetadata unit_md);

  private:
    using VecMetadata  = UnitMetadata::VecMetadata;
    using VecDbl       = UnitMetadata::VecDbl;
    using ShapeFace    = UnitMetadata::ShapeFace;
    using SFHash       = detail::HashPair<ShapeFace>;
    using SetShapeFace = FastHashSet<ShapeFace, SFHash>;
    using NodeId       = CSGTree::NodeId;

    //// DATA ////

    // Cell properties (all same size)
    std::vector<UnitRegion> regions_;
    std::vector<NodeId>     cell_nodes_;
    VecMetadata             cell_md_;
    VecDbl                  volumes_;
    CSGTree                 tree_;

    // Surface properties
    SurfaceContainer          surfaces_;
    std::vector<SetShapeFace> surface_md_;
    detail::SurfaceInserter   insert_surface_;

    // Exterior definition
    VecSenseShape implicit_interior_;
    BoundingBox   interior_bbox_;

    //// IMPLEMENTATION FUNCTIONS ////

    BoundingBox calc_bbox(const VecSenseShape& interior);
    VolumeId    add_volume(Sense                sense,
                           const VecSenseShape& interior,
                           ZOrder               zorder,
                           ObjectMetadata       md,
                           real_type            volume);

    void add_surface(SurfaceId surface, SPConstShape shape, std::string ext);
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

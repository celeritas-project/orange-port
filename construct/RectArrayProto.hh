//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RectArrayProto.hh
 * \brief RectArrayProto class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "ArrayProto.hh"

#include <vector>
#include "base/RegularGrid.hh"
#include "orange/query/ObjectMetadata.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Define a rectangular array of objects.
 *
 * The array's "interior" will always have a lower-left corner at (0,0,0). It
 * is the responsibility of the enclosing unit to use the \c calc_placement
 * methods to transform the array to the parent's coordinate system.
 *
 * The units in the array must be unrotated cuboids.
 */
class RectArrayProto final : public ArrayProto
{
  public:
    //! Construction arguments
    struct Params
    {
        ProtoVec3      units; //!< {i,j,k} array entries from lower left
        ObjectMetadata md;
    };

  public:
    // Constructor
    explicit RectArrayProto(Params params);

    // Default destructor
    ~RectArrayProto();

    //! Get the proto metadata
    const ObjectMetadata& metadata() const final { return md_; }

    // Construct a boundary definition for defining a hole
    const RegionVec& interior() const final;

    // Daughter-to-parent transformation for the array lower-left corner
    Transform calc_placement(SpanConstReal3 pos) const final;

    // Transform to place the center of unit `index` in parent `pos`
    Transform calc_placement(DimVector index, SpanConstReal3 pos) const final;

    // Construct a universe from this proto
    BuildResult build(BuildArgs args) const final;

  private:
    using Grid_t = RegularGrid;

    // Get the 'centrod' of an array daughter
    Real3 cuboid_proto_centroid(DimVector ijk) const;

    //// DATA ////

    ObjectMetadata            md_;
    Grid_t                    grid_;
    RegionVec                 boundary_;
    std::vector<SPConstProto> units_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

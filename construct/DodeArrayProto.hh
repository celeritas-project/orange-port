//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/DodeArrayProto.hh
 * \brief DodeArrayProto class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "ArrayProto.hh"

#include <vector>
#include "orange/BoundingBox.hh"
#include "orange/query/ObjectMetadata.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Define a dodecahedral array of units.
 */
class DodeArrayProto final : public ArrayProto
{
  public:
    struct Params
    {
        ProtoVec3      units; //!< {i,j,k} array entries from lower left
        ObjectMetadata md;
    };

  public:
    // Constructor
    explicit DodeArrayProto(Params params);

    // Default destructor
    ~DodeArrayProto();

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
    //// DATA ////

    ObjectMetadata md_;
    ProtoVec3      units_;
    real_type      apothem_;

    //// IMPLEMENTATION METHODS ////

    // Calc the origin of the given array coordinate
    Real3 calc_origin(const DimVector& coords) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

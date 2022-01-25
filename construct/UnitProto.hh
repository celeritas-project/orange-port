//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnitProto.hh
 * \brief UnitProto class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Proto.hh"

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "orange/Definitions.hh"
#include "orange/Transform.hh"
#include "orange/query/ObjectMetadata.hh"

namespace celeritas
{
class ArrayProto;
//---------------------------------------------------------------------------//
/*!
 * Construct a KENO-like unit ("general universe").
 *
 * The input types allow flexible construction:
 * - Media (cells) with optional matids and cell volumes
 * - Arrays truncated by a region definition vector (KENO-VI)
 * - Rectangular arrays with automatic boundaries (KENO-V):
 *   use Hole with \code zorder = ZOrder::media \endcode
 * - Well-connected holes (GG XML and VERA):
 *   use Hole with \code zorder = ZOrder::media \endcode and \c make_hole_shape
 *   to connect with other media/interior
 *
 * Setting the \c otf_error_checking option will *always* create a "masked"
 * tracker which is slower but can catch erroneous overlapping regions.
 */
class UnitProto final : public Proto
{
  public:
    //@{
    //! Public type aliases
    using matid_type   = geometria::matid_type;
    using Transform    = geometria::Transform;
    using SPConstArray = std::shared_ptr<const ArrayProto>;
    //@}

    //// INPUT TYPES ////

    //! A 'media' entry is a cell filled with a material
    struct Media
    {
        RegionVec      interior;
        matid_type     matid  = geometria::invalid_matid();
        real_type      volume = 0;
        ObjectMetadata md;

        // True if fully defined
        explicit inline operator bool() const;
    };

    //! An array is bounded by an interior vector.
    struct Array
    {
        SPConstProto   proto;     //!< Daughter array
        Transform      transform; //!< Daughter-to-parent
        RegionVec      interior;  //!< Array bounds
        ObjectMetadata md;

        // True if fully defined
        explicit inline operator bool() const;
    };

    //! A hole is another unit inserted into this one.
    struct Hole
    {
        SPConstProto   proto;                 //!< Daughter unit
        Transform      transform;             //!< Daughter-to-parent
        ZOrder         zorder = ZOrder::hole; //!< Overlap control
        ObjectMetadata md;

        // True if fully defined
        explicit inline operator bool() const;
    };

    //! The boundary defines the extents of the unit.
    struct Boundary
    {
        RegionVec      interior;
        bool           implicit_boundary = true; //!< Implicitly truncates
        ObjectMetadata md;

        // True if fully defined
        explicit inline operator bool() const;
    };

    struct Params
    {
        std::vector<Media> media;
        std::vector<Array> arrays;
        std::vector<Hole>  holes;
        bool               otf_error_checking = false;
        Boundary           boundary;
        ObjectMetadata     md;
    };

    //// STATIC UTILITIES ////

    // Construct a shape representing the inside of a hole
    static SPConstShape make_hole_shape(const Hole& hole);

  public:
    // Constructor
    explicit UnitProto(Params params);

    // Default destructor
    ~UnitProto();

    //! Get the proto metadata
    const ObjectMetadata& metadata() const final { return data_.md; }

    // Get the boundary definition for defining a hole in a higher level
    const RegionVec& interior() const final;

    // Construct a universe from this proto
    BuildResult build(BuildArgs args) const final;

  private:
    //// DATA ////

    Params data_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

#include "UnitProto.i.hh"

//---------------------------------------------------------------------------//

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/UniverseMetadata.hh
 * \brief UniverseMetadata class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <functional>
#include <string>
#include "orange/Definitions.hh"
#include "ObjectMetadata.hh"

namespace celeritas
{
class Tracker;
class Universe;

namespace detail
{
struct DescribeMetadata;
struct DescribeUniverseMetadata;
} // namespace detail
//---------------------------------------------------------------------------//
/*!
 * Abstract interface for managing user-facing metadata.
 *
 * The term "local" here means "inside this universe". All IDs are local.
 */
class UniverseMetadata
{
  public:
    //@{
    //! Public type aliases
    using volume_int  = VolumeId::size_type;
    using surface_int = SurfaceId::size_type;
    //@}

  public:
    // Virtual destructor for polymorphism
    virtual ~UniverseMetadata() = 0;

    //! Get metadata about this universe
    virtual const ObjectMetadata& metadata() const = 0;

    //! Local bounding box
    virtual BoundingBox bbox() const = 0;

    //! Number of local surfaces
    virtual surface_int num_surfaces() const = 0;

    //! Number of local cells
    virtual volume_int num_volumes() const = 0;

    //! Convert a local surface ID into a user-facing surface label
    virtual std::string id_to_label(SurfaceId) const = 0;

    //! Convert a local cell ID into a user-facing cell label
    virtual std::string id_to_label(VolumeId) const = 0;

    //! Get the volume of a cell
    virtual real_type volume(VolumeId) const = 0;

    //! Set the the volume of a cell
    virtual void set_volume(VolumeId id, real_type volume) = 0;

    //! Describe using data from the corresponding tracker
    virtual void describe(std::ostream& os, const Tracker& tracker) const = 0;
};

//---------------------------------------------------------------------------//
// INLINE FREE FUNCTIONS
//---------------------------------------------------------------------------//
// Return an object that calls 'describe' which can be streamed
inline detail::DescribeMetadata
to_stream(const UniverseMetadata&, const Tracker&);

//---------------------------------------------------------------------------//
// Return a streamable object for universe metadata plus daughter info
inline detail::DescribeUniverseMetadata
to_stream(const UniverseMetadata&,
          const Universe&,
          std::function<const ObjectMetadata&(UniverseId)>,
          std::function<void(std::ostream&, VolumeId)> = nullptr);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "UniverseMetadata.i.hh"
//---------------------------------------------------------------------------//

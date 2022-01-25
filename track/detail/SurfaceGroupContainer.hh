//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/SurfaceGroupContainer.hh
 * \brief SurfaceGroupContainer class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <cstddef>
#include <vector>
#include "Nemesis/containers/Span.hh"
#include "../Definitions.hh"

namespace celeritas
{
namespace detail
{
class UnitTrackingStorage;
//---------------------------------------------------------------------------//
/*!
 * Storage for groups of surface IDs sorted by Z order.
 *
 * When moving across a surface in the masked tracker, it's necessary to
 * test intersections for all surfaces with the same Z order to ensure there
 * are no overlaps. Each item in this container points to a set of surface IDs
 * for a given Z order.
 */
class SurfaceGroupContainer
{
  public:
    //@{
    //! Public type aliases
    using face_int = FaceId::size_type;
    //@}

    //! Group of surfaces
    struct value_type
    {
        zorder_int            zorder{0};
        face_int              num_intersections{0};
        span<const SurfaceId> faces{};
    };

    //@{
    //! Vector types
    using Storage_t      = std::vector<value_type>;
    using const_iterator = Storage_t::const_iterator;
    //@}

  public:
    // Empty constructor
    explicit SurfaceGroupContainer() = default;

    // Construct with regions
    explicit SurfaceGroupContainer(const UnitTrackingStorage& unit_storage);

    // Allow moving and move assignment, but prevent copies because of span
    SurfaceGroupContainer(SurfaceGroupContainer&&) = default;
    SurfaceGroupContainer& operator=(SurfaceGroupContainer&&) = default;
    SurfaceGroupContainer(const SurfaceGroupContainer&)       = delete;
    SurfaceGroupContainer& operator=(const SurfaceGroupContainer&) = delete;

    //// ACCESSORS ////

    //! Number of surface groups
    size_type size() const { return groups_.size(); }

    //@{
    //! Iterators
    const_iterator begin() const { return groups_.begin(); }
    const_iterator end() const { return groups_.end(); }
    //@}

  private:
    //// DATA ////
    Storage_t              groups_;
    std::vector<SurfaceId> surface_ids_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//

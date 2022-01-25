//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SurfaceContainer.hh
 * \brief SurfaceContainer class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <vector>
#include "base/OpaqueIdRange.hh"
#include "base/CompressedEnumUint.hh"
#include "Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Efficient collection of type-deleted surfaces.
 */
class SurfaceContainer
{
  public:
    //@{
    //! Public type aliases
    using size_type      = SurfaceId::size_type;
    using SurfaceIdRange = OpaqueIdRange<SurfaceId>;
    //@}

  public:
    // Default constructor with no surfaces
    SurfaceContainer() = default;

    // Reserve space for surfaces
    inline void reserve(size_type count);

    // Append a typed surface to the vector
    template<class S>
    inline void push_back(const S& surf);

    // Append a generic surface to the vector
    inline void push_back(GenericSurfaceRef surf);

    // Resize to a smaller value to erase elements
    inline void resize(size_type count);

    //// ACCESSORS ////

    // Get the "type" (internal enumeration) of a surface
    inline SurfaceType get_type(SurfaceId surface_id) const;

    // Whether the surface has the same type as the given one
    template<class S>
    inline bool is_type(SurfaceId surface_id) const;

    // Get a surface
    template<class S>
    inline S get(SurfaceId surface_id) const;

    // Get a surface view
    inline GenericSurfaceRef get_view(SurfaceId surface_id) const;

    //! Number of surfaces
    size_type size() const { return surfaces_.size(); }

    //! Whether no surfaces are present
    bool empty() const { return surfaces_.empty(); }

    // Range for looping over all valid IDs
    inline SurfaceIdRange all_ids() const;

  private:
    //// TYPES ////

    using PairSurfaceOffset = CompressedEnumUint<SurfaceType, size_type>;

    //// DATA ////

    //! Surface type, and offset index into coeff_
    std::vector<PairSurfaceOffset> surfaces_;

    //! Surface coefficients (all surfaces, flattened)
    std::vector<real_type> coeff_;

    //// FUNCTIONS ////

    inline size_type get_offset(SurfaceId) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "SurfaceContainer.i.hh"
//---------------------------------------------------------------------------//

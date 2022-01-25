//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/SurfaceContainer.i.hh
 * \brief SurfaceContainer inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "base/Macros.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
// CONSTRUCTION
//---------------------------------------------------------------------------//
/*!
 * Reserve space for surfaces
 */
void SurfaceContainer::reserve(size_type count)
{
    surfaces_.reserve(count);
    coeff_.reserve(count);
}

//---------------------------------------------------------------------------//
/*!
 * Add a surface (type-deleted) with the given coefficients
 */
template<class S>
void SurfaceContainer::push_back(const S& surface)
{
    this->push_back(surface.view());
}

//---------------------------------------------------------------------------//
/*!
 * Append a surface to the vector
 */
void SurfaceContainer::push_back(GenericSurfaceRef surf_view)
{
    // Add compressed surface ref: type and data offset
    surfaces_.emplace_back(surf_view.type, coeff_.size());

    // Extend stored data
    coeff_.insert(coeff_.end(), surf_view.data.begin(), surf_view.data.end());
}

//---------------------------------------------------------------------------//
/*!
 * Resize to a smaller value to erase elements.
 *
 * This is intended so that you can save the surface container size, add some
 * surfaces, and then "revert" the addition by restoring to the same size.
 *
 * It is *not* meaningful to increase the container size so that is prohibited.
 */
void SurfaceContainer::resize(size_type count)
{
    CELER_EXPECT(count <= this->size());

    if (count >= this->size())
        return;

    auto coeff_size = surfaces_[count].get_int();
    surfaces_.resize(count);
    coeff_.resize(coeff_size);
}

//---------------------------------------------------------------------------//
// ACCESS
//---------------------------------------------------------------------------//
/*!
 * Get the surface type enumeration for a surface ID
 */
CELER_FORCEINLINE_FUNCTION auto
SurfaceContainer::get_type(SurfaceId surface_id) const -> SurfaceType
{
    CELER_EXPECT(surface_id < surfaces_.size());
    return surfaces_[surface_id.unchecked_get()].get_enum();
}

//---------------------------------------------------------------------------//
/*!
 * Whether the surface has the same type as the given one
 */
template<class S>
bool SurfaceContainer::is_type(SurfaceId surface_id) const
{
    CELER_EXPECT(surface_id < surfaces_.size());
    return (S::surface_type() == surfaces_[surface_id.get()].get_enum());
}

//---------------------------------------------------------------------------//
/*!
 * Get a surface.
 */
template<class S>
CELER_FORCEINLINE_FUNCTION S SurfaceContainer::get(SurfaceId surface_id) const
{
    CELER_EXPECT(surface_id < surfaces_.size());
    return S(coeff_.data() + this->get_offset(surface_id));
}

//---------------------------------------------------------------------------//
/*!
 * Get a view to the data and type of a surface
 */
CELER_FORCEINLINE_FUNCTION auto
SurfaceContainer::get_view(SurfaceId surface_id) const -> GenericSurfaceRef
{
    CELER_EXPECT(surface_id < surfaces_.size());
    auto start = this->get_offset(surface_id);
    auto stop
        = surface_id.unchecked_get() + 1 < this->size()
              ? this->get_offset(SurfaceId{surface_id.unchecked_get() + 1})
              : coeff_.size();
    CELER_ASSERT(start <= stop && stop <= coeff_.size());

    GenericSurfaceRef result;
    result.data = {coeff_.data() + start, coeff_.data() + stop};
    result.type = this->get_type(surface_id);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get a range of all valid IDs
 */
auto SurfaceContainer::all_ids() const -> SurfaceIdRange
{
    return id_range(SurfaceId(this->size()));
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//
/*!
 * \brief Get the index of the start of the data for a surface
 */
CELER_FORCEINLINE_FUNCTION auto
SurfaceContainer::get_offset(SurfaceId surface_id) const -> size_type
{
    return surfaces_[surface_id.unchecked_get()].get_int();
}

//---------------------------------------------------------------------------//
} // namespace celeritas

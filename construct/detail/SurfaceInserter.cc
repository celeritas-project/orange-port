//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/SurfaceInserter.cc
 * \brief SurfaceInserter class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SurfaceInserter.hh"

#include <cmath>
#include "base/Future.hh"
#include "orange/Fuzziness.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * Construct with a reference to empty surfaces
 */
SurfaceInserter::SurfaceInserter(SurfaceContainer* surfaces)
    : surfaces_(surfaces)
{
    CELER_EXPECT(surfaces && surfaces->empty());
    // Calculate the number of base-10 digits of relative tolerance
    const real_type rel_err = fuzziness().surface_elision_rel();
    CELER_ASSERT(rel_err > 0 && rel_err < 1);
    deduplication_digits_ = static_cast<int>(std::round(-std::log10(rel_err)));
    CELER_ENSURE(deduplication_digits_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * Insert a generic surface.
 */
SurfaceId SurfaceInserter::operator()(GenericSurfaceRef surf_view)
{
    // Hash the surface type and surface coefficients
    HashKey hash = this->hash_surface(surf_view);

    SurfaceId result;

    auto hash_iter = surface_hash_.find(hash);
    if (hash_iter != surface_hash_.end())
    {
        result = hash_iter->second;
    }
    else
    {
        result = SurfaceId{surfaces_->size()};
        // Insert the new surface and its coeffients, and its name
        surfaces_->push_back(surf_view);
        // Add the new surface and its coefficients into our hash
        surface_hash_.insert({hash, result});
    }

    return result;
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * Hash a surface
 */
inline auto SurfaceInserter::hash_surface(GenericSurfaceRef surf_view) const
    -> HashKey
{
    SoftHasher hasher(deduplication_digits_);
    hasher << surf_view.type;
    // Since the array size is fixed for a particular surface type, don't
    // bother hashing it: just save the data
    hasher.save_array(surf_view.data.data(), surf_view.data.size());
    return hasher.hash();
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

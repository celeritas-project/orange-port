//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/SurfaceInserter.hh
 * \brief SurfaceInserter class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <vector>

#include "Nemesis/serialize/SoftHasher.hh"
#include "base/FastHashMap.hh"
#include "orange/surfaces/SurfaceContainer.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Insert deduplicated surfaces.
 *
 * This replaces the construction/deduplication functionality of the
 * \c gg::SurfaceContainer class.
 *
 * \note Hidden data usage:
 *   The constructor queries `fuzziness()` for the surface tolerance accuracy.
 *
 * \code
   SurfaceInserter insert_surface(surfaces);
   auto id = insert_surface(PlaneX(123));
   auto id2 = insert_surface(PlaneX(123.0000001)); // equals id
   \endcode
 */
class SurfaceInserter
{
  public:
    // Construct with reference to surfaces to build
    explicit SurfaceInserter(SurfaceContainer* surfaces_);

    // Add a new surface
    template<class S>
    inline SurfaceId operator()(const S& surface);

    // Append a generic surface view to the vector
    SurfaceId operator()(GenericSurfaceRef surf);

    //! Access the surfaces being inserted into
    const SurfaceContainer& surfaces() const { return *surfaces_; }

    //! Approx number of base-10 digits used for deduplication
    int deduplication_digits() const { return deduplication_digits_; }

  private:
    //// TYPE ALIASES ////

    using HashKey = SoftHasher::HashKey;

    //// DATA ////

    //! SurfaceContainer being built
    SurfaceContainer* surfaces_;

    //! Deduplication of similar surfaces during construction
    FastHashMap<HashKey, SurfaceId> surface_hash_;

    //! Approximate number of base-10 digits to preserve
    int deduplication_digits_;

    //// FUNCTIONS ////

    // Hash a surface
    HashKey hash_surface(GenericSurfaceRef surf) const;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "SurfaceInserter.i.hh"
//---------------------------------------------------------------------------//

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/TrackingError.hh
 * \brief TrackingError class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <stdexcept>
#include <vector>
#include "Nemesis/containers/Span.hh"
#include "base/Macros.hh"

//---------------------------------------------------------------------------//
/*!
 * Raise a TrackingError assertion if the condition is satisfied.
 */
#define ORANGE_TRACKING_ASSERT(COND, EXCEPTION_OBJ) \
    do                                              \
    {                                               \
        if (NEMESIS_UNLIKELY(!(COND)))              \
        {                                           \
            throw EXCEPTION_OBJ;                    \
        }                                           \
    } while (0)

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Base class for exceptions indicating a failure in tracking.
 *
 * Tracking errors do not generally contain the basic geometry state
 * information: it's assumed that ORANGEGeometry will catch these and rethrow
 * them as \c geometria::GeometryError with local state information translated
 * to globally useful names.
 */
class TrackingError : public std::runtime_error
{
  public:
    // Default constructors
    using std::runtime_error::runtime_error;
};

//---------------------------------------------------------------------------//
/*!
 * Raise whenever no valid region could be determined.
 */
class NoCellError : public TrackingError
{
  public:
    //! Construct with default error message
    NoCellError() : TrackingError("No cell is defined") {}
};

//---------------------------------------------------------------------------//
/*!
 * Raise whenever multiple valid regions are detected.
 *
 * This is used in the error detection mode of the \c MaskedUnitTracker .
 */
class OverlappingCellError : public TrackingError
{
  public:
    //@{
    //! Public type aliases
    using SpanConstCellId = span<const VolumeId>;
    //@}

  public:
    //! Construct with local IDs
    explicit OverlappingCellError(SpanConstCellId ids)
        : TrackingError("Multiple overlapping cells found")
        , local_vols_{ids.begin(), ids.end()}
    {
    }

    //! Access the list of overlapping cells
    SpanConstCellId local_vols() const { return make_span(local_vols_); }

  private:
    std::vector<VolumeId> local_vols_;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

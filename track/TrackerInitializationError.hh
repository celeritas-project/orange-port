//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/TrackerInitializationError.hh
 * \brief TrackerInitializationError class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <stdexcept>
#include <vector>
#include "../Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Low-level error raised when initialization failed.
 *
 * This should be converted to a higher-level error based on the current
 * universe ID, the associated metadata, and the information returned by this
 * exception.
 */
class TrackerInitializationError : public std::runtime_error
{
  public:
    //@{
    //! Public type aliases
    std::vector<VolumeId> VecCell;
    //@}

  public:
    // Constructor
    TrackerInitializationError(std::string description,
                               VecCell     nearby_ids,
                               const char* file,
                               int         line);

  private:
    //// DATA ////
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

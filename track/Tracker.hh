//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/Tracker.hh
 * \brief Tracker class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Definitions.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Interface class for tracking through a single 'level' of geometry.
 */
class Tracker
{
  public:
    // Constructor
    virtual ~Tracker() = 0;

    //// TRACKING ////

    // Find the local cell and possibly surface ID.
    virtual Initialization initialize(LocalState state) const = 0;

    // Calculate distance-to-intercept for the next surface
    virtual Intersection intersect(LocalState state) const = 0;

    // Calculate normal on the current surface
    virtual Real3 normal(LocalState state) const = 0;

    //// ACCESSORS ////

    //! Number of cells
    virtual size_type num_volumes() const = 0;

    //! Number of surfaces
    virtual size_type num_surfaces() const = 0;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

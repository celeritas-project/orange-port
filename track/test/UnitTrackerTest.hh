//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/test/UnitTrackerTest.hh
 * \brief UnitTrackerTest class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "TrackerTest.hh"

#include "base/Future.hh"
#include "orange/construct/ShapeContainer.hh"
#include "orange/query/ObjectMetadata.hh"
#include "orange/construct/UnitBuilder.hh"

namespace orange_test
{
//---------------------------------------------------------------------------//
/*!
 * Helper class for unit-based tracker construction and testing.
 */
class UnitTrackerTest : public TrackerTest
{
    using Base = TrackerTest;

  public:
    //@{
    //! Public type aliases
    using SPConstProvenance = std::shared_ptr<const Provenance>;
    using UnitBuilder       = celeritas::UnitBuilder;
    using ObjectMetadata    = celeritas::ObjectMetadata;
    using SPConstShape      = std::shared_ptr<const celeritas::PlacedShape>;
    //@}

    template<class Tracker_t>
    void make_tracker(UnitBuilder&& build, ObjectMetadata unit_md)
    {
        CELER_EXPECT(!tracker);

        // Construct tracker components and tracker
        auto built   = build(std::move(unit_md));
        auto tracker = make_unique<Tracker_t>(std::move(built.surfaces),
                                              std::move(built.regions));

        this->set_tracker(std::move(tracker), std::move(built.md));
    }

  protected:
    //// DATA ////

    celeritas::ShapeContainer shapes;
};

//---------------------------------------------------------------------------//
} // namespace orange_test

//---------------------------------------------------------------------------//

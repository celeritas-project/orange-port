//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/ProtoTest.hh
 * \brief ProtoTest class declaration
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include "Nemesis/gtest/Test.hh"

#include <memory>
#include <string>
#include "orange/construct/PlacedShape.hh"
#include "orange/query/ObjectMetadata.hh"

namespace celeritas
{
class Proto;

//---------------------------------------------------------------------------//
/*!
 * Helper functions for prot test construction.
 *
 * Long description or discussion goes here.
 */
class ProtoTest : public ::Test
{
  public:
    //@{
    //! Public type aliases
    using SPConstShape = std::shared_ptr<const PlacedShape>;
    using SPConstProto = std::shared_ptr<const Proto>;
    //@}

  public:
    // Make a single-material unit with a bounding shape
    SPConstProto make_simple_unit(ObjectMetadata      unit_md,
                                  PlacedShape::Params ps_params) const;

    // Build a proto with a single sphere radius 3, as "media"
    SPConstProto
    make_sphere_unit(ObjectMetadata unit_md, real_type radius = 3.0) const;
};

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//

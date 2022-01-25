//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/ProtoTest.cc
 * \brief ProtoTest class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ProtoTest.hh"

#include "../ShapeContainer.hh"
#include "../SphereShape.hh"
#include "../UnitProto.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Make a single-material unit with a bounding shape
 */
auto ProtoTest::make_simple_unit(ObjectMetadata      unit_md,
                                 PlacedShape::Params ps_params) const
    -> SPConstProto
{
    auto shape = std::make_shared<PlacedShape>(std::move(ps_params));
    UnitProto::Params params;
    {
        UnitProto::Media media;
        media.interior = {{inside, shape}};
        media.volume   = shape->shape()->volume();
        media.md       = ORANGE_MD_FROM_SOURCE("interior");
        media.matid    = 1;
        params.media.push_back(std::move(media));
    }
    {
        UnitProto::Boundary boundary;
        boundary.interior          = {{inside, shape}};
        boundary.implicit_boundary = false; // connected to media
        boundary.md                = ORANGE_MD_FROM_SOURCE("boundary");
        params.boundary            = boundary;
    }
    params.md = std::move(unit_md);
    return std::make_shared<UnitProto>(std::move(params));
}

//---------------------------------------------------------------------------//
/*!
 * Build a proto with a single sphere radius 3, as "media"
 */
auto ProtoTest::make_sphere_unit(ObjectMetadata md, real_type radius) const
    -> SPConstProto
{
    ShapeContainer shapes;

    auto sph = shapes.emplace<SphereShape>(
        ORANGE_MD_FROM_SOURCE("sphere"), Transform{}, radius);

    UnitProto::Params params;
    {
        UnitProto::Media media;
        media.interior = {{inside, sph}};
        media.matid    = 1;
        media.volume   = sph->shape()->volume();
        media.md       = ORANGE_MD_FROM_SOURCE("sphere medium");
        params.media.push_back(std::move(media));
    }
    {
        UnitProto::Boundary boundary;
        boundary.interior          = {{inside, sph}};
        boundary.implicit_boundary = false; // connected to media
        boundary.md     = ORANGE_MD_FROM_SOURCE("sphere unit boundary");
        params.boundary = boundary;
    }
    params.md = md;
    return std::make_shared<UnitProto>(std::move(params));
}

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/ShapeBuilder.i.hh
 * \brief ShapeBuilder inline method definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <utility>

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Build a surface with a custom name
 */
template<class S>
void ShapeBuilder::surface(std::string face_name, Sense sense, S&& surf)
{
    return this->construct(face_name, sense, std::forward<S>(surf));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constrain the pre-transformation bounding box size
 */
void ShapeBuilder::clip_bbox_lower(Axis axis, real_type value)
{
    this->local_bbox().clip_lower(axis, value);
}

//---------------------------------------------------------------------------//
/*!
 * Constrain the pre-transformation bounding box size
 */
void ShapeBuilder::clip_bbox_upper(Axis axis, real_type value)
{
    this->local_bbox().clip_upper(axis, value);
}

//---------------------------------------------------------------------------//
// PRIVATE INLINE HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Current cell's "non-transformed" bounding box
 */
auto ShapeBuilder::local_bbox() -> BoundingBox&
{
    CELER_EXPECT(!stack_.empty());
    return stack_.back().bbox;
}

//---------------------------------------------------------------------------//
/*!
 * Current cell's "global" reference frame bounding box
 */
auto ShapeBuilder::global_bbox() -> BoundingBox&
{
    CELER_EXPECT(!stack_.empty());
    return stack_.back().global_bbox;
}

//---------------------------------------------------------------------------//
/*!
 * Current cell's "global" reference frame bounding box
 */
auto ShapeBuilder::transform() const -> const Transform&
{
    CELER_EXPECT(!stack_.empty());
    return stack_.back().global_transform;
}

//---------------------------------------------------------------------------//
/*!
 * Insert a surface, accounting for rotation/reflection.
 */
template<class Surface>
void ShapeBuilder::construct(std::string face_name, Sense sense, Surface surf)
{
    // Clip pre-transformed bbox
    surf.clip(sense, this->local_bbox());

    const auto& t = this->transform();
    if (t.has_rotation())
    {
        return this->construct_transformed(
            std::move(face_name), sense, surf.transformed(t));
    }
    else if (t.has_translation())
    {
        return this->construct_transformed(
            std::move(face_name), sense, surf.translated(t));
    }
    else
    {
        return this->construct_transformed(std::move(face_name), sense, surf);
    }
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

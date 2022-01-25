//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/Definitions.hh
 * \brief Class definitions used by trackers
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC
 */
//---------------------------------------------------------------------------//
#pragma once

#include <type_traits>
#include <vector>
#include "../Definitions.hh"

namespace nemesis
{
template<typename T>
class Face;
}

namespace celeritas
{
// Forward declarations
class SenseContainer;

//---------------------------------------------------------------------------//
// TYPE ALIASES
//---------------------------------------------------------------------------//

//@{
//! State types
using FaceId     = OpaqueId<id::Face, std::uint_least32_t>;
using NextFace    = std::pair<FaceId, real_type>;
using VecNextFace = std::vector<NextFace>;
template<template<class T> class FaceTraits>
using ArrayFace
    = Face<FaceTraits<typename std::make_signed<FaceId::size_type>::type>>;
//@}

//---------------------------------------------------------------------------//
/*!
 * Properties of the local state.
 *
 * All variables (IDs, position, direction) are *local* to the given tracker.
 * Since this is given by \em copy, it is *not* expected to be modified, except
 * for the temporary storage references.
 *
 * NOTE: the maximum size of either of the temp vectors should just be the
 * maximum number of faces used in any one cell.
 */
struct LocalState
{
    const Real3&    pos;
    const Real3&    dir;
    VolumeId        cell;
    SurfaceId       surface;
    Sense           sense          = Sense::neg; // If on a surface
    SenseContainer* temp_senses    = nullptr;
    VecNextFace*    temp_face_dist = nullptr;
};

//---------------------------------------------------------------------------//
/*!
 * Cell ID and surface ID after initialization.
 *
 * Possible configurations for the initialization result ('X' means 'has
 * a valid ID', i.e. evaluates to true):
 *
 *  Cell  | Surface | Description
 * :----: | :-----: | :-------------------------------
 *        |         | Failed to find new cell
 *   X    |         | Initialized
 *   X    |   X     | Crossed surface into new cell
 *        |   X     | Initialized on a surface (reject)
 *
 */
struct Initialization
{
    VolumeId  cell;
    SurfaceId surface;
    Sense     sense = Sense::neg; // Sense if on a surface

    //! Whether initialization succeeded
    explicit operator bool() const { return static_cast<bool>(cell); }
};

//---------------------------------------------------------------------------//
/*!
 * Distance and next-surface information.
 *
 * The sense will be the post-crossing sense.
 */
//---------------------------------------------------------------------------//

struct Intersection
{
    SurfaceId surface;
    Sense     sense    = Sense::neg;
    real_type distance = no_intersection();

    //! Whether a next surface has been found
    explicit operator bool() const { return static_cast<bool>(surface); }

    template<class archiver>
    void serialize(archiver& ar)
    {
        // clang-format off
        ar & surface & sense & distance;
        // clang-format on
    }
};

//---------------------------------------------------------------------------//
/*!
 * Whether an initialization on a face was successful.
 *
 * - \c NOT_FOUND: not found at all, face is not set
 * - \c FOUND: found, but face may be null if not found on the face.
 */
//---------------------------------------------------------------------------//

struct FoundFace
{
    //// DATA ////

    bool   found{false};
    Sense  sense{Sense::neg}; //!< If on a face
    FaceId face{};

    //// METHODS ////

    //! Construct as 'not found'
    FoundFace() = default;

    //! Construct as 'found' but not on a face
    explicit FoundFace(bool found_) : found(found_) {}

    //! Construct as 'found', possibly on a face
    FoundFace(bool found_, Sense sense_, FaceId face_)
        : found(found_), sense(sense_), face(face_)
    {
    }

    explicit operator bool() const { return this->found; }
};

//---------------------------------------------------------------------------//
} // namespace celeritas

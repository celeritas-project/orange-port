//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/detail/SurfaceFunctors.hh
 * \brief Functor definitions used by SurfaceContainer via SurfaceAction
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include "base/Assert.hh"
#include "orange/Definitions.hh"
#include "orange/surfaces/Definitions.hh"
#include "../Definitions.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
//! Calculate the sense of a surface at a given position.
struct CalcSense
{
    SpanConstReal3 pos;

    template<class S>
    SignedSense operator()(S&& surf)
    {
        return surf.calc_sense(this->pos);
    }
};

//---------------------------------------------------------------------------//
//! Calculate the number of intersections
struct NumIntersections
{
    template<class S>
    size_type operator()(S&&)
    {
        return S::num_intersections();
    }
};

//---------------------------------------------------------------------------//
//! Calculate the outward normal at a position.
struct CalcNormal
{
    SpanConstReal3 pos;

    template<class S>
    Real3 operator()(S&& surf)
    {
        return surf.calc_normal(this->pos);
    }
};

//---------------------------------------------------------------------------//
//! Calculate the smallest distance from a point to the surface.
struct CalcSafetyDistance
{
    SpanConstReal3 pos;

    //! Operate on a surface
    template<class S>
    real_type operator()(S&& surf)
    {
        // Calculate outward normal
        Real3 dir = surf.calc_normal(this->pos);

        // If sense is "positive" (on or outside), flip direction to inward so
        // that the vector points toward the surface
        if (to_sense(surf.calc_sense(this->pos)) == Sense::pos)
            dir *= -1;

        real_type distance[S::num_intersections()];
        surf.calc_intersections(this->pos, dir, SurfaceState::off, distance);
        return distance[0];
    }
};

//---------------------------------------------------------------------------//
/*!
 * Fill an array with distances-to-intersection.
 *
 * This assumes that each call is to the next face index, starting with face
 * zero.
 */
class CalcIntersections
{
  public:
    //@{
    //! Public type aliases
    using NextFaceIter = VecNextFace::iterator;
    using face_int    = FaceId::size_type;
    //@}

  public:
    // Construct from the particle point, direction, and faces
    CalcIntersections(const SpanConstReal3&  pos,
                      const SpanConstReal3&  dir,
                      FaceId                 on_face,
                      celeritas::VecNextFace* face_dist)
        : pos_(pos)
        , dir_(dir)
        , on_face_idx_(on_face.unchecked_get())
        , cur_face_idx_(0)
        , cur_face_dist_(face_dist->begin())
        , end_face_dist_(face_dist->end())
    {
    }

    // Construct from the local state
    CalcIntersections(LocalState state, FaceId on_face)
        : pos_(state.pos)
        , dir_(state.dir)
        , on_face_idx_(on_face.unchecked_get())
        , cur_face_idx_(0)
        , cur_face_dist_(state.temp_face_dist->begin())
        , end_face_dist_(state.temp_face_dist->end())
    {
    }

    //! Operate on a surface
    template<class S>
    void operator()(S&& surf)
    {
        auto on_surface = (on_face_idx_ == cur_face_idx_) ? SurfaceState::on
                                                          : SurfaceState::off;

        // Calculate distance to surface along this direction
        real_type temp_distances[S::num_intersections()];
        surf.calc_intersections(pos_, dir_, on_surface, temp_distances);

        // Copy possible intersections and this surface to the output
        for (real_type dist : temp_distances)
        {
            CELER_ASSERT(cur_face_dist_ != end_face_dist_);
            CELER_ASSERT(dist >= 0);
            cur_face_dist_->first  = FaceId{cur_face_idx_};
            cur_face_dist_->second = dist;
            ++cur_face_dist_;
        }
        // Increment to next face
        ++cur_face_idx_;
    }

    NextFaceIter face_dist_iter() const { return cur_face_dist_; }
    face_int    face_idx() const { return cur_face_idx_; }

  private:
    //// DATA ////

    const SpanConstReal3 pos_;
    const SpanConstReal3 dir_;
    const face_int       on_face_idx_;
    face_int             cur_face_idx_;
    NextFaceIter          cur_face_dist_;
    NextFaceIter          end_face_dist_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/HexArrayTracker.cc
 * \brief HexArrayTracker class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "HexArrayTracker.hh"

#include <type_traits>

#include "base/Face.hh"
#include "base/GridLookup.hh"
#include "base/Range.hh"

#include "orange/Fuzziness.hh"
#include "orange/track/detail/SurfaceFunctors.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/construct/CSGCell.hh"
#include "orange/construct/detail/ShapeBuilder.hh"
#include "orange/construct/detail/SurfaceInserter.hh"
#include "orange/construct/PrismShape.hh"
#include "detail/Utils.hh"

using Axis::x;
using Axis::y;
using Axis::z;

using HexFace        = celeritas::HexArrayTracker::Face;
using HexOrientation = celeritas::HexArrayTracker::Orientation;

constexpr int U = static_cast<int>(HexFace::Axes::u);
constexpr int V = static_cast<int>(HexFace::Axes::v);
constexpr int W = static_cast<int>(HexFace::Axes::w);

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
//! Convert a surface state to a cuboid face
inline HexFace to_face(SurfaceId s)
{
    CELER_EXPECT(s < HexFace::num_faces());
    return HexFace(s.unchecked_get());
}

//---------------------------------------------------------------------------//
//! Convert a surface state to a cuboid face
inline SurfaceId from_face(HexFace s)
{
    CELER_EXPECT(s);
    return SurfaceId(s.unchecked_get());
}

//---------------------------------------------------------------------------//
/*!
 * Obtain the surface normal for the given face and orientation
 */
const Real3& hex_normal(HexFace face, HexOrientation orientation)
{
    CELER_EXPECT(face);
    CELER_EXPECT(orientation == HexOrientation::pointy_top
                 || orientation == HexOrientation::flat_top);
    static_assert(static_cast<int>(HexOrientation::pointy_top) == 0,
                  "bad value for hex orientation enum");

    // Hex face normals for both pointy-topped and flat-topped cases
    constexpr real_type half_sqrt_three   = half * constants::sqrt_three;
    static const Real3  hex_normals[2][8] = {
        // Normals for the hex (orientation: pointy topped)
        {
            {-1.0, 0.0, 0.0},               // -U
            {1.0, 0.0, 0.0},                // +U
            {-half, -half_sqrt_three, 0.0}, // -V
            {half, half_sqrt_three, 0.0},   // +V
            {half, -half_sqrt_three, 0.0},  // -W
            {-half, half_sqrt_three, 0.0},  // +W
            {0.0, 0.0, -1.0},               // -Z
            {0.0, 0.0, 1.0},                // +Z
        },
        // Normals for the hex (orientation: flat topped)
        {
            {-half_sqrt_three, -half, 0.0}, // -U
            {half_sqrt_three, half, 0.0},   // +U
            {0.0, -1.0, 0.0},               // -V
            {0.0, 1.0, 0.0},                // +V
            {half_sqrt_three, -half, 0.0},  // -W
            {-half_sqrt_three, half, 0.0},  // +W
            {0.0, 0.0, -1.0},               // -Z
            {0.0, 0.0, 1.0}                 // +Z
        }};

    return hex_normals[static_cast<int>(orientation)][face.get()];
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance to hex faces inside a hexagon.
 */
class HexFaceIntersections
{
  public:
    //@{
    //! Public type aliases
    using NextFaceIter = VecNextFace::iterator;
    using face_int    = FaceId::size_type;
    using SenseVector = Array<Sense, 6>;
    //@}

  public:
    // Construct from the particle point, direction, and faces
    HexFaceIntersections(const SpanConstReal3 pos,
                         const SpanConstReal3 dir,
                         HexFace              on_face,
                         VecNextFace*          face_dist,
                         const SenseVector&   senses,
                         HexOrientation       orientation)
        : pos_(pos)
        , dir_(dir)
        , on_face_idx_(on_face.unchecked_get())
        , cur_face_idx_(0)
        , cur_face_dist_(face_dist->begin())
        , end_face_dist_(face_dist->end())
        , senses_(senses)
        , orient_(orientation)
    {
        CELER_EXPECT(orient_ == HexOrientation::flat_top
                     || orient_ == HexOrientation::pointy_top);
        CELER_EXPECT(senses_.size() == 6);
    }

    //! Operate on a surface
    template<class S>
    void operator()(S surf)
    {
        CELER_EXPECT(S::num_intersections() == 1);

        auto on_surface = (on_face_idx_ == cur_face_idx_) ? SurfaceState::on
                                                          : SurfaceState::off;

        // Calculate distance to surface along this direction
        real_type distance;
        surf.calc_intersections(pos_, dir_, on_surface, &distance);

        // If 'on' or (outside == pos) or (inside == neg)
        SignedSense sense  = surf.calc_sense(pos_);
        bool        inside = (sense == SignedSense::on)
                      || ((senses_[cur_face_idx_] == pos)
                          == (sense == SignedSense::outside));

        // Check for corner case where particle is tracking
        // across shape plane intersect. Test for zero distance
        // on plane that we are outside and we are headed towards
        if (!inside)
        {
            real_type dot_product = dot_product(
                hex_normal(HexFace(cur_face_idx_), orient_), dir_);

            if (distance == 0 && dot_product > 0)
            {
                distance = celeritas::fuzziness().bump_rel();
            }
            else
            {
                distance = no_intersection();
            }
        }

        CELER_ASSERT(cur_face_dist_ != end_face_dist_);
        cur_face_dist_->first  = FaceId{cur_face_idx_};
        cur_face_dist_->second = distance;
        ++cur_face_idx_;
        ++cur_face_dist_;
    }

  private:
    //// DATA ////

    const SpanConstReal3 pos_;
    const SpanConstReal3 dir_;
    const face_int       on_face_idx_;
    face_int             cur_face_idx_;
    NextFaceIter          cur_face_dist_;
    NextFaceIter          end_face_dist_;
    const SenseVector&   senses_;
    HexOrientation       orient_;
};

//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Construct a hex array tracker.
 */
HexArrayTracker::HexArrayTracker(real_type      apothem,
                                 PlaneDimVector uv_dims,
                                 VecDbl         z_edges,
                                 HexOrientation orientation)
    : apothem_(apothem), z_edges_(std::move(z_edges)), orientation_(orientation)
{
    Insist(orientation == HexOrientation::flat_top
               || orientation == HexOrientation::pointy_top,
           "Missing hex array orientation");
    CELER_EXPECT(z_edges_.size() > 1);
    CELER_EXPECT(z_edges_.front() == 0.0);
    CELER_EXPECT(apothem_ > 0.0);

    // Calculate and set surface normals for first three hex faces
    // (set U/V/W vectors)
    for (int ax : {U, V, W})
    {
        auto norm        = hex_normal(HexFace(ax, true), orientation_);
        uvw_normals_[ax] = PlaneVector(norm[U], norm[V]);
    }

    // UVW vectors for 'find' function are normal divided by apothem
    uvw_ = uvw_normals_ / apothem_;

    // Store the distance between center points of the hexes: uv normals,
    // multiplied by 2 * apothem
    uv_span_[U] = uvw_normals_[U] * apothem * 2;
    uv_span_[V] = uvw_normals_[V] * apothem * 2;

    // Calculate and store the origin of the first internal hex (0,0)
    // translate from lower-left to middle of user hex (0,0), then
    // from middle of user hex (0,0) to middle of internal hex (0,0)

    // Make an indexer on the internal hex array (U, V, Z)
    // -> convert from column-major to row-major
    indexer_ = Indexer(DimVector(uv_dims[U], uv_dims[V], z_edges_.size() - 1));

    // Construct a representative lattice cell from infinite hex prisms
    // POINTY_TOP case requires 30 degree rotation
    real_type rotation
        = (orientation_ == HexOrientation::pointy_top ? half : 0.0);
    const real_type infty = std::numeric_limits<real_type>::infinity();
    PrismShape      shape(6, apothem_, rotation, -infty, infty);

    // Capture the hex surfaces for future intersection logic needs
    SurfaceContainer        temp_surfaces;
    detail::SurfaceInserter insert_surface(&temp_surfaces);
    CSGTree                 tree;
    detail::ShapeBuilder    build_shape(insert_surface, tree);

    // Build surfaces for representative cell
    build_shape.push(CSGCell::LOGIC_AND, {});
    shape.build(build_shape);
    build_shape.pop();

    const CSGNode& node = tree.at(build_shape().cell);
    CELER_ASSERT(!node.is_leaf());

    // Remap from HexFace order (UVWZ major, -+ minor) to shape's construction
    // order (counterclockwise)
    static const SurfaceId::size_type face_to_prism[] = {3, 0, 4, 1, 5, 2};
    surfaces_.reserve(temp_surfaces.size());

    const auto& daughters = node.daughters();
    CELER_ASSERT(senses_.size() == 6);
    CELER_ASSERT(daughters.size() == 6);
    CELER_ASSERT(temp_surfaces.size() == 6);
    for (auto i : range<SurfaceId::size_type>(daughters.size()))
    {
        auto prism_idx = face_to_prism[i];
        // inside daughter -> flip sense: pos wrt face is outside
        senses_[i] = flip_sense(daughters[prism_idx].first);
        // Remap surfaces from prism order to face order
        surfaces_.push_back(temp_surfaces.get_view(SurfaceId{prism_idx}));
    }

    CELER_ENSURE(surfaces_.size() == senses_.size());
}

//---------------------------------------------------------------------------//
/*!
 * Find the new local cell.
 *
 * If the particle is trying to move outside the grid, the cell ID will not
 * change. This is to allow for some fuzziness between the array and the
 * daughter cell.
 */
Initialization HexArrayTracker::initialize(LocalState state) const
{
    if (state.surface)
    {
        CELER_EXPECT(state.cell < num_volumes());
        return this->cross_surface(indexer_.index(state.cell.unchecked_get()),
                                   to_face(state.surface));
    }
    else
    {
        return this->initialize_interior(state);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface.
 */
Intersection HexArrayTracker::intersect(LocalState state) const
{
    CELER_EXPECT(state.cell < num_volumes());

    // Translate cell ID to internal hex coordinates
    auto uvz = indexer_.index(state.cell.unchecked_get());

    // Translate position so it's relative to a hex centered on the origin (for
    // our surfaces_ tracking)
    Real3 local_pos = make_vector(state.pos);
    {
        auto uv_center = this->centroid(PlaneDimVector{uvz[U], uvz[V]});
        local_pos[U] -= uv_center[U];
        local_pos[V] -= uv_center[V];
    }

    Face cur_face = state.surface ? to_face(state.surface) : Face{};

    // Allocate space for all possible hex-surface intersections in the current
    // cell (<= the number of hex faces, 6)
    state.temp_face_dist->resize(senses_.size());

    auto calc_intersections
        = make_surface_action(surfaces_,
                              HexFaceIntersections{local_pos,
                                                   state.dir,
                                                   cur_face,
                                                   state.temp_face_dist,
                                                   senses_,
                                                   orientation_});

    // Calculate distance-to-boundary for all hex faces
    for (SurfaceId surface : surfaces_.all_ids())
    {
        calc_intersections(surface);
    }

    auto crossing = std::min_element(state.temp_face_dist->begin(),
                                     state.temp_face_dist->end(),
                                     detail::CloserFace{});
    CELER_ASSERT(crossing != state.temp_face_dist->end());

    // Calculate distance to Z faces
    const real_type z_dir = state.dir[Z];
    if (z_dir != 0.0)
    {
        // Find the target edge space position
        Face target_face(Face::Axes::z, z_dir > 0);

        if (cur_face != target_face)
        {
            // Find the target edge space position
            int target_coord = int(uvz[Z]) + target_face.is_positive();
            CELER_ASSERT(static_cast<volume_int>(target_coord)
                         < z_edges_.size());
            real_type dist = (z_edges_[target_coord] - local_pos[Z]) / z_dir;
            if (dist > 0 && dist < crossing->second)
            {
                crossing->first  = FaceId(target_face.get());
                crossing->second = dist;
            }
        }
    }

    Intersection result;
    if (crossing->second != no_intersection())
    {
        result.surface  = SurfaceId(crossing->first.get());
        result.sense    = Sense::outside;
        result.distance = crossing->second;
    }

    CELER_ENSURE(!result || result.distance > 0);
    CELER_ENSURE(!result.surface
                 || result.surface.get() < this->num_surfaces());

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate normal to the current surface
 */
Real3 HexArrayTracker::normal(LocalState state) const
{
    CELER_EXPECT(state.surface);
    return hex_normal(to_face(state.surface), orientation_);
}

//---------------------------------------------------------------------------//
/*!
 * Cell id for given logical UVZ coordinate.
 */
VolumeId HexArrayTracker::volume_id(const DimVector& c) const
{
    CELER_EXPECT(c.all_lt(this->dims()));
    return VolumeId(indexer_.index(c));
}

//---------------------------------------------------------------------------//
/*!
 * Circumradius (radius of circle around the hex)
 */
real_type HexArrayTracker::circumradius() const
{
    constexpr real_type two_over_sqrt_three = 1.1547005383792517;
    return apothem_ * two_over_sqrt_three;
}

//---------------------------------------------------------------------------//
/*!
 * Find the logical index given a point
 *
 * It doesn't really matter if points on an edge belong to the left or the
 * right cell.. What we want to ensure is that if the point is outside the
 * array bounds (could happen because of floating point error), it gets lumped
 * into the closest actual cell.
 */
HexArrayTracker::DimVector HexArrayTracker::find(SpanConstReal3 pos) const
{
    Real3 point = pos.vector();

    // Dot product of point with basis vectors, truncated to signed int
    real_type a = std::floor(point[X] * uvw_[U][X] + point[Y] * uvw_[U][Y]);
    real_type b = std::floor(point[X] * uvw_[V][X] + point[Y] * uvw_[V][Y]);
    real_type c = std::floor(point[X] * uvw_[W][X] + point[Y] * uvw_[W][Y]);

    DimVector uvz;

    // See associated python notebooks for an explanation of this magic
    uvz[U] = static_cast<volume_int>(std::floor((1 + a - c) / 3));
    uvz[V] = static_cast<volume_int>(std::ceil((b + c) / 3));

    // Allow one hex over due to floating point error
    const DimVector& dims = this->dims();
    for (auto ax : {U, V})
    {
        if (uvz[ax] == static_cast<volume_int>(-1))
        {
            uvz[ax] = 0;
        }
        else if (uvz[ax] >= dims[ax])
        {
            uvz[ax] = dims[ax] - 1;
        }
    }

    // Look up Z coordinate
    uvz[Z] = grid_lookup_unbounded(z_edges_.begin(), z_edges_.end(), point[Z]);

    CELER_ENSURE(uvz.all_lt(dims));
    return uvz;
}

//---------------------------------------------------------------------------//
/*!
 * Return hex centroid \e (x,y,z) given \e (u,v,z) coordinates.
 */
auto HexArrayTracker::centroid(const DimVector& uvz) const -> Real3
{
    CELER_EXPECT(uvz.all_lt(this->dims()));

    auto xy = this->centroid(PlaneDimVector(uvz[U], uvz[V]));
    return Real3(
        xy[X], xy[Y], half * (z_edges_[uvz[Z]] + z_edges_[uvz[Z] + 1]));
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Initialize inside the array, *not* on a surface.
 *
 * We don't check for being on an edge. The result is *always* on the interior
 * of the array.
 */
Initialization
HexArrayTracker::initialize_interior(LocalState state) const
{
    CELER_EXPECT(!state.surface);

    Initialization result;
    result.cell = this->volume_id(this->find(state.pos));

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Initialize in the cell across the current surface.
 */
Initialization HexArrayTracker::cross_surface(DimVector uvz, Face face) const
{
    CELER_EXPECT(face);

    const DimVector saved_uvz = uvz;

    // Calculate the next cell coordinates and determine the entering face.
    // The out-of-bounds checks rely on unsigned integer overflow arithmetic.
    const auto& dims      = this->dims();
    bool        in_bounds = false;
    switch (face.axis_enum())
    {
        case Face::Axes::u:
        case Face::Axes::v:
            // Traveling along U or V faces
            uvz[face.axis()] += face.normal();
            // Check for being out of bounds
            in_bounds = uvz[face.axis()] < dims[face.axis()];
            break;
        case Face::Axes::w:
            // W uvw moves:  +W -> (-U, +V)
            uvz[U] -= face.normal();
            uvz[V] += face.normal();
            // Check for being out of bounds
            in_bounds = uvz[U] < dims[U] && uvz[V] < dims[V];
            break;
        case Face::Axes::z:
            // Z face
            uvz[Z] += face.normal();
            in_bounds = uvz[Z] < dims[Z];
            break;
        default:
            CELER_ASSERT_UNREACHABLE();
    }

    Initialization result;
    if (in_bounds)
    {
        face = face.opposite();
    }
    else
    {
        // Restore cell indices to put us back inside
        uvz = saved_uvz;
    }
    result.cell    = this->volume_id(uvz);
    result.surface = from_face(face);
    result.sense   = Sense::inside;

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Return hex centroid \e (x,y) given \e (u,v) coordinates.
 */
HexArrayTracker::PlaneVector
HexArrayTracker::centroid(const PlaneDimVector& uv) const
{
    CELER_EXPECT(uv[U] < this->dims()[U]);
    CELER_EXPECT(uv[V] < this->dims()[V]);

    PlaneVector center;
    for (int ax : {U, V})
    {
        center[ax] = uv[U] * uv_span_[U][ax] + uv[V] * uv_span_[V][ax];
    }

    return center;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

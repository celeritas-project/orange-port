//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file track/DodeArrayTracker.cc
 * \brief DodeArrayTracker class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "DodeArrayTracker.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <iomanip>

#include "base/Constants.hh"
#include "base/Definitions.hh"
#include "base/Face.hh"
#include "base/Range.hh"
#include "base/VectorFunctions.hh"
#include "orange/Fuzziness.hh"
#include "orange/track/detail/SurfaceFunctors.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/construct/CSGCell.hh"
#include "orange/construct/detail/ShapeBuilder.hh"
#include "orange/construct/detail/SurfaceInserter.hh"
#include "orange/construct/RhombicDodecahedronShape.hh"
#include "detail/Utils.hh"

using DodeFace       = celeritas::DodeArrayTracker::Face;
using IndexConstView = FixedView<const int, 3>;

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
using DodeNormals = Array<Real3, 12>;
using DodeCoord   = Array<int, 3>;

// Helper function for initializing dode_normal
DodeNormals calc_dode_normals()
{
    DodeNormals result;
    for (auto i : range(DodeFace::num_faces()))
    {
        result[i] = DodeFace{i}.calc_cartesian_normal();
    }
    return result;
}

//---------------------------------------------------------------------------//
//! Get the outward normal for the given face index
SpanConstReal3 dode_normal(DodeFace face)
{
    CELER_EXPECT(face);
    static const DodeNormals normals = calc_dode_normals();
    return make_fixed_view(normals[face.get()]);
}

//---------------------------------------------------------------------------//
//! Face index increment vectors
const DodeCoord& coordinate_update(DodeFace face)
{
    CELER_EXPECT(face);

    // Negative faces are multiplied by -1, but for negative abcd
    // (non-axis-aligned faces) the xy coords are one index lower.
    static const DodeCoord coords[] = {
        {-1, 0, 0},   // -x
        {1, 0, 0},    // +x
        {0, -1, 0},   // -y
        {0, 1, 0},    // +y
        {-1, -1, -1}, // -a
        {0, 0, 1},    // +a : +x+y+z
        {0, -1, -1},  // -b
        {-1, 0, 1},   // +b : -x+y+z
        {0, 0, -1},   // -c
        {-1, -1, 1},  // +c : -x-y+z
        {-1, 0, -1},  // -d
        {0, -1, 1},   // +d : +x-y+z
    };

    return coords[face.get()];
}

//---------------------------------------------------------------------------//
//! Convert a surface state to a cuboid face
inline DodeFace to_face(SurfaceId s)
{
    CELER_EXPECT(s < DodeFace::num_faces());
    return DodeFace(s.unchecked_get());
}

//---------------------------------------------------------------------------//
//! Convert a surface state to a cuboid face
inline SurfaceId to_surface(DodeFace s)
{
    CELER_EXPECT(s);
    return SurfaceId(s.unchecked_get());
}

//---------------------------------------------------------------------------//
/*!
 * Calculate surfaces inside a dodecahedron.
 */
class CalcInternalIntersections
{
  public:
    //@{
    //! Public type aliases
    using NextFaceIter = VecNextFace::iterator;
    using face_int    = FaceId::size_type;
    //@}

  public:
    // Construct from the particle point, direction, and faces
    CalcInternalIntersections(const SpanConstReal3                 pos,
                              const SpanConstReal3                 dir,
                              SurfaceId                            on_face,
                              VecNextFace*                          face_dist,
                              const DodeArrayTracker::SenseVector& senses)
        : pos_(pos)
        , dir_(dir)
        , on_face_idx_(on_face.unchecked_get())
        , cur_face_idx_(0)
        , cur_face_dist_(face_dist->begin())
        , end_face_dist_(face_dist->end())
        , senses_(senses)
    {
        CELER_EXPECT(senses_.size() == 12);
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

        SignedSense sense  = surf.calc_sense(pos_);
        bool        inside = (sense == SignedSense::on)
                      || ((senses_[cur_face_idx_] == pos)
                          == (sense == SignedSense::outside));

        // Check for corner case where particle is tracking
        // across shape plane intersect. Test for zero distance
        // on plane that we are outside and we are headed towards
        // TODO not sure if this is right: we might need to update
        // initialization logic so that 'outside' a face becomes 'on' a face
        if (!inside)
        {
            const auto& normal = dode_normal(DodeFace(cur_face_idx_));
            if (distance == 0 && dot_product(normal, dir_) > 0)
            {
                // Bump distance so it is chosen and tracking is consistent
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

    NextFaceIter face_dist_iter() const { return cur_face_dist_; }
    face_int    face_idx() const { return cur_face_idx_; }

  private:
    //// DATA ////

    const SpanConstReal3                 pos_;
    const SpanConstReal3                 dir_;
    const face_int                       on_face_idx_;
    face_int                             cur_face_idx_;
    NextFaceIter                          cur_face_dist_;
    NextFaceIter                          end_face_dist_;
    const DodeArrayTracker::SenseVector& senses_;
};
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Construct a rect array tracker from a cartesiain grid.
 */
DodeArrayTracker::DodeArrayTracker(real_type apothem, const DimVector& dims)
    : cell_indexer_(dims), apothem_(apothem)
{
    // Construct a representative lattice cell
    RhombicDodecahedronShape shape(apothem_);

    // Capture the surfaces for future intersection logic needs
    detail::SurfaceInserter insert_surface(&this->surfaces_);
    CSGTree                 tree;
    detail::ShapeBuilder    build_shape(insert_surface, tree);

    // Build surfaces for representative cell
    build_shape.push(CSGCell::LOGIC_AND, {});
    shape.build(build_shape);
    build_shape.pop();

    const CSGNode& node = tree.at(build_shape().cell);
    CELER_ASSERT(!node.is_leaf());
    const auto& daughters = node.daughters();
    CELER_ASSERT(daughters.size() == senses_.size());

    for (auto i : range(daughters.size()))
    {
        // inside daughter -> flip sense: pos wrt face is outside
        senses_[i] = flip_sense(daughters[i].first);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Find the new local cell.
 *
 * If the particle is trying to move outside the grid, the cell ID will not
 * change. This is to allow for some fuzziness between the array and the
 * daughter cell.
 *
 * TODO: if we initialize outside of the "valid" region, return
 * Initialization{} so that the higher level geometry can bump and try
 * somewhere else.
 */
Initialization DodeArrayTracker::initialize(LocalState state) const
{
    if (state.surface)
    {
        CELER_EXPECT(state.cell < num_volumes());
        return this->cross_surface(this->cell_to_coords(state.cell),
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
Intersection DodeArrayTracker::intersect(LocalState state) const
{
    CELER_EXPECT(state.cell);

    using Axis::x;
    using Axis::y;
    using Axis::z;
    using constants::sqrt_two;

    // Translate point into cell coordinates
    DimVector coords = this->cell_to_coords(state.cell);

    Real3 cell_pos(2 * apothem_ * coords[X],
                   2 * apothem_ * coords[Y],
                   apothem_ * sqrt_two * coords[Z]);

    // Correct X and Y for even Z array indices
    if (coords[Z] % 2 == 1)
    {
        cell_pos[X] += apothem_;
        cell_pos[Y] += apothem_;
    }
    // Transform state.pos into cell's coordinates
    // this allows surface intersects to perform correctly
    cell_pos = make_vector(state.pos) - cell_pos;

    // Allocate space for all possible surface intersections in the current
    // cell (>= the number of surfaces)
    state.temp_face_dist->resize(12);

    auto calc_intersections = make_surface_action(
        surfaces_,
        CalcInternalIntersections{
            cell_pos, state.dir, state.surface, state.temp_face_dist, senses_});

    // Find all surface intersection distances inside this cell
    for (SurfaceId surface : surfaces_.all_ids())
    {
        calc_intersections(surface);
    }
    CELER_ASSERT(calc_intersections.action().face_dist_iter()
                 == state.temp_face_dist->end());

    auto crossing = std::min_element(state.temp_face_dist->begin(),
                                     state.temp_face_dist->end(),
                                     detail::CloserFace{});
    CELER_ASSERT(crossing != state.temp_face_dist->end());

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
Real3 DodeArrayTracker::normal(LocalState state) const
{
    CELER_EXPECT(state.surface);
    return make_vector(dode_normal(to_face(state.surface)));
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Convert cell ID to ijk coords
 */
auto DodeArrayTracker::cell_to_coords(VolumeId cell) const -> DimVector
{
    CELER_EXPECT(cell && cell.get() < num_volumes());
    return cell_indexer_.index(cell.unchecked_get());
}

//---------------------------------------------------------------------------//
/*!
 * Convert IJK to cell ID
 */
auto DodeArrayTracker::coord_to_cell(const DimVector& ijk) const -> VolumeId
{
    CELER_EXPECT(ijk.all_lt(cell_indexer_.dims()));
    return VolumeId(cell_indexer_.index(ijk));
}

//---------------------------------------------------------------------------//
/*!
 * Initialize inside the array, *not* on a surface.
 *
 * We don't check for being on an edge. The result is *always* on the interior
 * of the array.
 */
Initialization
DodeArrayTracker::initialize_interior(LocalState state) const
{
    CELER_EXPECT(!state.surface);

    DimVector coords = this->find(state.pos);

    Initialization result;
    result.cell = this->coord_to_cell(coords);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Initialize in the cell across the current surface.
 *
 * The surface index and face will change if the cell
 * changes (e.g. from POSX to NEGX).
 */
Initialization
DodeArrayTracker::cross_surface(DimVector coords, Face face) const
{
    CELER_EXPECT(face);

    // Calculate the coordinate transform
    DodeCoord delta = coordinate_update(face);
    if (delta[Axis::z] != 0 && coords[Axis::z] % 2)
    {
        // Update for non-axis-aligned Z movement
        delta[Axis::x] += 1;
        delta[Axis::y] += 1;
    }

    // Apply the coordinate transform to the current coordinate, ensuring
    // clamped to array
    const auto& dims    = cell_indexer_.dims();
    bool        clamped = false;
    for (auto ax : {Axis::x, Axis::y, Axis::z})
    {
        auto orig_coord = coords[ax] + delta[ax];
        coords[ax]      = std::min(orig_coord, dims[ax] - 1);
        clamped         = clamped || (coords[ax] != orig_coord);
    }
    Initialization result;
    result.cell = this->coord_to_cell(coords);

    // Flip face only if we changed cells
    if (!clamped)
    {
        face = face.opposite();
    }
    result.surface = to_surface(face);
    result.sense   = Sense::inside;

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Find the logical index given a point
 *
 * The algorithm works as follows:
 *
 * 1.. Determine the X,Y index and closest even Z index
 * 2. Determine if the point is outside the even Z index
 *    a. Calculate the scalar projection of the axes normals and the point
 *    b. If the scalar projection is greater than the apothem, increment index
 *       according to the axes and sign of the projection
 */
DodeArrayTracker::DimVector DodeArrayTracker::find(SpanConstReal3 pos) const
{
    using Axes = DodeFace::Axes;
    using Axis::x;
    using Axis::y;
    using Axis::z;
    using constants::sqrt_two;

    Real3 point = pos.vector();

    const real_type xy_span = 2 * apothem_;
    const real_type height  = sqrt_two * xy_span;
    const real_type hheight = half * height;

    // Translate position to array coordinates
    point[X] += apothem_;
    point[Y] += apothem_;
    point[Z] += hheight;

    // Calculate initial cuboidal cell offset
    Real3 cell_pos;
    cell_pos[X] = point[X] / xy_span;
    cell_pos[Y] = point[Y] / xy_span;
    cell_pos[Z] = point[Z] / height;

    // Calculate initial nonnegative cuboidal X/Y/Z indices
    DodeCoord result;
    for (auto ax : {X, Y, Z})
    {
        result[ax] = std::max(0, static_cast<int>(std::floor(cell_pos[ax])));
    }

    // Calculate cuboidal array cell origin.
    // Cuboid spans from +/- apothem_ in x,y and +/- sqrt(2) * apothem_
    Real3 cuboidal_origin;
    cuboidal_origin[X] = result[X] * xy_span + apothem_;
    cuboidal_origin[Y] = result[Y] * xy_span + apothem_;
    cuboidal_origin[Z] = result[Z] * height + hheight;

    // Move cuboidal Z index into dodecahedral Z index
    result[Z] *= 2;

    // Update cell position relative to cuboidal origin
    cell_pos = point - cuboidal_origin;

    // Calculate scalar projections of point onto non-axis-aligned basis
    // vectors (a through d) to determine if its outside the actual cell's
    // boundaries. Test opposite-face pairs simultaneously by dotting with
    // positive normal and using absolute value. The winning direction must be
    // further away than the apothem.
    real_type max_proj = apothem_;
    DodeFace  max_face;

    for (auto test_ax : {Axes::a, Axes::b, Axes::c, Axes::d})
    {
        real_type proj
            = dot_product(dode_normal(DodeFace{test_ax, true}), cell_pos);
        real_type abs_proj = std::fabs(proj);

        if (abs_proj > max_proj)
        {
            max_proj = abs_proj;
            max_face = DodeFace{test_ax, proj > 0};
        }
    }

    if (max_face)
    {
        result += coordinate_update(max_face);
    }

    // Translate from maybe-out-of-bounds signed ints to in-bound unsigned
    auto result_uint = cell_indexer_.dims();

    // Clamp to in-bounds array cells
    for (auto ax : {X, Y, Z})
    {
        if (result[ax] < 0)
        {
            // Below the first coordinate
            CELER_ASSERT(result[ax] == -1);
            result_uint = 0;
        }
        else if (static_cast<size_type>(result[ax]) >= result_uint[ax])
        {
            // Above the last coordinate
            result_uint[ax] = result_uint[ax] - 1;
        }
        else
        {
            // Within normal range
            result_uint[ax] = static_cast<size_type>(result[ax]);
        }
    }
    return result_uint;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

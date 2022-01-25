//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/RectArrayProto.cc
 * \brief RectArrayProto class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RectArrayProto.hh"

#include <memory>
#include "base/Casts.hh"
#include "base/SoftEquivalence.hh"
#include "base/Definitions.hh"
#include "base/FixedViewArray.hh"
#include "base/Future.hh"
#include "base/Range.hh"
#include "orange/Fuzziness.hh"
#include "orange/query/RectArrayMetadata.hh"
#include "orange/track/RectArrayTracker.hh"
#include "CuboidShape.hh"
#include "PlacedShape.hh"

using make_unique;

namespace
{
//---------------------------------------------------------------------------//
/*!
 * Get the 'lower' corner of an array daughter.
 *
 * This is the translation value from the origin to the lower-left corner of a
 * cuboid shape.
 */
geometria::Real3 cuboid_proto_lower(const celeritas::Proto& proto)
{
    using namespace celeritas;
    using smart_cast;

    const Proto::RegionVec& interior = proto.interior();
    CELER_ASSERT(interior.size() == 1);
    CELER_ASSERT(interior.front().first == inside);

    // Extract the low-level shape and translate it as needed.
    const PlacedShape& ps  = *proto.interior().front().second;
    const CuboidShape& cub = smart_cast<const CuboidShape&>(*ps.shape());

    Real3 translation = cub.lower();
    translation += ps.transform().translation();

    return translation;
}

//---------------------------------------------------------------------------//
//! Bump the lower and upper bounds of a grid to reduce coincident surfaces.
RegularGrid make_bumped_grid(const RegularGrid& grid)
{
    CELER_EXPECT(grid);

    using Axis::x;
    using Axis::y;
    using Axis::z;

    const auto& fuzz = celeritas::fuzziness();

    // Copy vectors of edges
    auto edges = grid.edges();

    // Bump each one
    for (auto ax : {X, Y, Z})
    {
        real_type bump = std::max(
            fuzz.bump_abs(),
            fuzz.bump_rel() * (edges[ax].back() - edges[ax].front()));
        edges[ax].front() -= bump;
        edges[ax].back() += bump;
    }

    return {std::move(edges[X]), std::move(edges[Y]), std::move(edges[Z])};
}

} // namespace

namespace celeritas
{
//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
RectArrayProto::RectArrayProto(Params params) : md_(std::move(params.md))
{
    Insist(md_, "Missing array metadata");
    Insist(!params.units.empty(),
           "No units are defined for array '" << md_.name() << "'");

    const real_type tol = fuzziness().shape_enclosure_rel();

    // Calculate widths for each cell
    Grid_t::CellEdges widths;
    for (auto ax : {def::I, def::J, def::K})
    {
        widths[ax].assign(params.units.dims()[ax], 0.0);
    }

    for (auto i : range(widths[def::I].size()))
    {
        for (auto j : range(widths[def::J].size()))
        {
            for (auto k : range(widths[def::K].size()))
            {
                Grid_t::DimVector ijk{i, j, k};
                const auto&       sp_unit = params.units[ijk];
                if (!sp_unit)
                {
                    continue;
                }

                // Convert bounding shape of daughter to a cuboid
                const auto&        interior = sp_unit->interior();
                const CuboidShape* cuboid   = nullptr;
                if (interior.size() == 1 && interior.front().first == inside)
                {
                    const PlacedShape& placed = *interior.front().second;
                    CELER_VALIDATE(!placed.transform().has_rotation(),
                                   << "Entry {" << ijk << "} '"
                                   << sp_unit->metadata().name()
                                   << "' of array '" << md_.name()
                                   << "' cannot be rotated");
                    cuboid = dynamic_cast<const CuboidShape*>(
                        interior.front().second->shape().get());
                }
                CELER_VALIDATE(cuboid,
                               << "Entry {" << ijk << "} '"
                               << sp_unit->metadata().name() << "' of array '"
                               << md_.name() << "' is not a simple cuboid");

                // Compare or assign grid cell widths
                Real3 cell_widths = cuboid->upper() - cuboid->lower();
                for (auto ax : {def::I, def::J, def::K})
                {
                    auto* w = &widths[ax][ijk[ax]];
                    if (*w == 0.0)
                    {
                        // First time the cell has been encountered
                        *w = cell_widths[ax];
                    }
                    else
                    {
                        CELER_VALIDATE(soft_equiv(*w, cell_widths[ax], tol),
                                       << "Entry {" << ijk << "} '"
                                       << sp_unit->metadata().name()
                                       << "' of array '" << md_.name()
                                       << "' has an inconsistent width "
                                          "along the "
                                       << def::ijk_name(ax) << " axis");
                    }
                }
            }
        }
    }

    // Integrate widths to become edges, ensuring there's at least one cell per
    // IJK direction
    for (auto ax : {def::I, def::J, def::K})
    {
        real_type start = 0.0;
        for (auto& pos : widths[ax])
        {
            real_type cur_width = pos;
            CELER_VALIDATE(cur_width > 0,
                           << "Array '" << md_.name()
                           << "' has undefined units in a plane along the "
                           << def::ijk_name(ax) << " axis");
            pos = start;
            start += cur_width;
        }
        CELER_ASSERT(start > 0);
        widths[ax].push_back(start);
    }
    grid_ = Grid_t(std::move(widths[def::I]),
                   std::move(widths[def::J]),
                   std::move(widths[def::K]));

    // Transform IJK of input units to mesh cell indexing (KJI)
    units_.resize(grid_.num_volumes());
    for (auto i : range(grid_.num_cells_along(def::I)))
    {
        for (auto j : range(grid_.num_cells_along(def::J)))
        {
            for (auto k : range(grid_.num_cells_along(def::K)))
            {
                units_[grid_.index({i, j, k})]
                    = std::move(params.units[{i, j, k}]);
            }
        }
    }

    // Construct boundary shape
    PlacedShape::Params boundary_params;
    boundary_params.shape = std::make_shared<CuboidShape>(grid_.low_corner(),
                                                          grid_.high_corner());
    boundary_params.md    = md_;

    boundary_ = RegionVec(
        {{inside, std::make_shared<PlacedShape>(std::move(boundary_params))}});

    CELER_ENSURE(grid_);
    CELER_ENSURE(units_.size() == grid_.num_volumes());
    CELER_ENSURE(!boundary_.empty());
}

//---------------------------------------------------------------------------//
//! Default destructor
RectArrayProto::~RectArrayProto() = default;

//---------------------------------------------------------------------------//
/*!
 * Transformation to move the 'low' array corner to the specified point.
 *
 * This is the daughter-to-parent transform: destination in enclosing unit
 * (function argument) destination; and the low array corner is always the
 * origin. This method is used for KENO-V geometry as well as for VERA and XML
 * input.
 */
auto RectArrayProto::calc_placement(SpanConstReal3 pos) const -> Transform
{
    Real3 translation = make_vector(pos);
    // Don't need to subtract array "origin" since it's always (0, 0, 0)
    return Transform(translation);
}

//---------------------------------------------------------------------------//
/*!
 * Transformation to move the center of the specified unit.
 *
 * The transform is the daughter-to-parent transform that will place the
 * *center* of array unit \c ijk at the position \c pos in the parent unit.
 * This is not the same as the *origin* of the array unit.
 *
 * This method is needed for KENO-VI placement.
 */
auto RectArrayProto::calc_placement(DimVector ijk, SpanConstReal3 pos) const
    -> Transform
{
    CELER_EXPECT(ijk.all_lt(grid_.num_cells_dims()));
    CELER_EXPECT(units_[grid_.index(ijk)]);

    // Translation from the center of unit (it is required if the unit is
    // off-centered)
    Real3 translation = this->cuboid_proto_centroid(ijk);

    // Translation from array cell center to grid coordinates: -(hi + lo)/2
    translation -= (grid_.low_corner(ijk) + grid_.high_corner(ijk)) * half;

    // Array coordinates => enclosing parent unit
    translation += make_vector(pos);

    return translation;
}

//---------------------------------------------------------------------------//
/*!
 * Construct a universe from this proto
 */
auto RectArrayProto::build(BuildArgs args) const -> BuildResult
{
    CELER_EXPECT(args.implicit_boundary);

    BuildResult result;

    // Fill daughters
    for (auto k : range(grid_.num_cells_along(def::K)))
    {
        for (auto j : range(grid_.num_cells_along(def::J)))
        {
            for (auto i : range(grid_.num_cells_along(def::I)))
            {
                VolumeId            id(grid_.index({i, j, k}));
                const SPConstProto& sp_unit = units_[id.get()];
                if (!sp_unit)
                    continue;

                // Daughter => lower corner of array cell (0,0,0)
                Real3 translation = -cuboid_proto_lower(*sp_unit);

                // Low corner of array cell => array coordinates
                translation += grid_.low_corner({i, j, k});

                result.daughters.push_back(
                    {id, {sp_unit, Transform(translation)}});
            }
        }
    }

    // Create tracker with outer boundaries bumped outward to reduce coincident
    // surfaces
    result.tracker = make_unique<RectArrayTracker>(make_bumped_grid(grid_));

    // Create metadata
    RectArrayMetadata::Params md_params;
    md_params.unit = md_;
    md_params.dims = grid_.num_cells_dims();
    md_params.bbox = {grid_.low_corner(), grid_.high_corner()};
    result.md      = std::make_shared<RectArrayMetadata>(std::move(md_params));

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Get the boundary definition for defining a hole in a higher level
 */
auto RectArrayProto::interior() const -> const RegionVec&
{
    return boundary_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of an array daughter.
 *
 * Required translation is returned if the cuboid is off-center. Otherwise,
 * 0-translation '{0, 0, 0}' is returned.
 *
 */
Real3 RectArrayProto::cuboid_proto_centroid(DimVector ijk) const
{
    using namespace celeritas;
    using smart_cast;

    const auto& sp_unit = units_[grid_.index(ijk)];
    CELER_EXPECT(sp_unit);

    const Proto::RegionVec& interior = (*sp_unit).interior();
    CELER_ASSERT(interior.size() == 1);
    CELER_ASSERT(interior.front().first == inside);

    // Calculate the centroid of low-level shape and translate it as needed.
    const PlacedShape& ps  = *(*sp_unit).interior().front().second;
    const CuboidShape& cub = smart_cast<const CuboidShape&>(*ps.shape());

    Real3 centroid = (cub.lower() + cub.upper()) * half;
    centroid += ps.transform().translation();

    return centroid;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/UnitBuilder.cc
 * \brief UnitBuilder class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "UnitBuilder.hh"

#include "base/Range.hh"
#include "detail/ShapeBuilder.hh"
#include "PlacedShape.hh"

using celeritas::detail::ShapeBuilder;

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
UnitBuilder::UnitBuilder() : surfaces_(), insert_surface_(&surfaces_) {}

//---------------------------------------------------------------------------//
/*!
 * \brief Reserve space for regions
 */
void UnitBuilder::reserve(size_type num_regions)
{
    regions_.reserve(num_regions);
    cell_nodes_.reserve(num_regions);
    cell_md_.reserve(num_regions);
    volumes_.reserve(num_regions);
}

//---------------------------------------------------------------------------//
/*!
 * Define the exterior of the region by an implicit complement.
 *
 * The different zorder values have the following meanings:
 * - "MEDIA": in a well-connected universe so that no masking is needed to
 *   connect interior to exterior region
 * - "IMPLICIT_EXTERIOR": in a deeper-level universe where the exterior volume
 *   should be treated as the parent universe (excluded from search but not an
 *   error)
 * - "EXTERIOR": in a universe where the exterior volume has a higher
 *   precedence than other volumes (masked universe)
 */
void UnitBuilder::exterior(const VecSenseShape& interior,
                           ZOrder               zorder,
                           ObjectMetadata       md)
{
    CELER_EXPECT(!interior_bbox_);
    CELER_EXPECT(regions_.empty());
    CELER_EXPECT(zorder == ZOrder::media || zorder == ZOrder::implicit_exterior
                 || zorder == ZOrder::exterior);
    CELER_EXPECT(md);

    if (zorder != ZOrder::implicit_exterior)
    {
        // Create an actual outside cell
        this->add_volume(outside, interior, zorder, std::move(md), 0.0);

        interior_bbox_ = this->calc_bbox(interior);
    }
    else
    {
        // Add a temporary incomplete region
        UnitRegion incomplete;
        incomplete.zorder = zorder;
        regions_.push_back(incomplete);
        cell_nodes_.push_back({});
        cell_md_.emplace_back(std::move(md));
        volumes_.push_back(0.0);

        // Save interior for final build step
        implicit_interior_ = interior;
    }

    CELER_ENSURE(regions_.size() == 1 && cell_nodes_.size() == 1
                 && cell_md_.size() == 1 && volumes_.size() == 1);
}

//---------------------------------------------------------------------------//
/*!
 * Build a region from the given shape definition.
 *
 * This must happen *after* the exterior is built.
 */
VolumeId UnitBuilder::region(const VecSenseShape& interior,
                             ZOrder               zorder,
                             ObjectMetadata       md,
                             real_type            volume)
{
    CELER_EXPECT(zorder != ZOrder::invalid);
    CELER_EXPECT(md);
    CELER_EXPECT(volume >= 0);

    return this->add_volume(inside, interior, zorder, std::move(md), volume);
}

//---------------------------------------------------------------------------//
/*!
 * Complete construction
 */
auto UnitBuilder::operator()(ObjectMetadata unit_md) -> result_type
{
    CELER_EXPECT(unit_md);

    if (!implicit_interior_.empty())
    {
        // Save number of surfaces before adding exterior, so we can delete
        // implicit surfaces.
        const auto num_surfaces = surfaces_.size();
        // Add "exterior" cell with fake zorder
        this->add_volume(
            outside, implicit_interior_, ZOrder::media, cell_md_.front(), 0.0);
        auto outside_node = cell_nodes_.back();

        tree_.replace(outside_node, false);
        cell_nodes_.front() = outside_node;

        // The first (exterior) placeholder cell is replaced by the one we just
        // added. Remove the exterior temporrary.
        regions_.pop_back();
        cell_nodes_.pop_back();
        volumes_.pop_back();
        cell_md_.pop_back();

        // Calculate physical bounding box for the outside
        interior_bbox_ = this->calc_bbox(implicit_interior_);

        // Revert any temporary surfaces from the exterior that don't actually
        // matter.
        surfaces_.resize(num_surfaces);
        surface_md_.resize(num_surfaces);
    }
    else if (!interior_bbox_)
    {
        // If no exterior was given explicitly
        // TODO this might be only for unit tests, so we may want to eliminate
        interior_bbox_ = geometria::infinite_bbox();
    }

    CELER_ASSERT(interior_bbox_);

    CELER_ASSERT(cell_md_.size() == regions_.size()
                 && cell_nodes_.size() == regions_.size()
                 && volumes_.size() == regions_.size());

    // Construct interior regions
    for (auto i : range(cell_nodes_.size()))
    {
        CELER_ASSERT(cell_nodes_[i]);
        regions_[i].interior = tree_.build_cell(cell_nodes_[i]);
    }

    // Construct metadata
    CELER_ASSERT(surface_md_.size() == surfaces_.size());
    UnitMetadata::Params md_params;
    md_params.unit = std::move(unit_md);
    md_params.surfaces
        = UnitMetadata::TableShapeFace::from_container(std::move(surface_md_));
    md_params.cells     = std::move(cell_md_);
    md_params.bbox      = interior_bbox_;
    md_params.volumes   = std::move(volumes_);
    md_params.is_simple = this->is_simple();
    // TODO: add CSG tree to md params

    // Sort individual rows so the order is reproducible
    auto sort_shape_face = [](const ShapeFace& lhs, const ShapeFace& rhs) {
        return std::tie(lhs.first->name(), lhs.second)
               < std::tie(rhs.first->name(), rhs.second);
    };
    for (auto i : range(md_params.surfaces.num_rows()))
    {
        auto row_view = make_span(md_params.surfaces[i]);
        std::sort(row_view.begin(), row_view.end(), sort_shape_face);
    }

    result_type result;
    result.surfaces = std::move(surfaces_);
    result.regions  = std::move(regions_);
    result.md       = std::make_shared<UnitMetadata>(std::move(md_params));
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Whether all meaningful regions have the same Z order.
 *
 * If the first region (exterior) is "implicit", it's not going to be used for
 * tracking, and its zorder is *ignored*.
 *
 * This function is public primarily for unit testing. The UnitProto builder
 * can extract this important attribute from \c built.md->is_simple.
 */
bool UnitBuilder::is_simple() const
{
    CELER_EXPECT(!regions_.empty());
    auto region = regions_.begin();
    if (regions_.size() > 1 && region->zorder == ZOrder::implicit_exterior)
    {
        // Ignore the exterior
        ++region;
    }

    return region->zorder == regions_.back().zorder;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Calculate bounding box, adding surfaces to metadata.
 *
 * This is a simplified version of add_volume with an 'inside' sense.
 */
BoundingBox UnitBuilder::calc_bbox(const VecSenseShape& interior)
{
    // Construct "interior" for bounding box: surfaces should already be
    // added (note )
    CSGTree              temp_tree;
    detail::ShapeBuilder build_shape{insert_surface_, temp_tree};
    build_shape.push(CSGCell::LOGIC_AND, Transform{});
    for (const auto& sense_shape : interior)
    {
        build_shape.push(CSGCell::LOGIC_AND, sense_shape.second->transform());
        sense_shape.second->shape()->build(build_shape);
        build_shape.pop(sense_shape.first);
    }
    build_shape.pop();

    return build_shape().bbox;
}

//---------------------------------------------------------------------------//
/*!
 * Build a region from the given shape definition.
 */
VolumeId UnitBuilder::add_volume(Sense                sense,
                                 const VecSenseShape& interior,
                                 ZOrder               zorder,
                                 ObjectMetadata       md,
                                 real_type            volume)
{
    CELER_EXPECT(regions_.empty() || regions_.back().zorder >= zorder);
    CELER_EXPECT(
        std::all_of(interior.begin(), interior.end(), [](const Halfspace& ss) {
            return static_cast<bool>(ss.second);
        }));
    CELER_EXPECT(zorder != ZOrder::invalid);
    CELER_EXPECT(volume >= 0);
    CELER_EXPECT(md);

    detail::ShapeBuilder build_shape{insert_surface_, tree_};
    bool                 has_internal_surfaces = false;

    // Begin cell
    build_shape.push(CSGCell::LOGIC_AND, Transform{});

    for (const auto& sense_shape : interior)
    {
        // Add callback for associating shape and surface extension with
        // resulting shape ID
        auto mark_surface
            = [this, &sense_shape](SurfaceId id, std::string name) {
                  this->add_surface(id, sense_shape.second, std::move(name));
              };
        build_shape.set_surface_callback(mark_surface);

        // Apply the local shape's transformation, 'and'ing all daughters
        build_shape.push(CSGCell::LOGIC_AND, sense_shape.second->transform());

        // Tell unplaced shape to build surfaces
        sense_shape.second->shape()->build(build_shape);

        // Shapes 'AND' together all the surfaces they build
        build_shape.pop(sense_shape.first);

        // Update 'has internal'
        has_internal_surfaces = has_internal_surfaces
                                || sense_shape.first == outside
                                || !sense_shape.second->shape()->is_convex();
    }

    // Region is 'AND'-ed shapes, sense is usually inside unless 'exterior'
    build_shape.pop(sense);

    // Build CSG node and bbox
    auto built = build_shape();

    if (sense == outside && !tree_.at(built.cell).is_leaf())
    {
        // Outside of a multi-surface node, so likely to have internal
        // surface crossings
        has_internal_surfaces = true;
    }

    // Result
    UnitRegion region;
    region.zorder                = zorder;
    region.has_internal_surfaces = has_internal_surfaces;
    region.bbox                  = std::move(built.bbox);

    // Move to our lists of regions and metadata
    cell_nodes_.push_back(built.cell);
    regions_.push_back(std::move(region));
    cell_md_.push_back(std::move(md));
    volumes_.push_back(volume);

    CELER_ENSURE(regions_.size() == cell_md_.size());
    return VolumeId(regions_.size() - 1);
}

//---------------------------------------------------------------------------//
/*!
 * Add a surface.
 */
void UnitBuilder::add_surface(SurfaceId    surface,
                              SPConstShape shape,
                              std::string  ext)
{
    CELER_EXPECT(surface);
    CELER_EXPECT(shape);
    CELER_EXPECT(!ext.empty());
    CELER_EXPECT(surface.get() <= surface_md_.size());

    if (surface.get() >= surface_md_.size())
    {
        surface_md_.resize(surface.get() + 1);
    }
    surface_md_[surface.get()].insert({std::move(shape), std::move(ext)});
}

//---------------------------------------------------------------------------//
} // namespace celeritas

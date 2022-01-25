//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file query/UnitMetadata.cc
 * \brief UnitMetadata class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "UnitMetadata.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "base/Casts.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "orange/construct/CSGCell.hh"
#include "orange/construct/UnitRegion.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "orange/track/UnitTracker.hh"

namespace celeritas
{
namespace
{
//---------------------------------------------------------------------------//
// HELPER FUNCTIONS/CLASSES
//---------------------------------------------------------------------------//
struct Describe
{
    std::ostream& os;

    template<class S>
    void operator()(S surf)
    {
        this->os << surf;
    }
};

//---------------------------------------------------------------------------//
char zorder_to_char(ZOrder zo)
{
    switch (zo)
    {
        case ZOrder::invalid:
            return '?';
        case ZOrder::media:
            return 'M';
        case ZOrder::hole:
            return 'H';
        case ZOrder::implicit_exterior:
            return ' ';
        case ZOrder::exterior:
            return 'B';
    };
    return '?';
}

//---------------------------------------------------------------------------//
} // end namespace

//---------------------------------------------------------------------------//
/*!
 * Construct from cell/surface metadata
 */
UnitMetadata::UnitMetadata(Params params) : md_(std::move(params))
{
    CELER_EXPECT(md_.unit);
#ifdef REQUIRE_ON
    for (auto i : range(md_.surfaces.num_rows()))
    {
        CELER_EXPECT(!md_.surfaces[i].empty());
    }
#endif
    CELER_EXPECT(!md_.cells.empty());
    CELER_EXPECT(std::all_of(
        md_.cells.begin(), md_.cells.end(), [](const ObjectMetadata& obj_md) {
            return bool(obj_md);
        }));
    CELER_EXPECT(md_.bbox);
    CELER_EXPECT(md_.volumes.size() == md_.cells.size() || md_.volumes.empty());

    if (md_.volumes.empty())
    {
        md_.volumes.assign(md_.cells.size(), 0.0);
    }
}

//---------------------------------------------------------------------------//
//! Virtual destructor
UnitMetadata::~UnitMetadata() = default;

//---------------------------------------------------------------------------//
/*!
 * Get metadata about this universe
 */
const ObjectMetadata& UnitMetadata::metadata() const
{
    return md_.unit;
}

//---------------------------------------------------------------------------//
/*!
 * Local bounding box (extents) of this universe
 */
auto UnitMetadata::bbox() const -> BoundingBox
{
    return md_.bbox;
}

//---------------------------------------------------------------------------//
/*!
 * Number of local surfaces
 */
auto UnitMetadata::num_surfaces() const -> surface_int
{
    return md_.surfaces.size();
}

//---------------------------------------------------------------------------//
/*!
 * Number of local cells
 */
auto UnitMetadata::num_volumes() const -> volume_int
{
    return md_.cells.size();
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local surface ID into a user-facing surface label
 */
std::string UnitMetadata::id_to_label(SurfaceId id) const
{
    CELER_EXPECT(id < md_.surfaces.size());

    const ShapeFace& shape_face = this->surface_md(id).front();
    std::string      name       = shape_face.first->name();
    name.push_back('.');
    name += shape_face.second;
    return name;
}

//---------------------------------------------------------------------------//
/*!
 * Convert a local cell ID into a user-facing cell label
 */
std::string UnitMetadata::id_to_label(VolumeId id) const
{
    CELER_EXPECT(id < md_.cells.size());
    return md_.cells[id.get()].name();
}

//---------------------------------------------------------------------------//
/*!
 * Get the volume of a cell
 */
real_type UnitMetadata::volume(VolumeId id) const
{
    CELER_EXPECT(id < md_.volumes.size());
    real_type result = md_.volumes[id.get()];
    CELER_ENSURE(result >= 0);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Set the volume of a cell
 */
void UnitMetadata::set_volume(VolumeId id, real_type vol)
{
    CELER_EXPECT(id < md_.volumes.size());
    CELER_EXPECT(vol >= 0);

    md_.volumes[id.get()] = vol;
}

//---------------------------------------------------------------------------//
/*!
 * Describe using data from the corresponding tracker
 */
void UnitMetadata::describe(std::ostream& os, const Tracker& base_tracker) const
{
    const auto& tracker = smart_cast<const UnitTracker&>(base_tracker);
    using std::left;
    using std::setw;

    constexpr int max_cols = 60;
    std::string   long_sep(max_cols + 2, '=');
    long_sep.front() = ' ';
    long_sep.back()  = '\n';

    //// CELLS ////

    // Create separator for name
    size_type name_len = 4;
    {
        auto longest_object_md = std::max_element(
            md_.cells.begin(),
            md_.cells.end(),
            [](const ObjectMetadata& lhs, const ObjectMetadata& rhs) {
                return lhs.name().size() < rhs.name().size();
            });
        name_len = std::max(longest_object_md->name().size(), name_len);
    }
    std::string sep(name_len, '=');

    std::string zorder_sep = "";
    if (!this->is_simple())
    {
        zorder_sep = "= ";
    }

    // clang-format off
    os << "======= " << zorder_sep << sep << long_sep
       << "Cell    " << (this->is_simple() ? "" : "Z ")
                     << left << setw(name_len) << "Name"
                     << " Surface logic\n"
       << "======= " << zorder_sep << sep << long_sep;

    std::string newline_str(1 + 8 + zorder_sep.size() + name_len + 1, ' ');
    newline_str.front() = '\n';

    std::vector<std::string> internal_surfaces;

    for (VolumeId::size_type c : range(tracker.num_volumes()))
    {
        UnitRegion reg = tracker.get_region(VolumeId{c});
        const ObjectMetadata& vol_md = md_.cells[c];

        os << setw(7) << c;
        if (!this->is_simple())
        {
            os << ' ' << zorder_to_char(reg.zorder);
        }
        os << ' ' << setw(name_len) << vol_md.name() << ' ';

        // Write tokens
        size_type columns = 0;
        std::string str;
        bool prepend_space = false;
        for (CSGCell::logic_int token : reg.interior.logic())
        {
            str.resize(1);
            switch (token)
            {
                case CSGCell::LOGIC_TRUE: str.front() = '*'; break;
                case CSGCell::LOGIC_NOT: str.front() = '~'; break;
                case CSGCell::LOGIC_AND: str.front() = '&'; break;
                case CSGCell::LOGIC_OR : str.front() = '|'; break;
                default: str = std::to_string(token);
            }
            if (prepend_space)
            {
                ++columns;
            }
            if (columns + str.size() >= max_cols)
            {
                os << newline_str;
                columns = 0;
                prepend_space = false;
            }
            if (prepend_space)
            {
                os << ' ';
            }
            os << str;
            columns += str.size();
            prepend_space = (columns < max_cols);
        }
        os << newline_str << "(from ``" << *vol_md.provenance() << "``)\n";

        if (reg.has_internal_surfaces)
        {
            internal_surfaces.push_back(vol_md.name());
        }
    }
    os << "======= " << zorder_sep << sep << long_sep
       << "\n";
    // clang-format on

    if (!internal_surfaces.empty())
    {
        os << "Cells with reentrant surface tracking: \""
           << join(internal_surfaces.begin(), internal_surfaces.end(), "\", \"")
           << "\"\n\n";
        internal_surfaces.clear();
    }

    //// SURFACES ////

    // Construct surface labels and get name lengths
    name_len = 4;
    for (SurfaceId::size_type s : range(this->num_surfaces()))
    {
        for (const auto& shape_face : this->surface_md(SurfaceId{s}))
        {
            name_len = std::max(name_len,
                                shape_face.first->name().size() + 1
                                    + shape_face.second.size());
        }
    }
    sep.assign(name_len, '=');

    os << "======= " << sep << long_sep << "Surface " << left << setw(name_len)
       << "Name"
       << " Description\n"
       << "======= " << sep << long_sep;

    newline_str.assign(1 + 8, ' ');
    newline_str.front() = '\n';

    auto describe_surface
        = make_surface_action(tracker.surfaces(), Describe{os});

    for (SurfaceId::size_type s : range(tracker.num_surfaces()))
    {
        os << setw(7) << s << setw(name_len + 2) << ' ';
        describe_surface(SurfaceId{s});

        for (const auto& shape_face : this->surface_md(SurfaceId{s}))
        {
            std::string name = shape_face.first->name() + '.'
                               + shape_face.second;
            os << newline_str << setw(name_len) << name << " (from ``"
               << *shape_face.first->metadata().provenance() << "``)";
        }
        os << '\n';
    }
    os << "======= " << sep << long_sep << '\n';
}

//---------------------------------------------------------------------------//
} // namespace celeritas

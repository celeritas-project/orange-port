//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/test/ShapeTest.cc
 * \brief ShapeTest class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ShapeTest.hh"

#include <iomanip>
#include "base/FixedViewArray.hh"
#include "base/StringFunctions.hh"
#include "orange/surfaces/SurfaceAction.hh"
#include "../detail/SurfaceInserter.hh"
#include "../detail/ShapeBuilder.hh"
#include "../CSGTree.hh"
#include "../Shape.hh"

namespace
{
struct SurfaceTypeStr
{
    template<class S>
    const char* operator()(S surf) const
    {
        return celeritas::to_cstring(S::surface_type());
    }
};

struct SurfaceData
{
    template<class S>
    std::vector<real_type> operator()(S surf) const
    {
        auto generic_data = surf.view().data;
        return {generic_data.begin(), generic_data.end()};
    }
};
} // namespace

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Build a shape (with transformation)
 */
void ShapeTest::build(const Shape& shape, Transform t)
{
    this->surfaces = {};
    detail::SurfaceInserter         insert_surface(&this->surfaces);
    CSGTree                         tree;
    celeritas::detail::ShapeBuilder build_shape(insert_surface, tree);

    // Add callback for saving the surface name
    this->surface_names.clear();
    build_shape.set_surface_callback([this](SurfaceId, std::string name) {
        this->surface_names.push_back(std::move(name));
    });

    // Add the user-provided shape: inside all given surfaces
    build_shape.push(CSGCell::LOGIC_AND, t);
    shape.build(build_shape);
    build_shape.pop();

    // Save result
    auto result = build_shape();
    this->cell  = tree.build_cell(result.cell);
    this->bbox  = std::move(result.bbox);
}

//---------------------------------------------------------------------------//
/*!
 * Print expected results.
 */
void ShapeTest::print_expected() const
{
    using std::cout;
    using std::endl;

    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "EXPECT_VEC_EQ(CSGCell::from_string(\"" << cell
         << "\").logic(), this->cell.logic());\n\n"
         << "EXPECT_VEC_SOFT_EQ(Real3(" << make_fixed_view(bbox.lower())
         << "), this->bbox.lower());\n"
         << "EXPECT_VEC_SOFT_EQ(Real3(" << make_fixed_view(bbox.upper())
         << "), this->bbox.upper());\n\n"
         << "static const char* const expected_surface_names[] = "
         << to_string(surface_names) << ";\n"
         << "EXPECT_VEC_EQ(expected_surface_names, this->surface_names);\n\n";

    auto get_surface_data = make_surface_action(this->surfaces, SurfaceData{});
    auto get_surface_typestr
        = make_surface_action(this->surfaces, SurfaceTypeStr{});

    cout << "ASSERT_EQ(" << this->surfaces.size()
         << ", this->surfaces.size());\n";
    for (auto id : this->surfaces.all_ids())
    {
        auto        data      = get_surface_data(id);
        const char* surf_type = get_surface_typestr(id);
        cout << "EXPECT_VEC_SOFT_EQ(VecDbl({" << std::setprecision(14)
             << join(data.begin(), data.end(), ", ") << "}), "
             << "this->get_surface_data(SurfaceType::" << surf_type << ", "
             << id.get() << "));\n";
    }
    cout << "/*** END CODE ***/\n";
}

//---------------------------------------------------------------------------//
/*!
 * Get the surface data, checking for the correct surface type
 */
std::vector<real_type>
ShapeTest::get_surface_data(SurfaceType type, surfid_int idx) const
{
    SCOPED_TRACE("Testing surface " + std::to_string(idx));

    EXPECT_LT(idx, this->surfaces.size());

    SurfaceId id(idx);
    EXPECT_EQ(this->surfaces.get_type(id), type);
    auto get_surface_data = make_surface_action(this->surfaces, SurfaceData{});
    auto data             = get_surface_data(id);
    return data;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

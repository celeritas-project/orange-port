//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/XMLUtils.cc
 * \brief XMLUtils class definitions
 * \note   Copyright (c) 2021 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "XMLUtils.hh"

#include <unordered_map>
#include <TeuchosParameterList.hpp>
#include "Nemesis/database/PlistUtilities.hh"
#include "base/Assert.hh"
#include "base/Constants.hh"
#include "orange/TransformUtils.hh"
#include "../ConeShape.hh"
#include "../CuboidShape.hh"
#include "../CylinderSegmentShape.hh"
#include "../CylinderShape.hh"
#include "../ECylinderShape.hh"
#include "../EllipsoidShape.hh"
#include "../GeneralQuadricShape.hh"
#include "../HopperShape.hh"
#include "../ParallelepipedShape.hh"
#include "../PlacedShape.hh"
#include "../PlaneShape.hh"
#include "../PrismShape.hh"
#include "../RhombicDodecahedronShape.hh"
#include "../RightTetrahedronShape.hh"
#include "../RingShape.hh"
#include "../SlabShape.hh"
#include "../SphereShape.hh"
#include "../TriangularPrismShape.hh"
#include "../WedgeShape.hh"

using Teuchos::ParameterList;
using ArrInt       = Teuchos::Array<int>;
using ArrDbl       = Teuchos::Array<real_type>;
using ArrStr       = Teuchos::Array<std::string>;
using SPConstShape = std::shared_ptr<const celeritas::Shape>;

namespace celeritas
{
namespace detail
{
namespace
{
//---------------------------------------------------------------------------//
// ANONYMOUS HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Build a cuboid
 */
SPConstShape build_cuboid(const ParameterList& plist)
{
    Real3 lo{plist.get<real_type>("xmin"),
             plist.get<real_type>("ymin"),
             plist.get<real_type>("zmin")};
    Real3 hi{plist.get<real_type>("xmax"),
             plist.get<real_type>("ymax"),
             plist.get<real_type>("zmax")};
    return std::make_shared<CuboidShape>(lo, hi);
}

//---------------------------------------------------------------------------//
/*!
 * Build a sphere
 */
SPConstShape build_sphere(const ParameterList& plist)
{
    return std::make_shared<SphereShape>(plist.get<real_type>("radius"));
}

//---------------------------------------------------------------------------//
/*!
 * Build a cylinder
 */
SPConstShape build_cyl(const ParameterList& plist)
{
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad cylinder extents size " << extents.size());
    auto axis = static_cast<Axis>(plist.get<int>("axis"));

    return std::make_shared<CylinderShape>(
        axis, plist.get<real_type>("radius"), extents[0], extents[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a cylinder segment
 */
SPConstShape build_cylsegment(const ParameterList& plist)
{
    using constants::two_pi;

    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad cylsegment extents size " << extents.size());

    real_type begin_angle = plist.get<real_type>("angle") * two_pi;
    real_type arc         = plist.get<real_type>("arc") * two_pi;

    return std::make_shared<CylinderSegmentShape>(plist.get<real_type>("inner_"
                                                                       "radiu"
                                                                       "s"),
                                                  plist.get<real_type>("outer_"
                                                                       "radiu"
                                                                       "s"),
                                                  begin_angle,
                                                  arc,
                                                  extents[0],
                                                  extents[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a cylindrical shell
 */
SPConstShape build_ring(const ParameterList& plist)
{
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad ring extents size " << extents.size());

    return std::make_shared<RingShape>(plist.get<real_type>("inner_radius"),
                                       plist.get<real_type>("outer_radius"),
                                       extents[0],
                                       extents[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a prism
 */
SPConstShape build_prism(const ParameterList& plist)
{
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad prism extents size " << extents.size());

    return std::make_shared<PrismShape>(plist.get<int>("num_sides"),
                                        plist.get<real_type>("apothem"),
                                        plist.get<real_type>("rotfrac"),
                                        extents[0],
                                        extents[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a slab
 */
SPConstShape build_slab(const ParameterList& plist)
{
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad slab extents size " << extents.size());
    auto axis = static_cast<Axis>(plist.get<int>("axis"));

    return std::make_shared<SlabShape>(axis, extents[0], extents[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a plane
 */
SPConstShape build_plane(const ParameterList& plist)
{
    return std::make_shared<PlaneShape>(get_space_vector(plist, "normal"),
                                        get_space_vector(plist, "point"));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a wedge
 */
SPConstShape build_wedge(const ParameterList& plist)
{
    auto corner = plist.get<ArrDbl>("corner_pt");
    CELER_VALIDATE(corner.size() == 2,
                   << "Bad wedge corner size " << corner.size());

    return std::make_shared<WedgeShape>(
        plist.get<real_type>("width"), // x
                                       // length
        corner[0],
        corner[1],                     // x, y base point
        plist.get<real_type>("height") // z
                                       // length
    );
}

//---------------------------------------------------------------------------//
/*!
 * Build a cone
 */
SPConstShape build_cone(const ParameterList& plist)
{
    auto radii = plist.get<ArrDbl>("radii");
    CELER_VALIDATE(radii.size() == 2,
                   << "Bad cone radii size " << radii.size());
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad cone extents size " << extents.size());
    auto axis = static_cast<Axis>(plist.get<int>("axis"));

    using PairDbl = ConeShape::PairDbl;
    return std::make_shared<ConeShape>(
        axis, PairDbl{radii[0], radii[1]}, PairDbl{extents[0], extents[1]});
}

//---------------------------------------------------------------------------//
/*!
 * Build an ellipsoid
 */
SPConstShape build_ellipsoid(const ParameterList& plist)
{
    auto radii = plist.get<ArrDbl>("radii");
    CELER_VALIDATE(radii.size() == 3,
                   << "Bad ellipsoid radii size " << radii.size());

    return std::make_shared<EllipsoidShape>(
        Real3(radii[0], radii[1], radii[2]));
}

//---------------------------------------------------------------------------//
/*!
 * Build a hopper
 */
SPConstShape build_hopper(const ParameterList& plist)
{
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad hopper extents size " << extents.size());
    auto lower = plist.get<ArrDbl>("lower_pt");
    CELER_VALIDATE(lower.size() == 2,
                   << "Bad hopper lower size " << lower.size());
    auto upper = plist.get<ArrDbl>("upper_pt");
    CELER_VALIDATE(upper.size() == 2,
                   << "Bad hopper upper size " << upper.size());

    return std::make_shared<HopperShape>(
        upper[0], upper[1], extents[1], lower[0], lower[1], extents[0]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a right tetrahedron
 */
SPConstShape build_righttet(const ParameterList& plist)
{
    auto lengths = plist.get<ArrDbl>("lengths");
    CELER_VALIDATE(lengths.size() == 3,
                   << "Bad righttet lengths size " << lengths.size());

    return std::make_shared<RightTetrahedronShape>(
        lengths[0], lengths[1], lengths[2]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a triprism
 */
SPConstShape build_triprism(const ParameterList& plist)
{
    auto lengths = plist.get<ArrDbl>("lengths");
    CELER_VALIDATE(lengths.size() == 3,
                   << "Bad triprism lengths size " << lengths.size());

    return std::make_shared<TriangularPrismShape>(
        lengths[0], lengths[1], lengths[2]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a elliptical cylinder
 */
SPConstShape build_ecylinder(const ParameterList& plist)
{
    auto radii = plist.get<ArrDbl>("radii");
    CELER_VALIDATE(radii.size() == 2,
                   << "Bad cone radii size " << radii.size());
    auto extents = plist.get<ArrDbl>("extents");
    CELER_VALIDATE(extents.size() == 2,
                   << "Bad cone extents size " << extents.size());

    return std::make_shared<ECylinderShape>(
        ECylinderShape::PlaneVector{radii[0], radii[1]}, extents[0], extents[1]);
}

//---------------------------------------------------------------------------//
/*!
 * Build a rhombic dodecahedron
 */
SPConstShape build_rhombdod(const ParameterList& plist)
{
    auto apothem = plist.get<real_type>("apothem");
    return std::make_shared<RhombicDodecahedronShape>(apothem);
}

//---------------------------------------------------------------------------//
/*!
 * Build a parallelepiped
 */
SPConstShape build_ppiped(const ParameterList& plist)
{
    auto      length = get_space_vector(plist, "lengths");
    real_type psi    = plist.get<real_type>("psi");
    real_type theta  = plist.get<real_type>("theta");
    real_type phi    = plist.get<real_type>("phi");
    return std::make_shared<ParallelepipedShape>(length, psi, theta, phi);
}

//---------------------------------------------------------------------------//
/*!
 * Build a parallelepiped
 */
SPConstShape build_quadric(const ParameterList& plist)
{
    auto      second = get_space_vector(plist, "second");
    auto      cross  = get_space_vector(plist, "cross");
    auto      first  = get_space_vector(plist, "first");
    real_type scalar = plist.get<real_type>("scalar");

    return std::make_shared<GeneralQuadricShape>(second[0],
                                                 second[1],
                                                 second[2],
                                                 cross[0],
                                                 cross[1],
                                                 cross[2],
                                                 first[0],
                                                 first[1],
                                                 first[2],
                                                 scalar);
}
//---------------------------------------------------------------------------//
} // namespace

//---------------------------------------------------------------------------//
/*!
 * Construct metadata from a plist.
 *
 * This uses the "name" attribute that's assumed to be there, and checks for
 * the _provenance value that's exported by omnutils.
 */
ObjectMetadata build_md(const ParameterList& plist)
{
    ObjectMetadata::Params params;
    params.name        = plist.get<std::string>("name");
    params.description = get_default(plist, "description", std::string{});
    if (plist.isParameter("_provenance"))
    {
        params.provenance
            = Provenance::from_string(plist.get<std::string>("_provenance"));
    }
    else
    {
        params.provenance = Provenance::from_unknown();
    }
    return ObjectMetadata{std::move(params)};
}

//---------------------------------------------------------------------------//
/*!
 * Construct a shape from the given plist.
 */
std::shared_ptr<const PlacedShape> build_shape(const ParameterList& plist)
{
    using ShapeBuilder = SPConstShape (*)(const ParameterList&);

    //! Map of shape "type" -> build function
    static const std::unordered_map<std::string, ShapeBuilder> shape_builders = {
        {"cuboid", &build_cuboid},
        {"sphere", &build_sphere},
        {"cyl", &build_cyl},
        {"cylsegment", &build_cylsegment},
        {"ring", &build_ring},
        {"prism", &build_prism},
        {"slab", &build_slab},
        {"plane", &build_plane},
        {"wedge", &build_wedge},
        {"cone", &build_cone},
        {"ellipsoid", &build_ellipsoid},
        {"hopper", &build_hopper},
        {"righttet", &build_righttet},
        {"triprism", &build_triprism},
        {"ecylinder", &build_ecylinder},
        {"rhombdod", &build_rhombdod},
        {"ppiped", &build_ppiped},
        {"quadric", &build_quadric},
    };

    // Build the shape: use the "_type" parameter to determine what function to
    // call.
    auto iter = shape_builders.find(plist.get<std::string>("_type"));
    CELER_VALIDATE(iter != shape_builders.end(),
                   << "Invalid shape type '" << plist.get<std::string>("_type")
                   << "'");

    // Call member function pointer for the shape building method
    ShapeBuilder build_impl = iter->second;
    CELER_ASSERT(iter->second);
    auto shape = build_impl(plist);
    CELER_ASSERT(shape);

    PlacedShape::Params params;
    params.shape = std::move(shape);
    if (plist.isParameter("translate"))
    {
        params.transform.translation(get_space_vector(plist, "translate"));
    }
    if (plist.isParameter("rotate"))
    {
        SpaceMatrix mat = get_space_matrix(plist, "rotate");
        params.transform.rotation(mat);
    }
    params.md = build_md(plist);

    return std::make_shared<PlacedShape>(std::move(params));
}

//---------------------------------------------------------------------------//
/*!
 * Load a Real3 from a plist
 */
Real3 get_space_vector(const ParameterList& plist, const std::string& key)
{
    auto vals = plist.get<Teuchos::Array<real_type>>(key);
    Insist(vals.size() == 3,
           "Bad size " << vals.size() << " for '" << key << "' vector");
    return Real3(vals[0], vals[1], vals[2]);
}

//---------------------------------------------------------------------------//
/*!
 * Load a DimVector from a plist
 */
Array<def::size_type, 3>
get_dim_vector(const ParameterList& plist, const std::string& key)
{
    auto vals = plist.get<Teuchos::Array<int>>(key);
    Insist(vals.size() == 3,
           "Bad size " << vals.size() << " for '" << key << "' vector");
    return Array<def::size_type, 3>(vals[0], vals[1], vals[2]);
}

//---------------------------------------------------------------------------//
/*!
 * Load a SpaceMatrix from a plist
 */
SpaceMatrix get_space_matrix(const ParameterList& plist, const std::string& key)
{
    auto vals = plist.get<Teuchos::Array<real_type>>(key);
    Insist(vals.size() == 9,
           "Bad size " << vals.size() << " for '" << key << "' matrix");

    SpaceMatrix mat;

    // Load rotation matrix in row-major order and convert to column-major
    // order so that the mat is dimensioned [j][i]
    auto itr = vals.begin();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j, ++itr)
        {
            mat[j][i] = *itr;
        }
    }

    // Enforce unitariness of rotation matrix
    geometria::orthonormalize(mat);

    return mat;
}

//---------------------------------------------------------------------------//
/*!
 * Construct a Transform from a plist
 */
Transform build_transform(const ParameterList& plist)
{
    Transform transform;

    if (plist.isParameter("translate"))
    {
        transform.translation(get_space_vector(plist, "translate"));
    }

    if (plist.isParameter("rotate"))
    {
        SpaceMatrix mat = get_space_matrix(plist, "rotate");
        transform.rotation(mat);
    }
    return transform;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

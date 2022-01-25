//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/ShapeBuilder.cc
 * \brief ShapeBuilder class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ShapeBuilder.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>

#include "base/Assert.hh"
#include "Nemesis/comm/Logger.hh"
#include "base/Definitions.hh"
#include "base/Join.hh"
#include "base/Range.hh"
#include "base/StringFunctions.hh"
#include "base/VectorFunctions.hh"

#include "orange/Fuzziness.hh"
#include "orange/surfaces/PlaneAligned.hh"
#include "orange/surfaces/CylAligned.hh"
#include "orange/surfaces/Plane.hh"
#include "orange/surfaces/CylCentered.hh"
#include "orange/surfaces/CenteredSphere.hh"
#include "orange/surfaces/Sphere.hh"
#include "orange/surfaces/ConeAligned.hh"
#include "orange/surfaces/SimpleQuadric.hh"
#include "orange/surfaces/GeneralQuadric.hh"
#include "SurfaceInserter.hh"
// #include "PossiblyFudged.hh"

using soft_zero;
using Axis::x;
using Axis::y;
using Axis::z;
using geometria::infinite_bbox;
namespace celeritas
{
namespace detail
{
namespace
{
bool is_neg_inf(real_type v)
{
    return std::isinf(v) && std::signbit(v);
}
bool is_pos_inf(real_type v)
{
    return std::isinf(v) && !std::signbit(v);
}
} // namespace
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * Constructor.
 *
 * The top level starts with infinite bounding boxes and no transforms.
 */
ShapeBuilder::ShapeBuilder(SurfaceInserter& surfaces, CSGTree& tree)
    : insert_surface_(surfaces)
    , tree_(tree)
    , stack_{{CSGCell::LOGIC_AND,
              VecDaughter{},   // No daughters (all space)
              Transform{},     // No local transform
              infinite_bbox(), // All space
              Transform{},
              infinite_bbox()}}
    , num_planes_(0)
{
    const auto& fuzz    = ::celeritas::fuzziness();
    elision_abs_        = fuzz.surface_elision_abs();
    simplification_abs_ = fuzz.surface_simplification_abs();
}

//---------------------------------------------------------------------------//
/*!
 * Set callback for notifying of an added surface
 */
void ShapeBuilder::set_surface_callback(AddSurfaceCallback cb)
{
    added_surface_ = std::move(cb);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Join daughters under a new local tranform.
 *
 * This pushes the new transform on the *right* side of the existing transform.
 * In other words, if the existing transform is \f$ T_0 \f$ and the transform
 * being pushed is \f$ T_1 \f$, the transform that will be applied to new
 * surfaces is
 * \f[ T' = T_0 T_1 \f]
 *
 * \f[ T'' = T_0 T_1 T_2 = T' T_2\f]
 *
 * The \c Transform.transform method applies a new transform to the *left*, so
 * we must do global.transform(local) rather than local.transform(global)
 */
void ShapeBuilder::push(LogicToken conjunction, const Transform& local)
{
    // Construct new global transform
    Transform global = local;
    global.transform(this->transform());

    // Bounding box is 'infinite' and continually contracting if ANDing, or
    // 'empty' and continually expanding if ORing
    BoundingBox default_bbox
        = (conjunction == CSGCell::LOGIC_AND ? infinite_bbox() : BoundingBox{});

    // Push the transform and a new infinite bbox
    stack_.push_back(
        {conjunction, VecDaughter{}, local, default_bbox, global, default_bbox});
}

//---------------------------------------------------------------------------//
/*!
 * Pop the last pushed transform.
 *
 * This transforms the bounding box and clips the previous "no transform"
 * bbox.
 */
void ShapeBuilder::pop(Sense sense)
{
    CELER_EXPECT(stack_.size() > 1);

    // Save and pop the last transform and node
    auto frame = std::move(stack_.back());
    stack_.pop_back();

    // Pop the CSG node: empty back means infinite cell (unusual).
    auto node_id = tree_.emplace(std::move(frame.nodes), frame.conjunction);

    // Update parent frame nodes
    auto& parent_frame = stack_.back();
    parent_frame.nodes.push_back({sense, node_id});

    // Update parent bounding boxes
    if (sense == inside)
    {
        // Intersect the new 'local' and 'global' bboxes with that of the
        // just-built daughter cell
        auto transformed_bbox = frame.bbox.calc_transformed(frame.transform);

        if (parent_frame.conjunction == CSGCell::LOGIC_AND)
        {
            // Clip the current untranslated bbox
            parent_frame.bbox.clip(transformed_bbox);
            parent_frame.global_bbox.clip(frame.global_bbox);
        }
        else if (parent_frame.conjunction == CSGCell::LOGIC_OR)
        {
            parent_frame.bbox = parent_frame.bbox.calc_union(transformed_bbox);
            parent_frame.global_bbox
                = parent_frame.global_bbox.calc_union(frame.global_bbox);
        }
        else
        {
            CELER_ASSERT_UNREACHABLE();
        }
    }
    else if (parent_frame.conjunction == CSGCell::LOGIC_OR)
    {
        // Assume 'OR' with 'outside a volume' is infinite; not strictly true
        // for half-spaces but it's a conservative guess
        parent_frame.bbox        = infinite_bbox();
        parent_frame.global_bbox = infinite_bbox();
    }
}

//---------------------------------------------------------------------------//
// SURFACES
//---------------------------------------------------------------------------//
/*!
 * Build a sphere.
 */
void ShapeBuilder::sphere(Sense sense, SpanConstReal3 point, real_type radius)
{
    CELER_EXPECT(radius >= 0);

    // Create sphere, do rotate/translate about origin
    return this->construct("s", sense, Sphere(point, radius));
}

//---------------------------------------------------------------------------//
/*!
 * Build a plane normal to one axis (in the unrotated frame)
 *
 * Given sense is the sense of the new shape_Builder with respect to the new
 * surface:
 * for example, a cuboid is pos of the lower X coordinate and neg of the upper
 * X coordinate.
 *
 * This is used for cuboids and other objects.
 *
 * \param[in,out] surfaces SurfaceContainer object to add the plane to
 * \param[in]     axis     Unrotated axis on which the plane lies
 * \param[in]     Sense    Sense of the shape_Builder with respect to this
 * plane \param[in]     loc      Unrotated, untranslated plane location along
 * axis
 */
void ShapeBuilder::plane(Axis axis, Sense sense, real_type loc)
{
    CELER_EXPECT(axis < def::END_XYZ);

    if (sense == pos && is_neg_inf(loc))
    {
        // everything above -inf in local frame
        return;
    }
    else if (sense == neg && is_pos_inf(loc))
    {
        // everything below +inf
        return;
    }

    // Set name whether it's on the "plus" or "minus" side of the
    // shape, so the sense is negated
    std::string face_name{"??"};
    face_name[0] = (sense == neg ? 'p' : 'm');
    face_name[1] = to_cstring(axis)[0];

    switch (axis)
    {
        case X:
            return this->construct(std::move(face_name), sense, PlaneX(loc));
        case Y:
            return this->construct(std::move(face_name), sense, PlaneY(loc));
        case Z:
            return this->construct(std::move(face_name), sense, PlaneZ(loc));
        default:
            CELER_ASSERT_UNREACHABLE();
    }
}

//---------------------------------------------------------------------------//
/*!
 * Build a cylinder along the (untranslated) center
 *
 * \param[in]     axis     Center axis of the cylinder
 * \param[in]     Sense    neg for inside, pos for outside
 * \param[in]     radius   Radius of the cylinder
 */
void ShapeBuilder::cyl(Axis axis, Sense sense, real_type radius)
{
    CELER_EXPECT(axis < def::END_XYZ);

    // Set name whether it's on the "outside" or "inside"
    std::string face_name{"c??"};
    face_name[1] = (sense == neg ? 'o' : 'i');
    face_name[2] = to_cstring(axis)[0];

    switch (axis)
    {
        case X:
            return this->construct(std::move(face_name), sense, CCylX(radius));
        case Y:
            return this->construct(std::move(face_name), sense, CCylY(radius));
        case Z:
            return this->construct(std::move(face_name), sense, CCylZ(radius));
        default:
            CELER_ASSERT_UNREACHABLE();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a non-orthogonal, outward-facing plane.
 *
 * \param[in,out] surfaces SurfaceContainer object to add the plane to
 * \param[in]     normal   Outward-facing normal (unnormalized)
 * \param[in]     point    A point on the surface
 */
void ShapeBuilder::plane(SpanConstReal3 normal, SpanConstReal3 point)
{
    // Non-orthogonal plane counter
    std::string face_name{"p" + std::to_string(num_planes_++)};

    // Build transformed plane
    return this->construct(std::move(face_name), inside, Plane(normal, point));
}

//---------------------------------------------------------------------------//
/*!
 * Build a cone, given the origin (vanishing point) and tangent angle
 *        of opening.
 */
void ShapeBuilder::cone(Axis           axis,
                        Sense          sense,
                        SpanConstReal3 origin,
                        real_type      tangent)
{
    CELER_EXPECT(axis < def::END_XYZ);

    // Non-orthogonal plane counter
    std::string face_name{"k??"};
    face_name[1] = (sense == neg ? 'o' : 'i');
    face_name[2] = to_cstring(axis)[0];

    switch (axis)
    {
        case X:
            return this->construct(
                std::move(face_name), sense, ConeX(origin, tangent));
        case Y:
            return this->construct(
                std::move(face_name), sense, ConeY(origin, tangent));
        case Z:
            return this->construct(
                std::move(face_name), sense, ConeZ(origin, tangent));
        default:
            CELER_ASSERT_UNREACHABLE();
    }
}

//---------------------------------------------------------------------------//
/*!
 * Build a simple quadric (ellipsoid, hyperbola, ....)
 */
void ShapeBuilder::simple_quadric(SpanConstReal3 abc,
                                  SpanConstReal3 def,
                                  real_type      g,
                                  SpanConstReal3 origin)
{
    // Build and transform quadric
    return this->construct("sq", inside, SimpleQuadric(abc, def, g, origin));
}

//---------------------------------------------------------------------------//
/*!
 * Build a quadric using the general form equation
 *
 * I.e., \f[ ax2 + by2 + cz2 + dxy + exz + fyz + gx + hy + iz + j = 0) \f]
 */
void ShapeBuilder::general_quadric(SpanConstReal3 abc,
                                   SpanConstReal3 def,
                                   SpanConstReal3 ghi,
                                   real_type      j,
                                   Sense          sense)
{
    return this->construct("gq", sense, GeneralQuadric(abc, def, ghi, j));
}

//---------------------------------------------------------------------------//
/*!
 * Build and return the constructed bounding box and CSG node.
 */
auto ShapeBuilder::operator()() -> result_type
{
    CELER_EXPECT(stack_.size() == 1);
    CELER_EXPECT(stack_.front().nodes.size() == 1);

    result_type result;
    result.cell = tree_.emplace(stack_.front().nodes.front());
    result.bbox = stack_.front().global_bbox;
    result.bbox.clip(stack_.front().bbox);
    return result;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Insert a transformed surface
 *
 * This function allows for specialization on surface types for surface
 * simplification after a rotation.
 */
template<class S>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         S           surf)
{
    // Default: just do insertion
    return this->construct_final(std::move(face_name), sense, surf);
}

//---------------------------------------------------------------------------//
/*!
 * Insert a transformed ortho plane
 *
 * Ortho planes should be snapped to the origin when appropriate
 */
template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         PlaneX      surf)
{
    // Check for 'snap' threshold with plane
    real_type d = surf.position();
    if (soft_zero(d, elision_abs_))
    {
        d = 0;
    }

    return this->construct_final(std::move(face_name), sense, PlaneX(d));
}

template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         PlaneY      surf)
{
    // Check for 'snap' threshold with plane
    real_type d = surf.position();
    if (soft_zero(d, elision_abs_))
    {
        d = 0;
    }

    return this->construct_final(std::move(face_name), sense, PlaneY(d));
}

template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         PlaneZ      surf)
{
    // Check for 'snap' threshold with plane
    real_type d = surf.position();
    if (soft_zero(d, elision_abs_))
    {
        d = 0;
    }

    return this->construct_final(std::move(face_name), sense, PlaneZ(d));
}

//---------------------------------------------------------------------------//
#define GG_TRY_SIMPLIFY(CLS, SURF)                   \
    if (CLS::can_simplify(SURF))                     \
    {                                                \
        return this->construct_transformed(          \
            std::move(face_name), sense, CLS(SURF)); \
    }

/*!
 * Insert a transformed plane directly
 *
 * This function allows for specialization on surface types for surface
 * simplification after a rotation.
 */
template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         Plane       plane)
{
    const real_type tol = simplification_abs_;

    Real3     n(plane.normal());
    real_type d(plane.displacement());

    // Snap nearly-zero normals to zero
    for (auto ax : range(3u))
    {
        if (soft_zero(n[ax], tol))
            n[ax] = 0;
    }

    // To prevent opposite-value planes from being defined but not
    // deduplicated, ensure the first non-zero normal component is in the
    // positive half-space. This also takes care of flipping orthogonal planes
    // defined like {-x = 3}, translating them to { x = -3 }.
    for (auto ax : range(3u))
    {
        if (n[ax] > 0)
        {
            break;
        }
        else if (n[ax] < 0)
        {
            // Multiply any remaining nonzero axes by -1
            // (previous axes are zero so just skip them)
            // (don't multiply zeros by -1, otherwise we get -0 values)
            for (auto ax2 : range(ax, 3u))
            {
                if (n[ax2] != 0)
                    n[ax2] *= -1;
            }
            // Multiply nonzero distance by -1
            if (d != 0)
            {
                d *= -1;
            }
            sense = flip_sense(sense);
            break;
        }
    }

    // Construct updated plane with "snapped" normal directions and with a
    // positive normal, updating and renormalizing displacement in the process
    Plane p(n, d);

    // See if it's an orthogonal plane
    GG_TRY_SIMPLIFY(PlaneX, p);
    GG_TRY_SIMPLIFY(PlaneY, p);
    GG_TRY_SIMPLIFY(PlaneZ, p);

    // No simplification occurred. Before the final insertion, see if the
    // displacement needs to be snapped to the origin.
    if (soft_zero(d, elision_abs_))
    {
        p = Plane(n, 0);
    }

    return this->construct_final(std::move(face_name), sense, p);
}

//---------------------------------------------------------------------------//
/*!
 * Insert a sphere
 */
template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         Sphere      sph)
{
    GG_TRY_SIMPLIFY(CenteredSphere, sph);

    // No simplification occurred
    return this->construct_final(std::move(face_name), sense, sph);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert an axis-aligned cylinder
 */
template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         CylX        cyl)
{
    GG_TRY_SIMPLIFY(CCylX, cyl);

    // No simplification occurred
    return this->construct_final(std::move(face_name), sense, cyl);
}

template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         CylY        cyl)
{
    GG_TRY_SIMPLIFY(CCylY, cyl);

    // No simplification occurred
    return this->construct_final(std::move(face_name), sense, cyl);
}

template<>
void ShapeBuilder::construct_transformed(std::string face_name,
                                         Sense       sense,
                                         CylZ        cyl)
{
    GG_TRY_SIMPLIFY(CCylZ, cyl);

    // No simplification occurred
    return this->construct_final(std::move(face_name), sense, cyl);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert a simple quadric
 *
 * Simple quadrics can be automatically simplified to cylinders, cones, and the
 * general plane.
 */
template<>
void ShapeBuilder::construct_transformed(std::string   face_name,
                                         Sense         sense,
                                         SimpleQuadric sq)
{
    GG_TRY_SIMPLIFY(CylX, sq);
    GG_TRY_SIMPLIFY(CylY, sq);
    GG_TRY_SIMPLIFY(CylZ, sq);
    GG_TRY_SIMPLIFY(ConeX, sq);
    GG_TRY_SIMPLIFY(ConeY, sq);
    GG_TRY_SIMPLIFY(ConeZ, sq);
    GG_TRY_SIMPLIFY(Plane, sq);

    // No simplification occurred
    return this->construct_final(std::move(face_name), sense, sq);
}

//---------------------------------------------------------------------------//
/*!
 * Insert a general quadric
 *
 * General quadrics can be automatically simplified to other surface types.
 */
template<>
void ShapeBuilder::construct_transformed(std::string    face_name,
                                         Sense          sense,
                                         GeneralQuadric gq)
{
    // Try to simplify into a quadric
    GG_TRY_SIMPLIFY(SimpleQuadric, gq);

    // No simplification occurred
    return this->construct_final(std::move(face_name), sense, gq);
}

//---------------------------------------------------------------------------//
/*!
 * Insert a transformed surface
 *
 * This function allows for specialization on surface types for surface
 * simplification after a rotation.
 */
template<class S>
void ShapeBuilder::construct_final(std::string face_name, Sense sense, S surf)
{
    // Insert surface and add name to metadata
    auto surf_id = insert_surface_(surf);
    if (added_surface_)
    {
        added_surface_(surf_id, std::move(face_name));
    }

    // Add the constructed surface handle to the CSG tree
    // Flip the sense because "inside" a surface is false (determinant < 0) and
    // the resulting expression must be true.
    stack_.back().nodes.push_back({flip_sense(sense), tree_.emplace(surf_id)});

    // Clip post-transformed bbox, using a possibly de-duplicated surface
    // instead of the given surface.
    CELER_ASSERT(insert_surface_.surfaces().is_type<S>(surf_id));
    insert_surface_.surfaces().get<S>(surf_id).clip(sense, this->global_bbox());
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

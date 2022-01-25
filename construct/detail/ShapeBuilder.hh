//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/detail/ShapeBuilder.hh
 * \brief ShapeBuilder class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <functional>
#include <string>
#include <vector>
#include "orange/BoundingBox.hh"
#include "orange/Transform.hh"
#include "orange/Definitions.hh"
#include "orange/surfaces/Definitions.hh"
#include "../CSGCell.hh"
#include "../CSGTree.hh"

namespace celeritas
{
namespace detail
{
class SurfaceInserter;
//---------------------------------------------------------------------------//
/*!
 * Build a cell out of shapes.
 */
class ShapeBuilder
{
  public:
    //@{
    //! Public type aliases
    using surfid_int         = SurfaceId::size_type;
    using LogicToken         = CSGCell::LogicToken;
    using NodeId             = CSGTree::NodeId;
    using AddSurfaceCallback = std::function<void(SurfaceId, std::string)>;
    //@}

    struct result_type
    {
        NodeId      cell; //!< CSG Tree node from creating cell
        BoundingBox bbox; //!< Bounds in post-transform coordinate system
    };

  public:
    // Constructor
    ShapeBuilder(SurfaceInserter& surfaces, CSGTree& tree);

    // Set callback for notifying of an added  surface
    void set_surface_callback(AddSurfaceCallback cb);

    // Begin a new shape combining daughters with a local transform
    void push(LogicToken conjunction, const Transform& t);

    // End this shape, optionally negate by passing in Sense::outside
    void pop(Sense sense = inside);

    //// SURFACES ////

    // Build the given surface type
    template<class S>
    inline void surface(std::string face_name, Sense sense, S&& surf);

    // TODO: deprecate these methods in favor of manual surface insertion?

    // Build a sphere at the given point
    void sphere(Sense sense, SpanConstReal3 origin, real_type radius);

    // Build a plane normal to the given axis
    void plane(Axis axis, Sense sense, real_type loc);

    // Build a non-orthogonal, "outward"-facing plane
    void plane(SpanConstReal3 normal, SpanConstReal3 point);

    // Build a cylinder along the given axis
    void cyl(Axis axis, Sense sense, real_type radius);

    // Build a cone, given the vanishing point and tangent angle of opening
    void cone(Axis axis, Sense sense, SpanConstReal3 origin, real_type tangent);

    // Build an axis-aligned quadric
    void simple_quadric(SpanConstReal3 second_order,
                        SpanConstReal3 first_order,
                        real_type      constant,
                        SpanConstReal3 origin);

    // Build a quadric using the general equation form
    void general_quadric(SpanConstReal3 second_order,
                         SpanConstReal3 cross_terms,
                         SpanConstReal3 first_order,
                         real_type      constant,
                         Sense          sense);

    //// BBOX MODIFICATION ////

    // Reduce the pre-transformation bounding box size
    inline void clip_bbox_lower(Axis axis, real_type value);
    // Reduce the pre-transformation bounding box size
    inline void clip_bbox_upper(Axis axis, real_type value);

    //// COMPLETE ////

    // Build. Leaves this class in an empty/invalid state.
    result_type operator()();

  private:
    // Insert a surface (dispatches based on transform/translate)
    template<class Surface>
    inline void construct(std::string face_name, Sense sense, Surface surf);

    // Insert a transformed surface (has template overloads)
    template<class Surface>
    void
    construct_transformed(std::string face_name, Sense sense, Surface surf);

    // Insert a shape (no template overloads; actually does insertion)
    template<class Surface>
    void construct_final(std::string face_name, Sense sense, Surface surf);

    // Current cell's "non-transformed" bounding box
    inline BoundingBox& local_bbox();

    // Current cell's "global" reference frame bounding box
    inline BoundingBox& global_bbox();

    // Global transform
    inline const Transform& transform() const;

  private:
    //// TYPES ////

    // Vector of {sense, node ID}
    using VecDaughter = CSGNode::VecDaughter;

    struct LocalTransform
    {
        LogicToken  conjunction;      //!< How to combine daughters
        VecDaughter nodes;            //!< Currently added daughters
        Transform   transform;        //!< In this stack frame
        BoundingBox bbox;             //!< In this stack frame
        Transform   global_transform; //!< Combined with all lower frames
        BoundingBox global_bbox;      //!< In global frame
    };

    //// DATA ////

    //! SurfaceContainer inside the current unit
    SurfaceInserter& insert_surface_;

    //! CSG nodes inside the current unit
    CSGTree& tree_;

    // Callback for notifying of an added surface
    AddSurfaceCallback added_surface_;

    // Stack of bbox, transforms, etc. (depth = cell builder depth)
    std::vector<LocalTransform> stack_;

    // Global bounding box used for pre-transform surface truncation
    BoundingBox global_bbox_;

    // Tolerances
    real_type elision_abs_;
    real_type simplification_abs_;

    // Counter for plane surface labeling
    int num_planes_;
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
#include "ShapeBuilder.i.hh"
//---------------------------------------------------------------------------//

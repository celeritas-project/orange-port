//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/ShapeContainer.hh
 * \brief ShapeContainer class declaration
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
#pragma once

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include "PlacedShape.hh"
#include "Shape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Helper class for creating and storing shapes.
 *
 * This is mostly used for unit tests but could also be used for construction
 * helper classes.
 *
 * Differences from a regular map:
 * - Emplace syntax actually constructs the daughter shape *and* the placed
 *   shape in one go. The first two arguments must be the metadata/transform
 *   arguments to the placed shape. The 'key' is the name provided in the
 *   metadata.
 * - Insert uses the placed shape's name (from the metadata).
 * - Insertion must succeed (unique key given): Insist will throw if not.
 * - Insert/emplace return a const reference instead of an iterator/bool pair.
 * - operator[] is constant and will raise an Insist on an invalid key
 * - No mutable iterators: shapes are constant once inserted.
 */
class ShapeContainer
{
  public:
    //@{
    //! Public type aliases
    using key_type    = std::string;
    using mapped_type = std::shared_ptr<const PlacedShape>;

    using Map_t          = std::map<key_type, mapped_type>;
    using const_iterator = Map_t::const_iterator;
    using size_type      = Map_t::size_type;
    //@}

  public:
    // Constructor
    ShapeContainer() = default;

    // Add a shape, return insertion success
    inline const mapped_type& insert(mapped_type shape);

    //! Create and insert a placed shape with the given shape primitive
    // (inlined here to work around intel 14 compiler bug)
    template<class S, typename... Args>
    const mapped_type&
    emplace(ObjectMetadata md, Transform transform, Args&&... args)
    {
        static_assert(std::is_base_of<Shape, S>::value,
                      "Emplace template must be a Shape subclass");

        PlacedShape::Params params;
        params.md        = std::move(md);
        params.transform = std::move(transform);
        params.shape     = std::make_shared<S>(std::forward<Args>(args)...);

        return this->insert(std::make_shared<PlacedShape>(std::move(params)));
    }

    // Find the shape with the given name
    inline const_iterator find(const key_type& name) const;

    // Get the shape with the given name
    inline const mapped_type& at(const key_type& name) const;

    // Get the shape with the given name (it must exist!)
    inline const mapped_type& operator[](const key_type& name) const;

    //@{
    //! Accessors
    bool                  empty() const { return shapes_.empty(); }
    size_type             size() const { return shapes_.size(); }
    Map_t::const_iterator begin() const { return shapes_.begin(); }
    Map_t::const_iterator cbegin() const { return begin(); }
    Map_t::const_iterator end() const { return shapes_.end(); }
    Map_t::const_iterator cend() const { return end(); }
    //@}

  private:
    //// DATA ////

    std::map<key_type, mapped_type> shapes_;
};

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
std::ostream& operator<<(std::ostream& os, const ShapeContainer&);

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "ShapeContainer.i.hh"
//---------------------------------------------------------------------------//

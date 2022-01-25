//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file construct/PlacedShape.cc
 * \brief PlacedShape class definitions
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PlacedShape.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with fixed parameters
 */
PlacedShape::PlacedShape(Params params) : data_(std::move(params))
{
    CELER_EXPECT(data_.shape);

    if (!data_.md)
    {
        // Give an arbitrary but unique name
        // TODO: add a 'hash' method to shapes for deduplication and use that
        // instead?
        static int         shape_counter = 0;
        std::ostringstream os;
        os << "Untitled " << data_.shape->type() << ' ' << shape_counter;
        ObjectMetadata::Params md;
        md.name       = os.str();
        md.provenance = Provenance::from_unknown();
        data_.md      = ObjectMetadata(std::move(md));
    }

    CELER_ENSURE(this->shape());
    CELER_ENSURE(!this->name().empty());
    CELER_ENSURE(this->metadata());
}

//---------------------------------------------------------------------------//
/*!
 * Return a transformed version of this, moving to parent reference
 */
PlacedShape PlacedShape::transformed(const Transform& t) const
{
    Params params = data_;
    params.transform.transform(t);
    return PlacedShape{std::move(params)};
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Print a shape
 */
std::ostream& operator<<(std::ostream& os, const PlacedShape& shape)
{
    os << "{'" << shape.name() << "': " << *shape.shape() << " with transform "
       << shape.transform() << " from " << *shape.metadata().provenance()
       << '}';
    return os;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file surfaces/detail/SolveQuadratic.hh
 * \brief Quadric solvers
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#pragma once

#include "orange/Definitions.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//

// Find all positive (nonzero) roots for general quadric surfaces
real_type* solve_quadratic_general(
    bool on_surface, real_type a, real_type b_2, real_type c, real_type* x);

// Find all positive (nonzero) roots of x^2 + b*x + c = 0.
real_type* solve_quadratic(real_type b_2, real_type c, real_type* x);

// Find all positive roots of x^2 + b*x + c = 0 where x = 0 is a root ("on
// surface").
real_type* solve_degenerate_quadratic1(real_type b_2, real_type* x);

// Find all positive roots of a*x^2 + b*x + c = 0, for small/zero a ("along
// surface")
real_type*
solve_degenerate_quadratic2(real_type b_a_2, real_type c_a, real_type* x);

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas

//---------------------------------------------------------------------------//

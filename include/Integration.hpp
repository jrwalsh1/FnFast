//------------------------------------------------------------------------------
/// \file Integration.hpp
//
// Author(s):
//    Jon Walsh
//
// Copyright:
//    Copyright (C) 2015  LBL
//
//    This file is part of the EFTofLSS library. EFTofLSS is distributed under the
//    terms of the GNU General Public License version 3 (GPLv3), see the COPYING
//    file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Defintion of objects related to integration in CUBA
//------------------------------------------------------------------------------

#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include <vector>
#include <functional>

using namespace std;

//------------------------------------------------------------------------------
/**
 * \struct IntegralResult
 *
 * \brief Defines container to hold the results of an integral.
 *
 * Contains the integral result, error, and probability that the error is not robust
 */
//------------------------------------------------------------------------------
struct IntegralResult
{
   double result;       ///< result of the integral
   double error;        ///< error of the integral
   double prob;         ///< probability that the error is NOT a reliable estimate

   IntegralResult(double res, double err, double p) : result(res), error(err), prob(p) {}
};

#endif // INTEGRATION_HPP

//------------------------------------------------------------------------------
/// \file Integration.hpp
//
// Author(s):
//    Jon Walsh
//
// Copyright:
//    Copyright (C) 2015  LBL
//
//    This file is part of the FnFast library. FnFast is distributed under the
//    terms of the GNU General Public License version 3 (GPLv3), see the COPYING
//    file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Defintion of objects related to integration in CUBA
//------------------------------------------------------------------------------

#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include "cuba.h"

namespace fnfast {

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

//------------------------------------------------------------------------------
/**
 * \struct VEGASintegrator
 *
 * \brief Defines a simple interface to the Cuba VEGAS integration routine
 *
 * Returns the integration result into a IntegralResult container
 */
//------------------------------------------------------------------------------
struct VEGASintegrator
{
   int ndim;                  ///< number of dimensions in the integral
   double epsrel;             ///< relative accuracy desired
   int maxeval;               ///< maximum number of integrand evaluations
   int nstart;                ///< number of initial integrand evaluations
   int nincrease;             ///< number of integrand evaluations to increment by
   int nbatch;                ///< batch size for PS point sampling

   /// constructor
   VEGASintegrator(int numdim, double err = 1e-3, int neval = 500000, int numstart = 1000, int numincrease = 1000, int numbatch = 1000)
   : ndim(numdim), epsrel(err), maxeval(neval), nstart(numstart), nincrease(numincrease), nbatch(numbatch) {}

   /// integration function
   IntegralResult integrate(integrand_t integrand, void * userdata);
};

} // namespace fnfast

#endif // INTEGRATION_HPP

//------------------------------------------------------------------------------
/// \file Integration.cpp
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
//    Implementation of class Integration
//------------------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "Integration.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
IntegralResult VEGASintegrator::integrate(integrand_t integrand, void * userdata)
{
   // VEGAS integration parameters

   // PARAMETER: phase space dimensionality set by ndim
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // PARAMETER: relative precision set by epsrel
   const double epsabs = 0;
   // min, max number of points
   const int mineval = 0;
   // PARAMETER: maximum number of integrand calls set by maxeval
   // PARAMETER: starting number of points set by nstart
   // PARAMETER: number of additional pts sampled per iteration set by nincrease
   // PARAMETER: batch size to sample pts in set by nbatch
   // grid number
   // 1-10 saves the grid for another integration
   const int gridnum = 0;
   // file for the state of the integration
   const char *statefile = NULL;
   // spin
   void* spin = NULL;
   // random number seed
   const int vegasseed = 37;
   // flags:
   // bits 0&1: verbosity level
   // bit 2: whether or not to use only last sample (0 for all samps, 1 for last only)
   // bit 3: whether or not to use sharp edges in importance function (0 for no, 1 for yes)
   // bit 4: retain the state file (0 for no, 1 for yes)
   // bits 8-31: random number generator, also uses seed parameter:
   //    seed = 0: Sobol (quasi-random) used, ignores bits 8-31 of flags
   //    seed > 0, bits 8-31 of flags = 0: Mersenne Twister
   //    seed > 0, bits 8-31 of flags > 0: Ranlux
   // current flag setting: 1038 = 10000001110
   int flags = 1038;
   // number of regions, evaluations, fail code
   int neval, fail;

   // containers for output
   double integral[ncomp], error[ncomp], prob[ncomp];

   // run VEGAS
   Vegas(ndim, ncomp, integrand, userdata, nvec,
       epsrel, epsabs, flags, vegasseed,
       mineval, maxeval, nstart, nincrease, nbatch,
       gridnum, statefile, spin,
       &neval, &fail, integral, error, prob);

   // save the results in a container
   IntegralResult result(integral[0], error[0], prob[0]);

   return result;
}

} // namespace fnfast

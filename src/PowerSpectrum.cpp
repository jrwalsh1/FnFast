//------------------------------------------------------------------------------
/// \file PowerSpectrum.cpp
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
//    Implementation of class PowerSpectrum
//------------------------------------------------------------------------------

#include <iostream>

#include "PowerSpectrum.hpp"
#include "cuba.h"

//------------------------------------------------------------------------------
PowerSpectrum::PowerSpectrum(Order order)
: _order(order), _diagrams(DiagramSet2point(order)), _UVcutoff(10.)
{}

//------------------------------------------------------------------------------
double PowerSpectrum::tree(double k, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // set the external momenta
   ThreeVector k2(0, 0, k);
   MomentumMap<ThreeVector> momenta(unordered_map<MomentumLabel, ThreeVector> {{MomentumLabel::k1, -k2}, {MomentumLabel::k2, k2}});

   return _diagrams.value_tree(momenta, kernels, PL);
}

//------------------------------------------------------------------------------
double PowerSpectrum::oneLoop_excl(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const
{
   return _diagrams.value_oneLoop(mom, kernels, PL);
}

//------------------------------------------------------------------------------
double PowerSpectrum::twoLoop_excl(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const
{
   return _diagrams.value_twoLoop(mom, kernels, PL);
}

//------------------------------------------------------------------------------
IntegralResult PowerSpectrum::oneLoop(double k, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // integration method
   OneLoopIntegrator integrator(k, _UVcutoff, &kernels, PL, this);

   // Integration
   // VEGAS parameters
   // PS dimensionality
   // q: 2-dim
   const int ndim = 2;
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // absolute uncertainty (safeguard)
   const double epsrel = 1e-4;
   const double epsabs = 0;
   // min, max number of points
   const int mineval = 0;
   const int maxeval = 250000;
   // starting number of points
   const int nstart = 1000;
   // increment per iteration
   // number of additional pts sampled per iteration
   const int nincrease = 1000;
   // batch size to sample PS points in
   const int nbatch = 1000;
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
   Vegas(ndim, ncomp, oneLoop_integrand, &integrator, nvec,
       epsrel, epsabs, flags, vegasseed,
       mineval, maxeval, nstart, nincrease, nbatch,
       gridnum, statefile, spin,
       &neval, &fail, integral, error, prob);

   // save the results in a container
   IntegralResult result(integral[0], error[0], prob[0]);

   return result;
}

//------------------------------------------------------------------------------
pair<double, MomentumMap<ThreeVector>* const> PowerSpectrum::OneLoopIntegrator::generate_point(double qpts[2])
{
   // we sample q flat in spherical coordinates
   // q components
   double qmag = qpts[0] * qmax;
   double qcosth = 2 * qpts[1] - 1.;

   // jacobian
   // qmax from the magnitude integral,
   // 2 from the cos theta jacobian,
   // pick up a 2pi from the phi integral,
   // and a 1/(2pi)^3 from the measure
   double jacobian = qmag * qmag * qmax / (2 * pi*pi);

   // 3-vector for the loop momentum
   momenta[MomentumLabel::q] = ThreeVector(qmag * sqrt(1. - qcosth*qcosth), 0, qmag * qcosth);

   return pair<double, MomentumMap<ThreeVector>* const>(jacobian, &momenta);
}

//------------------------------------------------------------------------------
int PowerSpectrum::oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the integrator
   OneLoopIntegrator* integrator = static_cast<OneLoopIntegrator*>(userdata);

   // generate the PS points
   double qpts[2] = {xx[0], xx[1]};
   pair<double, MomentumMap<ThreeVector>* const> PSpoint = integrator->generate_point(qpts);

   // calculate the integrand
   double integrand = PSpoint.first * (integrator->powerspectrum->diagrams()->value_oneLoop(*(PSpoint.second), *(integrator->kernels), integrator->PL));

   ff[0] = integrand;

   return 0;
}

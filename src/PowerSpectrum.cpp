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

namespace fnfast {

//------------------------------------------------------------------------------
/*DAN*/
PowerSpectrum::PowerSpectrum(Order order) : _order(order), _diagrams(DiagramSet2pointSPT(_order)), _EFTdiagrams(DiagramSet2pointEFT(_EFTorder(_order))), _UVcutoff(10.) {}

//------------------------------------------------------------------------------
double PowerSpectrum::tree(double k, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // set the external momenta
   ThreeVector k2(0, 0, k);
   LabelMap<Momentum, ThreeVector> momenta {{Momentum::k1, -k2}, {Momentum::k2, k2}};

   return _diagrams.value_tree(momenta, kernels, PL);
}

//------------------------------------------------------------------------------
IntegralResult PowerSpectrum::oneLoop(double k, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // integration method
   LoopPhaseSpace phasespace(k, _UVcutoff, &kernels, PL, this);

   // VEGAS integration via cuba
   VEGASintegrator vegas(2);

   return vegas.integrate(oneLoop_integrand, &phasespace);
}
   
//------------------------------------------------------------------------------
/*DAN*/
double PowerSpectrum::treeEFT(double k, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // set the external momenta
   ThreeVector k2(0, 0, k);
   LabelMap<Momentum, ThreeVector> momenta {{Momentum::k1, -k2}, {Momentum::k2, k2}};
      
   return _EFTdiagrams.value_tree(momenta, kernels, PL);
}
   
//------------------------------------------------------------------------------
std::pair<double, LabelMap<Momentum, ThreeVector>* const> PowerSpectrum::LoopPhaseSpace::generate_point_oneLoop(std::vector<double> xpts)
{
   // we sample q flat in spherical coordinates
   // q components
   double qmag = xpts[0] * qmax;
   double qcosth = 2 * xpts[1] - 1.;

   // jacobian
   // qmax from the magnitude integral,
   // 2 from the cos theta jacobian,
   // pick up a 2pi from the phi integral,
   // and a 1/(2pi)^3 from the measure
   double jacobian = qmag * qmag * qmax / (2 * pi*pi);

   // 3-vector for the loop momentum
   momenta[Momentum::q] = ThreeVector(qmag * sqrt(1. - qcosth*qcosth), 0, qmag * qcosth);

   return std::pair<double, LabelMap<Momentum, ThreeVector>* const>(jacobian, &momenta);
}

//------------------------------------------------------------------------------
std::pair<double, LabelMap<Momentum, ThreeVector>* const> PowerSpectrum::LoopPhaseSpace::generate_point_twoLoop(std::vector<double> xpts)
{
   // we sample q1, q2 flat in spherical coordinates
   // q1, q2 components
   double q1mag = xpts[0] * qmax;
   double q1costh = 2 * xpts[1] - 1.;
   double q2mag = xpts[2] * qmax;
   double q2costh = 2 * xpts[3] - 1.;
   double q2phi = 2 * pi * xpts[4];

   // jacobian
   // each magnitude, cos theta gives 2 * qmax * qmag^2
   // the relative phi gives 2*pi, as does the free total phi integral
   // and a 1/(2pi)^6 from the measure
   double jacobian = q1mag*q1mag * q2mag*q2mag * qmax*qmax / (4 * pow(pi,4));

   // 3-vectors for the loop momenta
   momenta[Momentum::q] = ThreeVector(q1mag * sqrt(1. - q1costh*q1costh), 0, q1mag * q1costh);
   momenta[Momentum::q2] = ThreeVector(q2mag * sqrt(1. - q2costh*q2costh) * cos(q2phi), q2mag * sqrt(1. - q2costh*q2costh) * sin(q2phi), q2mag * q2costh);

   return std::pair<double, LabelMap<Momentum, ThreeVector>* const>(jacobian, &momenta);
}

//------------------------------------------------------------------------------
int PowerSpectrum::oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the integrator
   LoopPhaseSpace* phasespace = static_cast<LoopPhaseSpace*>(userdata);

   // generate the PS point
   std::vector<double> xpts;
   for (int i = 0; i < phasespace->ndim; i++) {
      xpts.push_back(xx[i]);
   }
   std::pair<double, LabelMap<Momentum, ThreeVector>* const> PSpoint = phasespace->generate_point_oneLoop(xpts);

   // calculate the integrand
   double integrand = PSpoint.first * (phasespace->powerspectrum->diagrams()->value_oneLoop(*(PSpoint.second), *(phasespace->kernels), phasespace->PL));

   ff[0] = integrand;

   return 0;
}

//------------------------------------------------------------------------------
int PowerSpectrum::twoLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the integrator
   LoopPhaseSpace* phasespace = static_cast<LoopPhaseSpace*>(userdata);

   // generate the PS point
   std::vector<double> xpts;
   for (int i = 0; i < phasespace->ndim; i++) {
      xpts.push_back(xx[i]);
   }
   std::pair<double, LabelMap<Momentum, ThreeVector>* const> PSpoint = phasespace->generate_point_twoLoop(xpts);

   // calculate the integrand
   double integrand = PSpoint.first * (phasespace->powerspectrum->diagrams()->value_twoLoop(*(PSpoint.second), *(phasespace->kernels), phasespace->PL));

   ff[0] = integrand;

   return 0;
}

} // namespace fnfast

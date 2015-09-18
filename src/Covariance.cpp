//------------------------------------------------------------------------------
/// \file Covariance.cpp
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
//    Implementation of class Covariance
//------------------------------------------------------------------------------

#include <iostream>

#include "Covariance.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
Covariance::Covariance(Order order)
: _order(order), _diagrams(DiagramSet4point(_order)), _UVcutoff(10.)
{}

//------------------------------------------------------------------------------
IntegralResult Covariance::tree(double k, double kprime, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // integration method
   PhaseSpace phasespace(k, kprime, _UVcutoff, &kernels, PL, this);
   phasespace.ndim = 1;

   // VEGAS integration via cuba
   VEGASintegrator vegas(phasespace.ndim);

   return vegas.integrate(tree_integrand, &phasespace);
}

//------------------------------------------------------------------------------
std::pair<double, LabelMap<Momentum, ThreeVector>* const> Covariance::PhaseSpace::generate_point_tree(std::vector<double> xpts)
{
   // angle
   double xth = 2*xpts[0] - 1;

   // jacobian = 2 from the costheta integral
   double jacobian = 2;

   // set the external momenta
   momenta[Momentum::k1] = ThreeVector(0, 0, k);
   momenta[Momentum::k2] = -momenta[Momentum::k1];
   momenta[Momentum::k3] = ThreeVector(kprime * sqrt(1. - xth*xth), 0, kprime * xth);
   momenta[Momentum::k4] = -momenta[Momentum::k3];

   return std::pair<double, LabelMap<Momentum, ThreeVector>* const>(jacobian, &momenta);
}

//------------------------------------------------------------------------------
int Covariance::tree_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the integrator
   PhaseSpace* phasespace = static_cast<PhaseSpace*>(userdata);

   // generate the PS point
   std::vector<double> xpts;
   for (int i = 0; i < phasespace->ndim; i++) {
      xpts.push_back(xx[i]);
   }
   std::pair<double, LabelMap<Momentum, ThreeVector>* const> PSpoint = phasespace->generate_point_tree(xpts);

   // calculate the integrand
   double integrand = PSpoint.first * (phasespace->covariance->diagrams()->value_tree(*(PSpoint.second), *(phasespace->kernels), phasespace->PL));

   ff[0] = integrand;

   return 0;
}

//------------------------------------------------------------------------------
IntegralResult Covariance::oneLoop(double k, double kprime, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // integration method
   PhaseSpace phasespace(k, kprime, _UVcutoff, &kernels, PL, this);
   phasespace.ndim = 4;

   // VEGAS integration via cuba
   VEGASintegrator vegas(phasespace.ndim);

   return vegas.integrate(oneLoop_integrand, &phasespace);
}

//------------------------------------------------------------------------------
std::pair<double, LabelMap<Momentum, ThreeVector>* const> Covariance::PhaseSpace::generate_point_oneLoop(std::vector<double> xpts)
{
   // we sample q flat in spherical coordinates
   // q components
   double qmag = xpts[0] * qmax;
   double qcosth = 2 * xpts[1] - 1.;
   double qphi = 2*pi * xpts[2];
   // angle
   double xth = 2*xpts[3] - 1;

   // jacobian
   // qmax from the magnitude integral,
   // 2 from the cos theta jacobian,
   // pick up a 2pi from the phi integral,
   // and a 1/(2pi)^3 from the measure
   // and a 2 from the costheta integral
   double jacobian = qmag * qmag * qmax / (pi*pi);

   // 3-vector for the loop momentum
   momenta[Momentum::q] = ThreeVector(qmag * sqrt(1. - qcosth*qcosth) * cos(qphi), qmag * sqrt(1. - qcosth*qcosth) * sin(qphi), qmag * qcosth);
   // set the external momenta
   momenta[Momentum::k1] = ThreeVector(0, 0, k);
   momenta[Momentum::k2] = -momenta[Momentum::k1];
   momenta[Momentum::k3] = ThreeVector(kprime * sqrt(1. - xth*xth), 0, kprime * xth);
   momenta[Momentum::k4] = -momenta[Momentum::k3];

   return std::pair<double, LabelMap<Momentum, ThreeVector>* const>(jacobian, &momenta);
}

//------------------------------------------------------------------------------
int Covariance::oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the integrator
   PhaseSpace* phasespace = static_cast<PhaseSpace*>(userdata);

   // generate the PS point
   std::vector<double> xpts;
   for (int i = 0; i < phasespace->ndim; i++) {
      xpts.push_back(xx[i]);
   }
   std::pair<double, LabelMap<Momentum, ThreeVector>* const> PSpoint = phasespace->generate_point_oneLoop(xpts);

   // calculate the integrand
   double integrand = PSpoint.first * (phasespace->covariance->diagrams()->value_oneLoop(*(PSpoint.second), *(phasespace->kernels), phasespace->PL));

   ff[0] = integrand;

   return 0;
}

} // namespace fnfast

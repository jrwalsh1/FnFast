//------------------------------------------------------------------------------
/// \file Bispectrum.cpp
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
//    Implementation of class Bispectrum
//------------------------------------------------------------------------------

#include <iostream>

#include "Bispectrum.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
/*DAN*/
Bispectrum::Bispectrum(Order order)
: _order(order), _diagrams(DiagramSet3pointSPT(_order)), _EFTdiagrams(DiagramSet3pointEFT(_EFTorder(_order))), _UVcutoff(10.)
{}

//------------------------------------------------------------------------------
double Bispectrum::tree(double k1, double k2, double theta12, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // set the external momenta
   ThreeVector k1vec(0, 0, k1);
   ThreeVector k2vec(k2 * sin(theta12), 0, k2 * cos(theta12));
   ThreeVector k3vec = -k1vec - k2vec;
   LabelMap<Momentum, ThreeVector> momenta {{Momentum::k1, k1vec}, {Momentum::k2, k2vec}, {Momentum::k3, k3vec}};

   return _diagrams.value_tree(momenta, kernels, PL);
}

//------------------------------------------------------------------------------
IntegralResult Bispectrum::oneLoop(double k1, double k2, double theta12, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // integration method
   LoopPhaseSpace phasespace(k1, k2, theta12, _UVcutoff, &kernels, PL, this);

   // VEGAS integration via cuba
   VEGASintegrator vegas(3);

   return vegas.integrate(oneLoop_integrand, &phasespace);
}
   
//------------------------------------------------------------------------------
/*DAN*/
double Bispectrum::treeEFT(double k1, double k2, double theta12, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // set the external momenta
   ThreeVector k1vec(0, 0, k1);
   ThreeVector k2vec(k2 * sin(theta12), 0, k2 * cos(theta12));
   ThreeVector k3vec = -k1vec - k2vec;
   LabelMap<Momentum, ThreeVector> momenta {{Momentum::k1, k1vec}, {Momentum::k2, k2vec}, {Momentum::k3, k3vec}};
      
   return _EFTdiagrams.value_tree(momenta, kernels, PL);
}

//------------------------------------------------------------------------------
std::pair<double, LabelMap<Momentum, ThreeVector>* const> Bispectrum::LoopPhaseSpace::generate_point_oneLoop(std::vector<double> xpts)
{
   // we sample q flat in spherical coordinates
   // q components
   double qmag = xpts[0] * qmax;
   double qcosth = 2 * xpts[1] - 1.;
   double qphi = 2*pi * xpts[2];

   // jacobian
   // qmax from the magnitude integral,
   // 2 from the cos theta jacobian,
   // pick up a 2pi from the phi integral,
   // and a 1/(2pi)^3 from the measure
   double jacobian = qmag * qmag * qmax / (2 * pi*pi);

   // 3-vector for the loop momentum
   momenta[Momentum::q] = ThreeVector(qmag * sqrt(1. - qcosth*qcosth) * cos(qphi), qmag * sqrt(1. - qcosth*qcosth) * sin(qphi), qmag * qcosth);

   return std::pair<double, LabelMap<Momentum, ThreeVector>* const>(jacobian, &momenta);
}

//------------------------------------------------------------------------------
int Bispectrum::oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
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
   double integrand = PSpoint.first * (phasespace->bispectrum->diagrams()->value_oneLoop(*(PSpoint.second), *(phasespace->kernels), phasespace->PL));

   ff[0] = integrand;

   return 0;
}

} // namespace fnfast

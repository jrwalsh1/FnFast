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
PowerSpectrum::PowerSpectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients) : _order(order), _PL(PL), _SPTkernels(new SPTkernels), _EFTkernels(new EFTkernels), _eftcoefficients(eftcoefficients), _UVcutoff(10.), _kBin(0), _W(NULL)
{
   // Diagram momenta
   _labels = { Momenta::k1, Momenta::k2, Momenta::q };

   // vertex kernels
   // SPT
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_SPT = {{Vertices::v1, _SPTkernels}, {Vertices::v2, _SPTkernels}};
   // counterterms
   _EFTkernels->set_coefficients(*_eftcoefficients);
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_EFTSPT = {{Vertices::v1, _EFTkernels}, {Vertices::v2, _SPTkernels}};

   // set up the diagrams, starting with the tree level
   // P11
   // propagators
   Propagator prop_P11_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
   // lines
   Line line_P11_12(Vertices::v1, Vertices::v2, prop_P11_k2);
   vector<Line> lines_P11 {line_P11_12};
   // define the diagram
   Diagram* P11 = new Diagram(lines_P11, kernels_SPT, _PL);

   // define the tree diagrams
   _tree = {P11};
   _diagrams[Graphs::P11] = P11;
    
   if (_order == kOneLoop) {
      // P31
      // propagators
      Propagator prop_P31_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_P31_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
      // lines
      Line line_P31_11(Vertices::v1, Vertices::v1, prop_P31_q);
      Line line_P31_12(Vertices::v1, Vertices::v2, prop_P31_k2);
      vector<Line> lines_P31 {line_P31_11, line_P31_12};
      // define the diagram
      Diagram* P31 = new Diagram(lines_P31, kernels_SPT, _PL);

      // P22
      // propagators
      Propagator prop_P22_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_P22_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      // lines
      Line line_P22_12a(Vertices::v1, Vertices::v2, prop_P22_q);
      Line line_P22_12b(Vertices::v1, Vertices::v2, prop_P22_qk2);
      vector<Line> lines_P22 {line_P22_12a, line_P22_12b};
      // define the diagram
      Diagram* P22 = new Diagram(lines_P22, kernels_SPT, _PL);
       
      // define the loop diagrams
      _loop = {P31, P22};
      _diagrams[Graphs::P31] = P31;
      _diagrams[Graphs::P22] = P22;
       
      // P31x
      // propagators
      Propagator prop_P31x_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
      // lines
      Line line_P31x_12(Vertices::v1, Vertices::v2, prop_P31x_k2);
      vector<Line> lines_P31x {line_P31x_12};
      // define the diagram
      Diagram* P31x = new Diagram(lines_P31x, kernels_EFTSPT, _PL);
      P31x->set_perms(P31->get_perms());

      // define the counterterm diagrams
      _cterms = {P31x};
      _diagrams[Graphs::P31x] = P31x;
       
   }
}

//------------------------------------------------------------------------------
double PowerSpectrum::tree(double k)
{
   // set the external momenta
   ThreeVector k2(0, 0, k);
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2}, {Momenta::k2, k2}});

   double value = 0;
   // sum the tree level diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double PowerSpectrum::loopSPT_excl(ThreeVector k, ThreeVector q)
{
   // set the external momenta
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k}, {Momenta::k2, k}, {Momenta::q, q}});
    
   double value = 0;
   // sum the loop diagrams
   for (size_t i = 0; i < _loop.size(); i++) {
      value += _loop[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
IntegralResult PowerSpectrum::loopSPT(double k)
{
   // options passed into the integration
   LoopIntegrationOptions data;
   data.k = k;
   data.powerspectrum = this;
   double qmax = 10;
   LoopPhaseSpace loopPS(qmax);
   data.loopPS = &loopPS;

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
   Vegas(ndim, ncomp, loop_integrand, &data, nvec,
       epsrel, epsabs, flags, vegasseed,
       mineval, maxeval, nstart, nincrease, nbatch,
       gridnum, statefile, spin,
       &neval, &fail, integral, error, prob);

   // save the results in a container
   IntegralResult result(integral[0], error[0], prob[0]);

   return result;
}

//------------------------------------------------------------------------------
double PowerSpectrum::ctermsEFT(double k)
{
   // set the external momenta
   ThreeVector k2(0, 0, k);
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2}, {Momenta::k2, k2}});

   double value = 0;
   for (size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double PowerSpectrum::LoopPhaseSpace::setPS(double qpts[2])
{
   // we sample q flat in spherical coordinates
   // q components
   double qmag = qpts[0] * _qmax;
   double qcosth = 2 * qpts[1] - 1.;

   // jacobian
   // qmax from the magnitude integral,
   // 2 from the cos theta jacobian,
   // pick up a 2pi from the phi integral,
   // and a 1/(2pi)^3 from the measure
   _jacobian = qmag * qmag * _qmax / (2 * pi*pi);

   // 3-vector for the loop momentum
   _q = ThreeVector(qmag * sqrt(1. - qcosth*qcosth), 0, qmag * qcosth);

   return _jacobian;
}

//------------------------------------------------------------------------------
int PowerSpectrum::loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the options
   LoopIntegrationOptions* data = static_cast<LoopIntegrationOptions*>(userdata);

   // external momentum magnitude
   double k = data->k;

   // define the variables needed for the PS point
   double qpts[2] = {xx[0], xx[1]};

   // set the PS point and return the integrand
   double jacobian = data->loopPS->setPS(qpts);
   double integrand = 0;
   if (jacobian > 0) {
      ThreeVector q = data->loopPS->q();
      ThreeVector k2(0, 0, k);
      integrand = data->powerspectrum->loopSPT_excl(k2, q);
   }

   // loop calculation
   ff[0] = jacobian * integrand;

   return 0;
}

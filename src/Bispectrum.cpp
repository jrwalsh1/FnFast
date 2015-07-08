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
#include "cuba.h"

//------------------------------------------------------------------------------
Bispectrum::Bispectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients) : _order(order), _PL(PL), _SPTkernels(new SPTkernels), _EFTkernels(new EFTkernels), _eftcoefficients(eftcoefficients), _UVcutoff(10.), _kBin(0), _W(NULL)
{
   // Diagram momenta
   _labels = { Momenta::k1, Momenta::k2, Momenta::k3, Momenta::q };

   // vertex kernels
   // SPT
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_SPT = {{Vertices::v1, _SPTkernels}, {Vertices::v2, _SPTkernels}, {Vertices::v3, _SPTkernels}};
   // counterterms
   _EFTkernels->set_coefficients(*_eftcoefficients);
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_EFTSPT = {{Vertices::v1, _EFTkernels}, {Vertices::v2, _SPTkernels}, {Vertices::v3, _SPTkernels}};

   // set up the diagrams, starting with the tree level
   // B211
   // propagators
   Propagator prop_B211_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
   Propagator prop_B211_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
   // lines
   Line line_B211_12(Vertices::v1, Vertices::v2, prop_B211_k2);
   Line line_B211_13(Vertices::v1, Vertices::v3, prop_B211_k3);
   vector<Line> lines_B211 {line_B211_12, line_B211_13};
   // define the diagram
   Diagram* B211 = new Diagram(lines_B211, kernels_SPT, _PL);

   // define the tree diagrams
   _tree = {B211};
   _diagrams[Graphs::B211] = B211;

   if (_order == kOneLoop) {
      // B411
      // propagators
      Propagator prop_B411_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B411_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
      Propagator prop_B411_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B411_11(Vertices::v1, Vertices::v1, prop_B411_q);
      Line line_B411_12(Vertices::v1, Vertices::v2, prop_B411_k2);
      Line line_B411_13(Vertices::v1, Vertices::v3, prop_B411_k3);
      vector<Line> lines_B411 {line_B411_11, line_B411_12, line_B411_13};
      // define the diagram
      Diagram* B411 = new Diagram(lines_B411, kernels_SPT, _PL);

      // B321a
      // propagators
      Propagator prop_B321a_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B321a_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_B321a_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B321a_11(Vertices::v1, Vertices::v1, prop_B321a_q);
      Line line_B321a_12(Vertices::v1, Vertices::v2, prop_B321a_k23);
      Line line_B321a_23(Vertices::v2, Vertices::v3, prop_B321a_k3);
      vector<Line> lines_B321a {line_B321a_11, line_B321a_12, line_B321a_23};
      // define the diagram
      Diagram* B321a = new Diagram(lines_B321a, kernels_SPT, _PL);

      // B321b
      // propagators
      Propagator prop_B321b_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B321b_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      Propagator prop_B321b_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B321b_12a(Vertices::v1, Vertices::v2, prop_B321b_q);
      Line line_B321b_12b(Vertices::v1, Vertices::v2, prop_B321b_qk2);
      Line line_B321b_13(Vertices::v1, Vertices::v3, prop_B321b_k3);
      vector<Line> lines_B321b {line_B321b_12a, line_B321b_12b, line_B321b_13};
      // define the diagram
      Diagram* B321b = new Diagram(lines_B321b, kernels_SPT, _PL);

      // B222
      // propagators
      Propagator prop_B222_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B222_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}, {Momenta::k2, -1}});
      Propagator prop_B222_qk23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}, {Momenta::k3, 1}});
      // lines
      Line line_B222_12(Vertices::v1, Vertices::v2, prop_B222_q);
      Line line_B222_23(Vertices::v2, Vertices::v3, prop_B222_qk2);
      Line line_B222_13(Vertices::v1, Vertices::v3, prop_B222_qk23);
      vector<Line> lines_B222 {line_B222_12, line_B222_23, line_B222_13};
      // define the diagram
      Diagram* B222 = new Diagram(lines_B222, kernels_SPT, _PL);
       
      // define the loop diagrams
      _loop = {B411, B321a, B321b, B222};
      _diagrams[Graphs::B411] = B411;
      _diagrams[Graphs::B321a] = B321a;
      _diagrams[Graphs::B321b] = B321b;
      _diagrams[Graphs::B222] = B222;
       
      // B411x
      // propagators
      Propagator prop_B411x_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
      Propagator prop_B411x_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B411x_12(Vertices::v1, Vertices::v2, prop_B411x_k2);
      Line line_B411x_13(Vertices::v1, Vertices::v3, prop_B411x_k3);
      vector<Line> lines_B411x {line_B411x_12, line_B411x_13};
      // define the diagram
      Diagram* B411x = new Diagram(lines_B411x, kernels_EFTSPT, _PL);
      B411x->set_perms(B411->get_perms());
       
      // B321bx
      // propagators
      Propagator prop_B321ax_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_B321ax_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B321ax_12(Vertices::v1, Vertices::v2, prop_B321ax_k23);
      Line line_B321ax_23(Vertices::v2, Vertices::v3, prop_B321ax_k3);
      vector<Line> lines_B321ax {line_B321ax_12, line_B321ax_23};
      // define the diagram
      Diagram* B321ax = new Diagram(lines_B321ax, kernels_EFTSPT, _PL);
      B321ax->set_perms(B321b->get_perms());
       
      // define the counterterm diagrams
      _cterms = {B411x,B321ax};
      _diagrams[Graphs::B411x] = B411x;
      _diagrams[Graphs::B321ax] = B321ax;
   }
}

//------------------------------------------------------------------------------
double Bispectrum::tree(double k1, double k2, double costheta12)
{
   // set the external momenta
   ThreeVector pk1(0, 0, k1);
   ThreeVector pk2(k2 * sqrt(1. - costheta12*costheta12), 0, k2 * costheta12);
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, pk1}, {Momenta::k2, pk2}, {Momenta::k3, -pk1-pk2}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Bispectrum::loopSPT_excl(ThreeVector k1, ThreeVector k2, ThreeVector q)
{
   // set the external momenta
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, -k1-k2}, {Momenta::q, q}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _loop.size(); i++) {
      value += _loop[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Bispectrum::loopSPT(double k1, double k2, double costheta12)
{
   // options passed into the integration
   LoopIntegrationOptions data;
   data.k1 = k1;
   data.k2 = k2;
   data.costheta12 = costheta12;
   data.bispectrum = this;
   double qmax = 10;
   LoopPhaseSpace loopPS(qmax);
   data.loopPS = &loopPS;

   // Integration
   // VEGAS parameters
   // PS dimensionality
   // q: 3-dim
   const int ndim = 3;
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // absolute uncertainty (safeguard)
   const double epsrel = 1e-4;
   const double epsabs = 1e-12;
   // min, max number of points
   const int mineval = 0;
   const int maxeval = 100000;
//   const int mineval = 10;
//   const int maxeval = 100000;
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
   // cubature rule degree
//   int key = 13;
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
   int nregions, neval, fail;

   // containers for output
   double integral[ncomp], error[ncomp], prob[ncomp];

   // run VEGAS
   Vegas(ndim, ncomp, loop_integrand, &data, nvec,
       epsrel, epsabs, flags, vegasseed,
       mineval, maxeval, nstart, nincrease, nbatch,
       gridnum, statefile, spin,
       &neval, &fail, integral, error, prob);

   /*
   int nnew = 2000;
   int nmin = 100;
   double flatness = 5;

   Suave(ndim, ncomp, loop_integrand, &data, nvec,
      epsrel, epsabs, flags, vegasseed,
      mineval, maxeval, nnew, nmin, flatness,
      statefile, spin, &nregions,
      &neval, &fail, integral, error, prob);
   */
   /*
   // run Cuhre
   Cuhre(ndim, ncomp, loop_integrand, &data, nvec,
       epsrel, epsabs, flags,
       mineval, maxeval,
       key, statefile, spin, &nregions,
       &neval, &fail, integral, error, prob);
   */

   cout << "integral, error, prob = " << integral[0] << ", " << error[0] << ", " << prob[0] << endl;

   return integral[0];
}

//------------------------------------------------------------------------------
double Bispectrum::ctermsEFT(double k1, double k2, double costheta12)
{
   // set the external momenta
   ThreeVector pk1(0, 0, k1);
   ThreeVector pk2(k2 * sqrt(1. - costheta12*costheta12), 0, k2 * costheta12);
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, pk1}, {Momenta::k2, pk2}, {Momenta::k3, -pk1-pk2}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Bispectrum::LoopPhaseSpace::setPS(double qpts[3])
{
   // we sample q flat in spherical coordinates, setting k1 along the z-axis, k2 in the x-z plane
   // q components
   double qmag = qpts[0] * _qmax;
   double qcosth = 2 * qpts[1] - 1.;
   double qphi = 2*pi * qpts[2];

   // jacobian
   // qmax from the magnitude integral,
   // 2 from the cos theta jacobian,
   // pick up a 2pi from the phi integral,
   // and a 1/(2pi)^3 from the measure
   _jacobian = qmag * qmag * _qmax / (2 * pi*pi);

   // 3-vector for the loop momentum
   _q = ThreeVector(qmag * sqrt(1. - qcosth*qcosth) * cos(qphi), qmag * sqrt(1. - qcosth*qcosth) * sin(qphi), qmag * qcosth);

   return _jacobian;
}

//------------------------------------------------------------------------------
int Bispectrum::loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the options
   LoopIntegrationOptions* data = static_cast<LoopIntegrationOptions*>(userdata);

   // external momentum magnitude
   double k1 = data->k1;
   double k2 = data->k2;
   double costheta12 = data->costheta12;

   // define the variables needed for the PS point
   double qpts[3] = {xx[0], xx[1], xx[2]};

   // set the PS point and return the integrand
   double jacobian = data->loopPS->setPS(qpts);
   double integrand = 0;
   if (jacobian > 0) {
      ThreeVector q = data->loopPS->q();
      ThreeVector pk1(0, 0, k1);
      ThreeVector pk2(k2 * sqrt(1. - costheta12*costheta12), 0, k2 * costheta12);
      integrand = data->bispectrum->loopSPT_excl(pk1, pk2, q);
   }

   // loop calculation
   ff[0] = jacobian * integrand;

   return 0;
}

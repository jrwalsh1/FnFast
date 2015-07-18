//------------------------------------------------------------------------------
/// \file Trispectrum.cpp
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
//    Implementation of class Trispectrum
//------------------------------------------------------------------------------

#include <iostream>

#include "Trispectrum.hpp"
#include "cuba.h"

//------------------------------------------------------------------------------
Trispectrum::Trispectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients) : _order(order), _PL(PL), _SPTkernels(new SPTkernels), _EFTkernels(new EFTkernels), _eftcoefficients(eftcoefficients), _UVcutoff(10.), _kBin(0), _W(NULL), _seed(37)
{
   // Diagram momenta
   _labels = { Momenta::k1, Momenta::k2, Momenta::k3, Momenta::k4, Momenta::q };

   // vertex kernels
   // SPT
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_SPT = {{Vertices::v1, _SPTkernels}, {Vertices::v2, _SPTkernels}, {Vertices::v3, _SPTkernels}, {Vertices::v4, _SPTkernels}};
   // counterterms
   _EFTkernels->set_coefficients(*_eftcoefficients);
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_EFTSPT = {{Vertices::v1, _EFTkernels}, {Vertices::v2, _SPTkernels}, {Vertices::v3, _SPTkernels}, {Vertices::v4, _SPTkernels}};

   // set up the diagrams, starting with the tree level
   // T3111
   // propagators
   Propagator prop_T3111_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
   Propagator prop_T3111_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
   Propagator prop_T3111_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
   // lines
   Line line_T3111_12(Vertices::v1, Vertices::v2, prop_T3111_k2);
   Line line_T3111_13(Vertices::v1, Vertices::v3, prop_T3111_k3);
   Line line_T3111_14(Vertices::v1, Vertices::v4, prop_T3111_k4);
   vector<Line> lines_T3111 {line_T3111_12, line_T3111_13, line_T3111_14};
   // define the diagram
   Diagram* T3111 = new Diagram(lines_T3111, kernels_SPT, _PL);

   // T2211
   // propagators
   Propagator prop_T2211_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
   Propagator prop_T2211_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
   Propagator prop_T2211_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
   // lines
   Line line_T2211_12(Vertices::v1, Vertices::v2, prop_T2211_k23);
   Line line_T2211_23(Vertices::v2, Vertices::v3, prop_T2211_k3);
   Line line_T2211_14(Vertices::v1, Vertices::v4, prop_T2211_k4);
   vector<Line> lines_T2211 {line_T2211_12, line_T2211_23, line_T2211_14};
   // define the diagram
   Diagram* T2211 = new Diagram(lines_T2211, kernels_SPT, _PL);

   // define the tree diagrams
   _tree = {T3111, T2211};
   _diagrams[Graphs::T3111] = T3111;
   _diagrams[Graphs::T2211] = T2211;

   if (_order == kOneLoop) {
      // T5111
      // propagators
      Propagator prop_T5111_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T5111_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
      Propagator prop_T5111_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T5111_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T5111_11(Vertices::v1, Vertices::v1, prop_T5111_q);
      Line line_T5111_12(Vertices::v1, Vertices::v2, prop_T5111_k2);
      Line line_T5111_13(Vertices::v1, Vertices::v3, prop_T5111_k3);
      Line line_T5111_14(Vertices::v1, Vertices::v4, prop_T5111_k4);
      vector<Line> lines_T5111 {line_T5111_11, line_T5111_12, line_T5111_13, line_T5111_14};
      // define the diagram
      Diagram* T5111 = new Diagram(lines_T5111, kernels_SPT, _PL);

      // T4211a
      // propagators
      Propagator prop_T4211a_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T4211a_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_T4211a_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T4211a_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T4211a_11(Vertices::v1, Vertices::v1, prop_T4211a_q);
      Line line_T4211a_12(Vertices::v1, Vertices::v2, prop_T4211a_k23);
      Line line_T4211a_23(Vertices::v2, Vertices::v3, prop_T4211a_k3);
      Line line_T4211a_14(Vertices::v1, Vertices::v4, prop_T4211a_k4);
      vector<Line> lines_T4211a {line_T4211a_11, line_T4211a_12, line_T4211a_23, line_T4211a_14};
      // define the diagram
      Diagram* T4211a = new Diagram(lines_T4211a, kernels_SPT, _PL);

      // T4211b
      // propagators
      Propagator prop_T4211b_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T4211b_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      Propagator prop_T4211b_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T4211b_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T4211b_12a(Vertices::v1, Vertices::v2, prop_T4211b_q);
      Line line_T4211b_12b(Vertices::v1, Vertices::v2, prop_T4211b_qk2);
      Line line_T4211b_13(Vertices::v1, Vertices::v3, prop_T4211b_k3);
      Line line_T4211b_14(Vertices::v1, Vertices::v4, prop_T4211b_k4);
      vector<Line> lines_T4211b {line_T4211b_12a, line_T4211b_12b, line_T4211b_13, line_T4211b_14};
      // define the diagram
      Diagram* T4211b = new Diagram(lines_T4211b, kernels_SPT, _PL);

      // T3311a
      // propagators
      Propagator prop_T3311a_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T3311a_k234(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3311a_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T3311a_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3311a_11(Vertices::v1, Vertices::v1, prop_T3311a_q);
      Line line_T3311a_12(Vertices::v1, Vertices::v2, prop_T3311a_k234);
      Line line_T3311a_23(Vertices::v2, Vertices::v3, prop_T3311a_k3);
      Line line_T3311a_24(Vertices::v2, Vertices::v4, prop_T3311a_k4);
      vector<Line> lines_T3311a {line_T3311a_11, line_T3311a_12, line_T3311a_23, line_T3311a_24};
      // define the diagram
      Diagram* T3311a = new Diagram(lines_T3311a, kernels_SPT, _PL);

      // T3311b
      // propagators
      Propagator prop_T3311b_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T3311b_qk23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_T3311b_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T3311b_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3311b_12a(Vertices::v1, Vertices::v2, prop_T3311b_q);
      Line line_T3311b_12b(Vertices::v1, Vertices::v2, prop_T3311b_qk23);
      Line line_T3311b_23(Vertices::v2, Vertices::v3, prop_T3311b_k3);
      Line line_T3311b_14(Vertices::v1, Vertices::v4, prop_T3311b_k4);
      vector<Line> lines_T3311b {line_T3311b_12a, line_T3311b_12b, line_T3311b_23, line_T3311b_14};
      // define the diagram
      Diagram* T3311b = new Diagram(lines_T3311b, kernels_SPT, _PL);

      // T3221a
      // propagators
      Propagator prop_T3221a_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T3221a_k234(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3221a_k34(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3221a_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3221a_11(Vertices::v1, Vertices::v1, prop_T3221a_q);
      Line line_T3221a_12(Vertices::v1, Vertices::v2, prop_T3221a_k234);
      Line line_T3221a_23(Vertices::v2, Vertices::v3, prop_T3221a_k34);
      Line line_T3221a_34(Vertices::v3, Vertices::v4, prop_T3221a_k4);
      vector<Line> lines_T3221a {line_T3221a_11, line_T3221a_12, line_T3221a_23, line_T3221a_34};
      // define the diagram
      Diagram* T3221a = new Diagram(lines_T3221a, kernels_SPT, _PL);

      // T3221b
      // propagators
      Propagator prop_T3221b_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T3221b_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      Propagator prop_T3221b_k34(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3221b_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3221b_12a(Vertices::v1, Vertices::v2, prop_T3221b_q);
      Line line_T3221b_12b(Vertices::v1, Vertices::v2, prop_T3221b_qk2);
      Line line_T3221b_13(Vertices::v1, Vertices::v3, prop_T3221b_k34);
      Line line_T3221b_34(Vertices::v3, Vertices::v4, prop_T3221b_k4);
      vector<Line> lines_T3221b {line_T3221b_12a, line_T3221b_12b, line_T3221b_13, line_T3221b_34};
      // define the diagram
      Diagram* T3221b = new Diagram(lines_T3221b, kernels_SPT, _PL);

      // T3221c
      // propagators
      Propagator prop_T3221c_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T3221c_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}, {Momenta::k2, -1}});
      Propagator prop_T3221c_qk23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_T3221c_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3221c_12(Vertices::v1, Vertices::v2, prop_T3221c_q);
      Line line_T3221c_23(Vertices::v2, Vertices::v3, prop_T3221c_qk2);
      Line line_T3221c_13(Vertices::v1, Vertices::v3, prop_T3221c_qk23);
      Line line_T3221c_14(Vertices::v1, Vertices::v4, prop_T3221c_k4);
      vector<Line> lines_T3221c {line_T3221c_12, line_T3221c_23, line_T3221c_13, line_T3221c_14};
      // define the diagram
      Diagram* T3221c = new Diagram(lines_T3221c, kernels_SPT, _PL);

      // T2222
      // propagators
      Propagator prop_T2222_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_T2222_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}, {Momenta::k2, -1}});
      Propagator prop_T2222_qk23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}, {Momenta::k2, -1}, {Momenta::k3, -1}});
      Propagator prop_T2222_qk234(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}, {Momenta::k3, 1}, {Momenta::k4, 1}});
      // lines
      Line line_T2222_12(Vertices::v1, Vertices::v2, prop_T2222_q);
      Line line_T2222_23(Vertices::v2, Vertices::v3, prop_T2222_qk2);
      Line line_T2222_34(Vertices::v3, Vertices::v4, prop_T2222_qk23);
      Line line_T2222_14(Vertices::v1, Vertices::v4, prop_T2222_qk234);
      vector<Line> lines_T2222 {line_T2222_12, line_T2222_23, line_T2222_34, line_T2222_14};
      // define the diagram
      Diagram* T2222 = new Diagram(lines_T2222, kernels_SPT, _PL);

      // define the loop diagrams
      _loop = {T5111, T4211a, T4211b, T3311a, T3311b, T3221a, T3221b, T3221c, T2222};
      _diagrams[Graphs::T5111] = T5111;
      _diagrams[Graphs::T4211a] = T4211a;
      _diagrams[Graphs::T4211b] = T4211b;
      _diagrams[Graphs::T3311a] = T3311a;
      _diagrams[Graphs::T3311b] = T3311b;
      _diagrams[Graphs::T3221a] = T3221a;
      _diagrams[Graphs::T3221b] = T3221b;
      _diagrams[Graphs::T3221c] = T3221c;
      _diagrams[Graphs::T2222] = T2222;

      // T5111x
      // propagators
      Propagator prop_T5111x_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
      Propagator prop_T5111x_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T5111x_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T5111x_12(Vertices::v1, Vertices::v2, prop_T5111x_k2);
      Line line_T5111x_13(Vertices::v1, Vertices::v3, prop_T5111x_k3);
      Line line_T5111x_14(Vertices::v1, Vertices::v4, prop_T5111x_k4);
      vector<Line> lines_T5111x {line_T5111x_12, line_T5111x_13, line_T5111x_14};
      // define the diagram
      Diagram* T5111x = new Diagram(lines_T5111x, kernels_EFTSPT, _PL);
      T5111x->set_perms(T5111->get_perms());
       
      // T4211ax
      // propagators
      Propagator prop_T4211ax_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_T4211ax_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T4211ax_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T4211ax_12(Vertices::v1, Vertices::v2, prop_T4211ax_k23);
      Line line_T4211ax_23(Vertices::v2, Vertices::v3, prop_T4211ax_k3);
      Line line_T4211ax_14(Vertices::v1, Vertices::v4, prop_T4211ax_k4);
      vector<Line> lines_T4211ax {line_T4211ax_12, line_T4211ax_23, line_T4211ax_14};
      // define the diagram
      Diagram* T4211ax = new Diagram(lines_T4211ax, kernels_EFTSPT, _PL);
      T4211ax->set_perms(T4211a->get_perms());

      // T3311ax
      // propagators
      Propagator prop_T3311ax_k234(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3311ax_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      Propagator prop_T3311ax_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3311ax_12(Vertices::v1, Vertices::v2, prop_T3311ax_k234);
      Line line_T3311ax_23(Vertices::v2, Vertices::v3, prop_T3311ax_k3);
      Line line_T3311ax_24(Vertices::v2, Vertices::v4, prop_T3311ax_k4);
      vector<Line> lines_T3311ax {line_T3311ax_12, line_T3311ax_23, line_T3311ax_24};
      // define the diagram
      Diagram* T3311ax = new Diagram(lines_T3311ax, kernels_EFTSPT, _PL);
      T3311ax->set_perms(T3311a->get_perms());

      // T3221ax
      // propagators
      Propagator prop_T3221ax_k234(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3221ax_k34(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}, {Momenta::k4, 1}});
      Propagator prop_T3221ax_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
      // lines
      Line line_T3221ax_12(Vertices::v1, Vertices::v2, prop_T3221ax_k234);
      Line line_T3221ax_23(Vertices::v2, Vertices::v3, prop_T3221ax_k34);
      Line line_T3221ax_34(Vertices::v3, Vertices::v4, prop_T3221ax_k4);
      vector<Line> lines_T3221ax {line_T3221ax_12, line_T3221ax_23, line_T3221ax_34};
      // define the diagram
      Diagram* T3221ax = new Diagram(lines_T3221ax, kernels_EFTSPT, _PL);
      T3221ax->set_perms(T3221a->get_perms());

      // define the counterterm diagrams
      _cterms = {T5111x,T4211ax,T3311ax,T3221ax};
      _diagrams[Graphs::T5111x] = T5111x;
      _diagrams[Graphs::T4211ax] = T4211ax;
      _diagrams[Graphs::T3311ax] = T3311ax;
      _diagrams[Graphs::T3221ax] = T3221ax;
   }
}

//******************************************************************************
// tree level:
// - full trispectrum
// - covariance limit, differential in angle
// - covariance limit, integrated over angle
//******************************************************************************

//------------------------------------------------------------------------------
double Trispectrum::tree(ThreeVector k1, ThreeVector k2, ThreeVector k3)
{
   // set the external momenta
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, -k1-k2-k3}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Trispectrum::cov_tree(double k, double kp, double costheta)
{
   // set the external momenta
   ThreeVector k1(0, 0, k);
   ThreeVector k2(0, 0, -k);
   ThreeVector k3(kp * sqrt(1 - costheta*costheta), 0, kp * costheta);
   ThreeVector k4(-kp * sqrt(1 - costheta*costheta), 0, -kp * costheta);
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, k4}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
IntegralResult Trispectrum::cov_tree(double k, double kp)
{
   // options passed into the integration
   AngularIntegrationOptions data;
   data.k = k;
   data.kp = kp;
   data.trispectrum = this;

   // Integration
   // CUBA parameters
   // PS dimensionality
   // theta: 1-dim
   const int ndim = 1;
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // absolute uncertainty (safeguard)
   const double epsrel = 1e-4;
   const double epsabs = 0;
   // min, max number of points
   const int mineval = 0;
   const int maxeval = 100000;
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
   const int vegasseed = _seed;
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
   Vegas(ndim, ncomp, tree_angular_integrand, &data, nvec,
       epsrel, epsabs, flags, vegasseed,
       mineval, maxeval, nstart, nincrease, nbatch,
       gridnum, statefile, spin,
       &neval, &fail, integral, error, prob);

   // save the results in a container
   IntegralResult result(integral[0], error[0], prob[0]);

   return result;
}

//------------------------------------------------------------------------------
int Trispectrum::tree_angular_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the options
   AngularIntegrationOptions* data = static_cast<AngularIntegrationOptions*>(userdata);

   // external momenta magnitudes
   double k = data->k;
   double kp = data->kp;

   // define the variables needed for the PS point
   double costheta = 2 * xx[0] - 1.;
   double jacobian = 2;

   // set the PS point and return the integrand
   double integrand = data->trispectrum->cov_tree(k, kp, costheta);

   // loop calculation
   ff[0] = jacobian * integrand;

   return 0;
}

//******************************************************************************
// loop level:
// - full trispectrum, SPT loop fully differential in loop momentum
// - full trispectrum, EFT counterterms fully differential in loop momentum
// - covariance limit, SPT loop integrated over loop momentum, but differential in angle
// - covariance limit, EFT counterterms differential in angle
// - covariance limit, SPT loop + EFT counterterms integrated over angle, loop momentum (SPT)
//******************************************************************************

//------------------------------------------------------------------------------
double Trispectrum::loopSPT_excl(ThreeVector& k1, ThreeVector& k2, ThreeVector& k3, ThreeVector& q)
{
   // set the external momenta
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, -k1-k2-k3}, {Momenta::q, q}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _loop.size(); i++) {
      value += _loop[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Trispectrum::ctermsEFT(ThreeVector k1, ThreeVector k2, ThreeVector k3)
{
   // set the external momenta
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, -k1-k2-k3}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
IntegralResult Trispectrum::cov_loopSPT(double k, double kp, double costheta)
{
   // options passed into the integration
   LoopIntegrationOptions data;
   data.k = k;
   data.kp = kp;
   data.costheta = costheta;
   data.trispectrum = this;
   double qmax = 10;
   LoopPhaseSpace loopPS(qmax);
   data.loopPS = &loopPS;

   // Integration
   // CUBA parameters
   // PS dimensionality
   // q: 3-dim
   const int ndim = 3;
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // absolute uncertainty (safeguard)
   const double epsrel = 1e-4;
   const double epsabs = 0;
   // min, max number of points
   const int mineval = 0;
   const int maxeval = 25000;
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
   const int vegasseed = _seed;
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
double Trispectrum::LoopPhaseSpace::set_loopPS(double qpts[3], double x12)
{
   // ----- COMMENT/UNCOMMENT TO SWITCH TO NON-ORTHOGONAL COORDINATES FOR PHASE SPACE SAMPLING -----
   // we sample q in non-orthogonal coordinates centered on the two principal directions
   // sample the coordinates
   double qmag = qpts[0] * _qmax;
   double s = qpts[1];
   double alpha = 2*pi * qpts[2];
   // convert to spherical coordinates to define q
   double theta12 = acos(x12);
   double qcosth = sqrt(abs(1 - s*s)) * cos(0.5*theta12 - alpha);
   double qphi = acos((1. / sqrt(abs(1 - qcosth*qcosth))) * sqrt(abs(1 - s*s)) * sin(0.5*theta12 - alpha));
   // ----- END OF BLOCK FOR NON-ORTHOGONAL COORDINATES -----


   // ----- COMMENT/UNCOMMENT TO SWITCH TO SPHERICAL COORDINATES FOR PHASE SPACE SAMPLING -----
   /*
   // we sample q flat in spherical coordinates, setting k along the z-axis, kp in the x-z plane
   // q components
   double qmag = qpts[0] * _qmax;
   double qcosth = 2 * qpts[1] - 1.;
   double qphi = 2*pi * qpts[2];
   */
   // ----- END OF BLOCK FOR SPHERICAL COORDINATES -----

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
int Trispectrum::loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the options
   LoopIntegrationOptions* data = static_cast<LoopIntegrationOptions*>(userdata);

   // external momentum magnitude
   double k = data->k;
   double kp = data->kp;
   double costheta = data->costheta;

   // define the variables needed for the PS point
   double qpts[3] = {xx[0], xx[1], xx[2]};

   // set the PS point and return the integrand
   double jacobian = data->loopPS->set_loopPS(qpts, costheta);
   double integrand = 0;
   if (jacobian > 0) {
      ThreeVector q = data->loopPS->q();
      ThreeVector k1(0, 0, k);
      ThreeVector k2(0, 0, -k);
      ThreeVector k3(kp * sqrt(1 - costheta*costheta), 0, kp * costheta);
      integrand = data->trispectrum->loopSPT_excl(k1, k2, k3, q);
   }

   // loop calculation
   ff[0] = jacobian * integrand;

   return 0;
}

//------------------------------------------------------------------------------
double Trispectrum::cov_ctermsEFT(double k, double kp, double costheta)
{
   // set the external momenta
   ThreeVector k1(0, 0, k);
   ThreeVector k2(0, 0, -k);
   ThreeVector k3(kp * sqrt(1 - costheta*costheta), 0, kp * costheta);
   ThreeVector k4(-kp * sqrt(1 - costheta*costheta), 0, -kp * costheta);
   DiagramMomenta momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, k4}});
    
   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
IntegralResult Trispectrum::cov_loopSPT(double k, double kp)
{
   // options passed into the integration
   LoopIntegrationOptions data;
   data.k = k;
   data.kp = kp;
   data.trispectrum = this;
   double qmax = 10;
   LoopPhaseSpace loopPS(qmax);
   data.loopPS = &loopPS;

   // Integration
   // CUBA parameters
   // PS dimensionality
   // q + costheta: 4-dim
   const int ndim = 4;
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // absolute uncertainty (safeguard)
   const double epsrel = 1e-4;
   const double epsabs = 0;
   // min, max number of points
   const int mineval = 0;
   const int maxeval = 100000;
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
   const int vegasseed = _seed;
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
   Vegas(ndim, ncomp, angular_loop_integrand, &data, nvec,
       epsrel, epsabs, flags, vegasseed,
       mineval, maxeval, nstart, nincrease, nbatch,
       gridnum, statefile, spin,
       &neval, &fail, integral, error, prob);

   // save the results in a container
   IntegralResult result(integral[0], error[0], prob[0]);

   return result;
}

//------------------------------------------------------------------------------
int Trispectrum::angular_loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the options
   LoopIntegrationOptions* data = static_cast<LoopIntegrationOptions*>(userdata);

   // external momentum magnitude
   double k = data->k;
   double kp = data->kp;

   // define the variables needed for the PS point
   double qpts[3] = {xx[0], xx[1], xx[2]};
   double costheta = 2 * xx[3] - 1.;

   // set the PS point and return the integrand
   double jacobian_costh = 2;
   double jacobian_loop = data->loopPS->set_loopPS(qpts, costheta);
   double integrand = 0;
   if (jacobian_loop > 0) {
      ThreeVector q = data->loopPS->q();
      ThreeVector k1(0, 0, k);
      ThreeVector k2(0, 0, -k);
      ThreeVector k3(kp * sqrt(1 - costheta*costheta), 0, kp * costheta);
      integrand = data->trispectrum->loopSPT_excl(k1, k2, k3, q);
   }

   // loop calculation
   ff[0] = jacobian_costh * jacobian_loop * integrand;

   return 0;
}

//------------------------------------------------------------------------------
IntegralResult Trispectrum::cov_ctermsEFT(double k, double kp){
   
   // options passed into the integration
   AngularIntegrationOptions data;
   data.k = k;
   data.kp = kp;
   data.trispectrum = this;
   
   // Integration
   // CUBA parameters
   // PS dimensionality
   // theta: 1-dim
   const int ndim = 1;
   // number of computations
   const int ncomp = 1; // only 1 computation
   // number of points sent to the integrand per invocation
   const int nvec = 1; // no vectorization
   // absolute uncertainty (safeguard)
   const double epsrel = 1e-4;
   const double epsabs = 0;
   // min, max number of points
   const int mineval = 0;
   const int maxeval = 100000;
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
   const int vegasseed = _seed;
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
   Vegas(ndim, ncomp, cterms_angular_integrand, &data, nvec,
         epsrel, epsabs, flags, vegasseed,
         mineval, maxeval, nstart, nincrease, nbatch,
         gridnum, statefile, spin,
         &neval, &fail, integral, error, prob);

   // save the results in a container
   IntegralResult result(integral[0], error[0], prob[0]);

   return result;
}

//------------------------------------------------------------------------------
int Trispectrum::cterms_angular_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
   // get the options
   AngularIntegrationOptions* data = static_cast<AngularIntegrationOptions*>(userdata);
   
   // external momenta magnitudes
   double k = data->k;
   double kp = data->kp;
   
   // define the variables needed for the PS point
   double costheta = 2 * xx[0] - 1.;
   double jacobian = 2;
   
   // set the PS point and return the integrand
   double integrand = data->trispectrum->cov_ctermsEFT(k, kp, costheta);
   
   // loop calculation
   ff[0] = jacobian * integrand;
   
   return 0;
}
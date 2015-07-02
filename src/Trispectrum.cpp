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

//------------------------------------------------------------------------------
Trispectrum::Trispectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients) : _order(order), _PL(PL), _SPTkernels(new SPTkernels), _EFTkernels(new EFTkernels), _eftcoefficients(eftcoefficients), _UVcutoff(10.), _kBin(0), _W(NULL)
{
   // Diagram momenta
   _labels = { Momenta::k1, Momenta::k2, Momenta::k3, Momenta::k4, Momenta::q };
   _momenta=DiagramMomenta(_labels);
    
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
      Propagator prop_T3221c_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
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
      Propagator prop_T2222_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      Propagator prop_T2222_qk23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}, {Momenta::k3, 1}});
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

//------------------------------------------------------------------------------
double Trispectrum::treeLevel_value(ThreeVector k2, ThreeVector k3, ThreeVector k4)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3-k4}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, k4}, {Momenta::q, ThreeVector(1,0,0)}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(_momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Trispectrum::oneLoopSPT_value(ThreeVector k2, ThreeVector k3, ThreeVector k4, ThreeVector q)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3-k4}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, k4}, {Momenta::q, q}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _loop.size(); i++) {
      value += _loop[i]->value_IRreg(_momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Trispectrum::oneLoopSPT_value(ThreeVector k2, ThreeVector k3, ThreeVector k4)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3-k4}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, k4}, {Momenta::q, ThreeVector(1,0,0)}});
    
   double value = 0;
   // Integration

   return value;
}

//------------------------------------------------------------------------------
double Trispectrum::oneLoopCterms_value(ThreeVector k2, ThreeVector k3, ThreeVector k4)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3-k4}, {Momenta::k2, k2}, {Momenta::k3, k3},  {Momenta::k4, k4}, {Momenta::q, ThreeVector(1,0,0)}});
    
   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(_momenta);
   }
   return value;
}
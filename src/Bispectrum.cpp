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

//------------------------------------------------------------------------------
Bispectrum::Bispectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients) : _order(order), _PL(PL), _SPTkernels(new SPTkernels), _EFTkernels(new EFTkernels), _eftcoefficients(eftcoefficients), _UVcutoff(10.), _kBin(0), _W(NULL)
{
   // Diagram momenta
   _labels = { Momenta::k1, Momenta::k2, Momenta::k3, Momenta::q };
   _momenta=DiagramMomenta(_labels);

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
      Propagator prop_B321a_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      Propagator prop_B321a_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B321a_12a(Vertices::v1, Vertices::v2, prop_B321a_q);
      Line line_B321a_12b(Vertices::v1, Vertices::v2, prop_B321a_qk2);
      Line line_B321a_13(Vertices::v1, Vertices::v3, prop_B321a_k3);
      vector<Line> lines_B321a {line_B321a_12a, line_B321a_12b, line_B321a_13};
      // define the diagram
      Diagram* B321a = new Diagram(lines_B321a, kernels_SPT, _PL);

      // B321b
      // propagators
      Propagator prop_B321b_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B321b_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_B321b_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B321b_11(Vertices::v1, Vertices::v1, prop_B321b_q);
      Line line_B321b_12(Vertices::v1, Vertices::v2, prop_B321b_k23);
      Line line_B321b_23(Vertices::v2, Vertices::v3, prop_B321b_k3);
      vector<Line> lines_B321b {line_B321b_11, line_B321b_12, line_B321b_23};
      // define the diagram
      Diagram* B321b = new Diagram(lines_B321b, kernels_SPT, _PL);

      // B222
      // propagators
      Propagator prop_B222_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B222_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
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
      Propagator prop_B321bx_k23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}, {Momenta::k3, 1}});
      Propagator prop_B321bx_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
      // lines
      Line line_B321bx_12(Vertices::v1, Vertices::v2, prop_B321bx_k23);
      Line line_B321bx_23(Vertices::v2, Vertices::v3, prop_B321bx_k3);
      vector<Line> lines_B321bx {line_B321bx_12, line_B321bx_23};
      // define the diagram
      Diagram* B321bx = new Diagram(lines_B321bx, kernels_EFTSPT, _PL);
      B321bx->set_perms(B321b->get_perms());
       
      // define the counterterm diagrams
      _cterms = {B411x,B321bx};
      _diagrams[Graphs::B411x] = B411x;
      _diagrams[Graphs::B321bx] = B321bx;
   }
}

//------------------------------------------------------------------------------
double Bispectrum::treeLevel_value(ThreeVector k2, ThreeVector k3)
{
   // set the external momenta
    _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::q, ThreeVector(1,0,0)}});

   double value = 0;
   // sum over diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(_momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Bispectrum::oneLoopSPT_value(ThreeVector k2, ThreeVector k3, ThreeVector q)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::q, q}});

   double value = 0;
   // sum over diagrams
   for(unsigned int i=0; i<_loop.size(); i++) {
      value += _loop[i]->value_IRreg(_momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double Bispectrum::oneLoopSPT_value(ThreeVector k2, ThreeVector k3)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::q, ThreeVector(1,0,0)}});
 
   double value = 0;
   // Integration

   return value;
}

//------------------------------------------------------------------------------
double Bispectrum::oneLoopCterms_value(ThreeVector k2, ThreeVector k3)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k2-k3}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::q, ThreeVector(1,0,0)}});

   double value = 0;
   // sum over diagrams
   for(size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(_momenta);
   }
   return value;
}


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

//------------------------------------------------------------------------------
PowerSpectrum::PowerSpectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients) : _order(order), _PL(PL), _SPTkernels(new SPTkernels), _EFTkernels(new EFTkernels), _eftcoefficients(eftcoefficients), _UVcutoff(10.), _kBin(0), _W(NULL)
{
   // Diagram momenta
    _labels = { Momenta::k1, Momenta::k2, Momenta::q };
   _momenta = DiagramMomenta(_labels);
    
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
      Line line_P22_11(Vertices::v1, Vertices::v2, prop_P22_q);
      Line line_P22_12(Vertices::v1, Vertices::v2, prop_P22_qk2);
      vector<Line> lines_P22 {line_P22_11, line_P22_12};
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
double PowerSpectrum::treeLevel_value(ThreeVector k)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k}, {Momenta::k2, k}, {Momenta::q, ThreeVector(1,0,0)}});

   double value = 0;
   // sum the tree level diagrams
   for (size_t i = 0; i < _tree.size(); i++) {
      value += _tree[i]->value_IRreg(_momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double PowerSpectrum::oneLoopSPT_value(ThreeVector k, ThreeVector q)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k}, {Momenta::k2, k}, {Momenta::q, q}});
    
   double value = 0;
   // sum the loop diagrams
   for (size_t i = 0; i < _loop.size(); i++) {
      value += _loop[i]->value_IRreg(_momenta);
   }
   return value;
}

//------------------------------------------------------------------------------
double PowerSpectrum::oneLoopSPT_value(ThreeVector k)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k}, {Momenta::k2, k}, {Momenta::q, ThreeVector(1,0,0)}});

   double value = 0;
   // Integration

   return value;
}

//------------------------------------------------------------------------------
double PowerSpectrum::oneLoopCterms_value(ThreeVector k)
{
   // set the external momenta
   _momenta.set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> {{Momenta::k1, -k}, {Momenta::k2, k}, {Momenta::q, ThreeVector(1,0,0)}});

   double value = 0;
   for (size_t i = 0; i < _cterms.size(); i++) {
      value += _cterms[i]->value_IRreg(_momenta);
   }
   return value;
}

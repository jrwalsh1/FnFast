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
Bispectrum::Bispectrum(Order order, LinearPowerSpectrumBase* PL) : _order(order), _PL(PL), _SPTkernels(new SPTkernels)
{
   // vertex kernels
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_SPT = {{Vertices::v1, _SPTkernels}, {Vertices::v2, _SPTkernels}, {Vertices::v3, _SPTkernels}};

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
   Diagram B211(lines_B211, kernels_SPT, _PL);

   // define the tree diagrams
   _tree = {B211};

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
      Diagram B411(lines_B411, kernels_SPT, _PL);

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
      Diagram B321a(lines_B321a, kernels_SPT, _PL);

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
      Diagram B321b(lines_B321b, kernels_SPT, _PL);

      // B222
      // propagators
      Propagator prop_B222_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
      Propagator prop_B222_qk2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}});
      Propagator prop_B222_qk23(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, -1}, {Momenta::k2, 1}, {Momenta::k3, 1}});
      // lines
      Line line_B222_12(Vertices::v1, Vertices::v1, prop_B222_q);
      Line line_B222_23(Vertices::v2, Vertices::v3, prop_B222_qk2);
      Line line_B222_13(Vertices::v1, Vertices::v3, prop_B222_qk23);
      vector<Line> lines_B222 {line_B222_12, line_B222_23, line_B222_13};
      // define the diagram
      Diagram B222(lines_B222, kernels_SPT, _PL);

      // define the loop diagrams
      _loop = {B411, B321a, B321b, B222};
   }
}

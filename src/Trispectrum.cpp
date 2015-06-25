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
Trispectrum::Trispectrum(Order order, LinearPowerSpectrumBase* PL) : _order(order), _PL(PL), _SPTkernels(new SPTkernels)
{
   // vertex kernels
   unordered_map<Vertices::VertexLabel, KernelBase*> kernels_SPT = {{Vertices::v1, _SPTkernels}, {Vertices::v2, _SPTkernels}, {Vertices::v3, _SPTkernels}, {Vertices::v4, _SPTkernels}};

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
   Diagram T3111(lines_T3111, kernels_SPT, _PL);

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
   Diagram T2211(lines_T2211, kernels_SPT, _PL);

   // define the tree diagrams
   _tree = {T3111, T2211};

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
      Diagram T5111(lines_T5111, kernels_SPT, _PL);

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
      Diagram T4211a(lines_T4211a, kernels_SPT, _PL);

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
      Diagram T4211b(lines_T4211b, kernels_SPT, _PL);

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
      Diagram T3311a(lines_T3311a, kernels_SPT, _PL);

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
      Diagram T3311b(lines_T3311b, kernels_SPT, _PL);

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
      Diagram T3221a(lines_T3221a, kernels_SPT, _PL);

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
      Diagram T3221b(lines_T3221b, kernels_SPT, _PL);

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
      Diagram T3221c(lines_T3221c, kernels_SPT, _PL);

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
      Diagram T2222(lines_T2222, kernels_SPT, _PL);

      // define the loop diagrams
      _loop = {T5111, T4211a, T4211b, T3311a, T3311b, T3221a, T3221b, T3221c, T2222};
   }
}

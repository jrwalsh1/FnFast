//------------------------------------------------------------------------------
/// \file DiagramSet2point.cpp
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
//    Implementation of class DiagramSet2point
//------------------------------------------------------------------------------

#include <iostream>

#include "DiagramSet2point.hpp"

//------------------------------------------------------------------------------
DiagramSet2point::DiagramSet2point(Order order) : DiagramSetBase(order)
{
   // External momentum labels of the diagrams
   _extmomlabels = { MomentumLabel::k1, MomentumLabel::k2};

   // set up the diagrams, starting with the tree level
   // P11
   // propagators
   Propagator prop_P11_k2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_P11_12(VertexLabel::v1, VertexLabel::v2, prop_P11_k2);
   vector<Line> lines_P11 {line_P11_12};
   // define the diagram
   DiagramTree* P11 = new DiagramTree(lines_P11);

   // define the tree diagrams
   _tree = {P11};
   _diagrams[Graphs_2point::P11] = P11;

   // add the one loop graphs (if requested)
   if ((_order == Order::kOneLoop) || (_order == Order::kTwoLoop)) {
      // P31
      // propagators
      Propagator prop_P31_q(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_P31_k2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_P31_11(VertexLabel::v1, VertexLabel::v1, prop_P31_q);
      Line line_P31_12(VertexLabel::v1, VertexLabel::v2, prop_P31_k2);
      vector<Line> lines_P31 {line_P31_11, line_P31_12};
      // define the diagram
      DiagramOneLoop* P31 = new DiagramOneLoop(lines_P31);

      // P22
      // propagators
      Propagator prop_P22_q(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_P22_qk2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kMinus}, {MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_P22_12a(VertexLabel::v1, VertexLabel::v2, prop_P22_q);
      Line line_P22_12b(VertexLabel::v1, VertexLabel::v2, prop_P22_qk2);
      vector<Line> lines_P22 {line_P22_12a, line_P22_12b};
      // define the diagram
      DiagramOneLoop* P22 = new DiagramOneLoop(lines_P22);

      // define the loop diagrams
      _oneLoop = {P31, P22};
      _diagrams[Graphs_2point::P31] = P31;
      _diagrams[Graphs_2point::P22] = P22;
   }

   // add the two loop graphs (if requested)
   if (_order == Order::kTwoLoop) {
      // P51
      // propagators
      Propagator prop_P51_q(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_P51_q2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q2, Propagator::LabelFlow::kPlus}});
      Propagator prop_P51_k2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_P51_11a(VertexLabel::v1, VertexLabel::v1, prop_P51_q);
      Line line_P51_11b(VertexLabel::v1, VertexLabel::v1, prop_P51_q2);
      Line line_P51_12(VertexLabel::v1, VertexLabel::v2, prop_P51_k2);
      vector<Line> lines_P51 {line_P51_11a, line_P51_11b, line_P51_12};
      // define the diagram
      DiagramTwoLoop* P51 = new DiagramTwoLoop(lines_P51);

      // P42
      // propagators
      Propagator prop_P42_q(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_P42_q2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q2, Propagator::LabelFlow::kPlus}});
      Propagator prop_P42_qk2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kMinus}, {MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_P42_11(VertexLabel::v1, VertexLabel::v1, prop_P42_q2);
      Line line_P42_12a(VertexLabel::v1, VertexLabel::v2, prop_P42_q);
      Line line_P42_12b(VertexLabel::v1, VertexLabel::v2, prop_P42_qk2);
      vector<Line> lines_P42 {line_P42_11, line_P42_12a, line_P42_12b};
      // define the diagram
      DiagramTwoLoop* P42 = new DiagramTwoLoop(lines_P42);

      // P33a
      // propagators
      Propagator prop_P33a_q(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_P33a_q2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q2, Propagator::LabelFlow::kPlus}});
      Propagator prop_P33a_k2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_P33a_11(VertexLabel::v1, VertexLabel::v1, prop_P33a_q);
      Line line_P33a_22(VertexLabel::v2, VertexLabel::v2, prop_P33a_q2);
      Line line_P33a_12(VertexLabel::v1, VertexLabel::v2, prop_P33a_k2);
      vector<Line> lines_P33a {line_P33a_11, line_P33a_22, line_P33a_12};
      // define the diagram
      DiagramTwoLoop* P33a = new DiagramTwoLoop(lines_P33a);

      // P33b
      // propagators
      Propagator prop_P33b_q(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_P33b_q2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q2, Propagator::LabelFlow::kPlus}});
      Propagator prop_P33b_qq2k2(MomentumMap<Propagator::LabelFlow> {{MomentumLabel::q, Propagator::LabelFlow::kMinus}, {MomentumLabel::q2, Propagator::LabelFlow::kMinus}, {MomentumLabel::k2, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_P33b_12a(VertexLabel::v1, VertexLabel::v2, prop_P33b_q);
      Line line_P33b_12b(VertexLabel::v1, VertexLabel::v2, prop_P33b_q2);
      Line line_P33b_12c(VertexLabel::v1, VertexLabel::v2, prop_P33b_qq2k2);
      vector<Line> lines_P33b {line_P33b_12a, line_P33b_12b, line_P33b_12c};
      // define the diagram
      DiagramTwoLoop* P33b = new DiagramTwoLoop(lines_P33b);

      // define the two loop diagrams
      _twoLoop = {P51, P42, P33a, P33b};
      _diagrams[Graphs_2point::P51] = P51;
      _diagrams[Graphs_2point::P42] = P42;
      _diagrams[Graphs_2point::P33a] = P33a;
      _diagrams[Graphs_2point::P33b] = P33b;
   }
}

//------------------------------------------------------------------------------
/// \file DiagramSet3point.cpp
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
//    Implementation of class DiagramSet3point
//------------------------------------------------------------------------------

#include <iostream>

#include "DiagramSet3point.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramSet3point::DiagramSet3point(Order order) : DiagramSetBase(order)
{
   // External momentum labels of the diagrams
   _extmomlabels = { Momentum::k1, Momentum::k2, Momentum::k3};

   // set up the diagrams, starting with the tree level
   // B211
   // propagators
   Propagator prop_B211_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}});
   Propagator prop_B211_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_B211_12(Vertex::v1, Vertex::v2, prop_B211_k2);
   Line line_B211_13(Vertex::v1, Vertex::v3, prop_B211_k3);
   std::vector<Line> lines_B211 {line_B211_12, line_B211_13};
   // define the diagram
   DiagramTree* B211 = new DiagramTree(lines_B211);

   // define the tree diagrams
   _tree = {B211};
   _diagrams[Graphs_3point::B211] = B211;

   if (_order == Order::kOneLoop) {
      // B411
      // propagators
      Propagator prop_B411_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_B411_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}});
      Propagator prop_B411_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_B411_11(Vertex::v1, Vertex::v1, prop_B411_q);
      Line line_B411_12(Vertex::v1, Vertex::v2, prop_B411_k2);
      Line line_B411_13(Vertex::v1, Vertex::v3, prop_B411_k3);
      std::vector<Line> lines_B411 {line_B411_11, line_B411_12, line_B411_13};
      // define the diagram
      DiagramOneLoop* B411 = new DiagramOneLoop(lines_B411);

      // B321a
      // propagators
      Propagator prop_B321a_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_B321a_k23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_B321a_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_B321a_11(Vertex::v1, Vertex::v1, prop_B321a_q);
      Line line_B321a_12(Vertex::v1, Vertex::v2, prop_B321a_k23);
      Line line_B321a_23(Vertex::v2, Vertex::v3, prop_B321a_k3);
      std::vector<Line> lines_B321a {line_B321a_11, line_B321a_12, line_B321a_23};
      // define the diagram
      DiagramOneLoop* B321a = new DiagramOneLoop(lines_B321a);

      // B321b
      // propagators
      Propagator prop_B321b_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_B321b_qk2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}});
      Propagator prop_B321b_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_B321b_12a(Vertex::v1, Vertex::v2, prop_B321b_q);
      Line line_B321b_12b(Vertex::v1, Vertex::v2, prop_B321b_qk2);
      Line line_B321b_13(Vertex::v1, Vertex::v3, prop_B321b_k3);
      std::vector<Line> lines_B321b {line_B321b_12a, line_B321b_12b, line_B321b_13};
      // define the diagram
      DiagramOneLoop* B321b = new DiagramOneLoop(lines_B321b);

      // B222
      // propagators
      Propagator prop_B222_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_B222_qk2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}, {Momentum::k2, Propagator::LabelFlow::kMinus}});
      Propagator prop_B222_qk23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_B222_12(Vertex::v1, Vertex::v2, prop_B222_q);
      Line line_B222_23(Vertex::v2, Vertex::v3, prop_B222_qk2);
      Line line_B222_13(Vertex::v1, Vertex::v3, prop_B222_qk23);
      std::vector<Line> lines_B222 {line_B222_12, line_B222_23, line_B222_13};
      // define the diagram
      DiagramOneLoop* B222 = new DiagramOneLoop(lines_B222);

      // define the loop diagrams
      _oneLoop = {B411, B321a, B321b, B222};
      _diagrams[Graphs_3point::B411] = B411;
      _diagrams[Graphs_3point::B321a] = B321a;
      _diagrams[Graphs_3point::B321b] = B321b;
      _diagrams[Graphs_3point::B222] = B222;
   }
}

} // namespace fnfast

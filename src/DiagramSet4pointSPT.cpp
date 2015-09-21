//------------------------------------------------------------------------------
/// \file DiagramSet4pointSPT.cpp
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
//    Implementation of class DiagramSet4pointSPT
//------------------------------------------------------------------------------

#include <iostream>

#include "DiagramSet4pointSPT.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramSet4pointSPT::DiagramSet4pointSPT(Order order, LabelMap<Vertex, KernelType> kerneltypes)
: DiagramSetBase(order), _kerneltypes(kerneltypes)
{
   // External momentum labels of the diagrams
   _extmomlabels = { Momentum::k1, Momentum::k2, Momentum::k3, Momentum::k4};

   // set the vertex and kernel types
   _vertextypes = LabelMap<Vertex, VertexType>({{Vertex::v1, VertexType::type1}, {Vertex::v2, VertexType::type1}, {Vertex::v3, VertexType::type1}, {Vertex::v4, VertexType::type1}});

   // set up the diagrams, starting with the tree level
   // T3111
   // propagators
   Propagator prop_T3111_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}});
   Propagator prop_T3111_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T3111_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_T3111_12(Vertex::v1, Vertex::v2, prop_T3111_k2);
   Line line_T3111_13(Vertex::v1, Vertex::v3, prop_T3111_k3);
   Line line_T3111_14(Vertex::v1, Vertex::v4, prop_T3111_k4);
   std::vector<Line> lines_T3111 {line_T3111_12, line_T3111_13, line_T3111_14};
   // define the diagram
   DiagramTree* T3111 = new DiagramTree(lines_T3111, _vertextypes, _kerneltypes);

   // T2211
   // propagators
   Propagator prop_T2211_k23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T2211_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T2211_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_T2211_12(Vertex::v1, Vertex::v2, prop_T2211_k23);
   Line line_T2211_23(Vertex::v2, Vertex::v3, prop_T2211_k3);
   Line line_T2211_14(Vertex::v1, Vertex::v4, prop_T2211_k4);
   std::vector<Line> lines_T2211 {line_T2211_12, line_T2211_23, line_T2211_14};
   // define the diagram
   DiagramTree* T2211 = new DiagramTree(lines_T2211, _vertextypes, _kerneltypes);

   // define the tree diagrams
   _tree = {T3111, T2211};
   _diagrams = LabelMap<Graphs_4point, DiagramBase*> {{Graphs_4point::T3111, T3111}, {Graphs_4point::T2211, T2211}};

   if (_order == Order::kOneLoop) {
      // T5111
      // propagators
      Propagator prop_T5111_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T5111_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}});
      Propagator prop_T5111_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T5111_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T5111_11(Vertex::v1, Vertex::v1, prop_T5111_q);
      Line line_T5111_12(Vertex::v1, Vertex::v2, prop_T5111_k2);
      Line line_T5111_13(Vertex::v1, Vertex::v3, prop_T5111_k3);
      Line line_T5111_14(Vertex::v1, Vertex::v4, prop_T5111_k4);
      std::vector<Line> lines_T5111 {line_T5111_11, line_T5111_12, line_T5111_13, line_T5111_14};
      // define the diagram
      DiagramOneLoop* T5111 = new DiagramOneLoop(lines_T5111, _vertextypes, _kerneltypes);

      // T4211a
      // propagators
      Propagator prop_T4211a_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T4211a_k23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T4211a_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T4211a_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T4211a_11(Vertex::v1, Vertex::v1, prop_T4211a_q);
      Line line_T4211a_12(Vertex::v1, Vertex::v2, prop_T4211a_k23);
      Line line_T4211a_23(Vertex::v2, Vertex::v3, prop_T4211a_k3);
      Line line_T4211a_14(Vertex::v1, Vertex::v4, prop_T4211a_k4);
      std::vector<Line> lines_T4211a {line_T4211a_11, line_T4211a_12, line_T4211a_23, line_T4211a_14};
      // define the diagram
      DiagramOneLoop* T4211a = new DiagramOneLoop(lines_T4211a, _vertextypes, _kerneltypes);

      // T4211b
      // propagators
      Propagator prop_T4211b_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T4211b_qk2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}});
      Propagator prop_T4211b_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T4211b_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T4211b_12a(Vertex::v1, Vertex::v2, prop_T4211b_q);
      Line line_T4211b_12b(Vertex::v1, Vertex::v2, prop_T4211b_qk2);
      Line line_T4211b_13(Vertex::v1, Vertex::v3, prop_T4211b_k3);
      Line line_T4211b_14(Vertex::v1, Vertex::v4, prop_T4211b_k4);
      std::vector<Line> lines_T4211b {line_T4211b_12a, line_T4211b_12b, line_T4211b_13, line_T4211b_14};
      // define the diagram
      DiagramOneLoop* T4211b = new DiagramOneLoop(lines_T4211b, _vertextypes, _kerneltypes);

      // T3311a
      // propagators
      Propagator prop_T3311a_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3311a_k234(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3311a_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3311a_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T3311a_11(Vertex::v1, Vertex::v1, prop_T3311a_q);
      Line line_T3311a_12(Vertex::v1, Vertex::v2, prop_T3311a_k234);
      Line line_T3311a_23(Vertex::v2, Vertex::v3, prop_T3311a_k3);
      Line line_T3311a_24(Vertex::v2, Vertex::v4, prop_T3311a_k4);
      std::vector<Line> lines_T3311a {line_T3311a_11, line_T3311a_12, line_T3311a_23, line_T3311a_24};
      // define the diagram
      DiagramOneLoop* T3311a = new DiagramOneLoop(lines_T3311a, _vertextypes, _kerneltypes);

      // T3311b
      // propagators
      Propagator prop_T3311b_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3311b_qk23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3311b_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3311b_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T3311b_12a(Vertex::v1, Vertex::v2, prop_T3311b_q);
      Line line_T3311b_12b(Vertex::v1, Vertex::v2, prop_T3311b_qk23);
      Line line_T3311b_23(Vertex::v2, Vertex::v3, prop_T3311b_k3);
      Line line_T3311b_14(Vertex::v1, Vertex::v4, prop_T3311b_k4);
      std::vector<Line> lines_T3311b {line_T3311b_12a, line_T3311b_12b, line_T3311b_23, line_T3311b_14};
      // define the diagram
      DiagramOneLoop* T3311b = new DiagramOneLoop(lines_T3311b, _vertextypes, _kerneltypes);

      // T3221a
      // propagators
      Propagator prop_T3221a_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221a_k234(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221a_k34(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221a_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T3221a_11(Vertex::v1, Vertex::v1, prop_T3221a_q);
      Line line_T3221a_12(Vertex::v1, Vertex::v2, prop_T3221a_k234);
      Line line_T3221a_23(Vertex::v2, Vertex::v3, prop_T3221a_k34);
      Line line_T3221a_34(Vertex::v3, Vertex::v4, prop_T3221a_k4);
      std::vector<Line> lines_T3221a {line_T3221a_11, line_T3221a_12, line_T3221a_23, line_T3221a_34};
      // define the diagram
      DiagramOneLoop* T3221a = new DiagramOneLoop(lines_T3221a, _vertextypes, _kerneltypes);

      // T3221b
      // propagators
      Propagator prop_T3221b_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221b_qk2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221b_k34(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221b_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T3221b_12a(Vertex::v1, Vertex::v2, prop_T3221b_q);
      Line line_T3221b_12b(Vertex::v1, Vertex::v2, prop_T3221b_qk2);
      Line line_T3221b_13(Vertex::v1, Vertex::v3, prop_T3221b_k34);
      Line line_T3221b_34(Vertex::v3, Vertex::v4, prop_T3221b_k4);
      std::vector<Line> lines_T3221b {line_T3221b_12a, line_T3221b_12b, line_T3221b_13, line_T3221b_34};
      // define the diagram
      DiagramOneLoop* T3221b = new DiagramOneLoop(lines_T3221b, _vertextypes, _kerneltypes);

      // T3221c
      // propagators
      Propagator prop_T3221c_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221c_qk2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}, {Momentum::k2, Propagator::LabelFlow::kMinus}});
      Propagator prop_T3221c_qk23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
      Propagator prop_T3221c_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T3221c_12(Vertex::v1, Vertex::v2, prop_T3221c_q);
      Line line_T3221c_23(Vertex::v2, Vertex::v3, prop_T3221c_qk2);
      Line line_T3221c_13(Vertex::v1, Vertex::v3, prop_T3221c_qk23);
      Line line_T3221c_14(Vertex::v1, Vertex::v4, prop_T3221c_k4);
      std::vector<Line> lines_T3221c {line_T3221c_12, line_T3221c_23, line_T3221c_13, line_T3221c_14};
      // define the diagram
      DiagramOneLoop* T3221c = new DiagramOneLoop(lines_T3221c, _vertextypes, _kerneltypes);

      // T2222
      // propagators
      Propagator prop_T2222_q(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}});
      Propagator prop_T2222_qk2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}, {Momentum::k2, Propagator::LabelFlow::kMinus}});
      Propagator prop_T2222_qk23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kPlus}, {Momentum::k2, Propagator::LabelFlow::kMinus}, {Momentum::k3, Propagator::LabelFlow::kMinus}});
      Propagator prop_T2222_qk234(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::q, Propagator::LabelFlow::kMinus}, {Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
      // lines
      Line line_T2222_12(Vertex::v1, Vertex::v2, prop_T2222_q);
      Line line_T2222_23(Vertex::v2, Vertex::v3, prop_T2222_qk2);
      Line line_T2222_34(Vertex::v3, Vertex::v4, prop_T2222_qk23);
      Line line_T2222_14(Vertex::v1, Vertex::v4, prop_T2222_qk234);
      std::vector<Line> lines_T2222 {line_T2222_12, line_T2222_23, line_T2222_34, line_T2222_14};
      // define the diagram
      DiagramOneLoop* T2222 = new DiagramOneLoop(lines_T2222, _vertextypes, _kerneltypes);

      // define the loop diagrams
      _oneLoop = {T5111, T4211a, T4211b, T3311a, T3311b, T3221a, T3221b, T3221c, T2222};
      _diagrams = LabelMap<Graphs_4point, DiagramBase*> {{Graphs_4point::T3111, T3111}, {Graphs_4point::T2211, T2211}, {Graphs_4point::T5111, T5111}, {Graphs_4point::T4211a, T4211a}, {Graphs_4point::T4211b, T4211b},
            {Graphs_4point::T3311a, T3311a}, {Graphs_4point::T3311b, T3311b}, {Graphs_4point::T3221a, T3221a}, {Graphs_4point::T3221b, T3221b}, {Graphs_4point::T3221c, T3221c}, {Graphs_4point::T2222, T2222}};
   }
}

} // namespace fnfast

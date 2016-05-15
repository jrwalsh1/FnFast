//------------------------------------------------------------------------------
/// \file DiagramSet4pointEFT.cpp
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
//    Implementation of class DiagramSet4pointEFT
//------------------------------------------------------------------------------

#include <iostream>

#include "DiagramSet4pointEFT.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramSet4pointEFT::DiagramSet4pointEFT(Order order, LabelMap<Vertex, KernelType> kerneltypes)
: DiagramSetBase(order), _kerneltypes(kerneltypes)
{
   // External momentum labels of the diagrams
   _extmomlabels = { Momentum::k1, Momentum::k2, Momentum::k3, Momentum::k4};

   // set the vertex and kernel types
   _vertextypes = LabelMap<Vertex, VertexType>({{Vertex::v1, VertexType::type1}, {Vertex::v2, VertexType::type2}, {Vertex::v3, VertexType::type2}, {Vertex::v4, VertexType::type2}});

   // get perms from SPT diagram
   /*DAN*/
   DiagramSet4pointSPT SPTdiagrams(Order::kOneLoop);
   std::vector<LabelMap<Momentum, Momentum> > perms;
   
   // set up the diagrams, starting with the tree level
   // T5111x
   // propagators
   Propagator prop_T5111x_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}});
   Propagator prop_T5111x_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T5111x_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_T5111x_12(Vertex::v1, Vertex::v2, prop_T5111x_k2);
   Line line_T5111x_13(Vertex::v1, Vertex::v3, prop_T5111x_k3);
   Line line_T5111x_14(Vertex::v1, Vertex::v4, prop_T5111x_k4);
   std::vector<Line> lines_T5111x {line_T5111x_12, line_T5111x_13, line_T5111x_14};
   // define the diagram
   DiagramTree* T5111x = new DiagramTree(lines_T5111x, _vertextypes, _kerneltypes);
   /*DAN*/
   perms = SPTdiagrams[Graphs_4point::T5111]->get_perms();
   T5111x->set_perms(perms);
   

   // T4211ax
   // propagators
   Propagator prop_T4211ax_k23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T4211ax_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T4211ax_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_T4211ax_12(Vertex::v1, Vertex::v2, prop_T4211ax_k23);
   Line line_T4211ax_23(Vertex::v2, Vertex::v3, prop_T4211ax_k3);
   Line line_T4211ax_14(Vertex::v1, Vertex::v4, prop_T4211ax_k4);
   std::vector<Line> lines_T4211ax {line_T4211ax_12, line_T4211ax_23, line_T4211ax_14};
   // define the diagram
   DiagramTree* T4211ax = new DiagramTree(lines_T4211ax, _vertextypes, _kerneltypes);
   /*DAN*/
   perms = SPTdiagrams[Graphs_4point::T4211a]->get_perms();
   T4211ax->set_perms(perms);

   // T3311ax
   // propagators
   Propagator prop_T3311ax_k234(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
   Propagator prop_T3311ax_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}});
   Propagator prop_T3311ax_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_T3311ax_12(Vertex::v1, Vertex::v2, prop_T3311ax_k234);
   Line line_T3311ax_23(Vertex::v2, Vertex::v3, prop_T3311ax_k3);
   Line line_T3311ax_24(Vertex::v2, Vertex::v4, prop_T3311ax_k4);
   std::vector<Line> lines_T3311ax {line_T3311ax_12, line_T3311ax_23, line_T3311ax_24};
   // define the diagram
   DiagramTree* T3311ax = new DiagramTree(lines_T3311ax, _vertextypes, _kerneltypes);
   /*DAN*/
   perms = SPTdiagrams[Graphs_4point::T3311a]->get_perms();
   T3311ax->set_perms(perms);

   // T3221ax
   // propagators
   Propagator prop_T3221ax_k234(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}, {Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
   Propagator prop_T3221ax_k34(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::kPlus}, {Momentum::k4, Propagator::LabelFlow::kPlus}});
   Propagator prop_T3221ax_k4(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k4, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_T3221ax_12(Vertex::v1, Vertex::v2, prop_T3221ax_k234);
   Line line_T3221ax_23(Vertex::v2, Vertex::v3, prop_T3221ax_k34);
   Line line_T3221ax_34(Vertex::v3, Vertex::v4, prop_T3221ax_k4);
   std::vector<Line> lines_T3221ax {line_T3221ax_12, line_T3221ax_23, line_T3221ax_34};
   // define the diagram
   DiagramTree* T3221ax = new DiagramTree(lines_T3221ax, _vertextypes, _kerneltypes);
   /*DAN*/
   perms = SPTdiagrams[Graphs_4point::T3221a]->get_perms();
   T3221ax->set_perms(perms);

   // define the tree diagrams
   _tree = {T5111x, T4211ax, T3311ax, T3221ax};
   _diagrams = LabelMap<Graphs_4point, DiagramBase*> {{Graphs_4point::T5111x, T5111x}, {Graphs_4point::T4211ax, T4211ax}, {Graphs_4point::T3311ax, T3311ax}, {Graphs_4point::T3221ax, T3221ax}};
}

} // namespace fnfast

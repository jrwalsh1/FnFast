//------------------------------------------------------------------------------
/// \file DiagramSet3pointEFT.cpp
//
// Author(s):
//    Jon Walsh
//
// Copyright:
//    Copyright (C) 2015  LBL
//
//    This file is part of the FnFast library. FnFast is distributed under the
//    terms of the GNU General Public License version 3 (GPLv3), see the COPYING
//    file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Implementation of class DiagramSet3pointEFT
//------------------------------------------------------------------------------

#include <iostream>

#include "DiagramSet3pointEFT.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramSet3pointEFT::DiagramSet3pointEFT(Order order, LabelMap<Vertex, KernelType> kerneltypes)
: DiagramSetBase(order), _kerneltypes(kerneltypes)
{
   // External momentum labels of the diagrams
   _extmomlabels = { Momentum::k1, Momentum::k2, Momentum::k3};

   // set the vertex and kernel types
   _vertextypes = LabelMap<Vertex, VertexType>({{Vertex::v1, VertexType::type1}, {Vertex::v2, VertexType::type2}, {Vertex::v3, VertexType::type2}});

   // set up the diagrams, starting with the tree level
   // B411x
   // propagators
   Propagator prop_B411x_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::Plus}});
   Propagator prop_B411x_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::Plus}});
   // lines
   Line line_B411x_12(Vertex::v1, Vertex::v2, prop_B411x_k2);
   Line line_B411x_13(Vertex::v1, Vertex::v3, prop_B411x_k3);
   std::vector<Line> lines_B411x {line_B411x_12, line_B411x_13};
   // define the diagram
   DiagramTree* B411x = new DiagramTree(lines_B411x, _vertextypes, _kerneltypes);

   // B321ax
   // propagators
   Propagator prop_B321ax_k23(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::Plus}, {Momentum::k3, Propagator::LabelFlow::Plus}});
   Propagator prop_B321ax_k3(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k3, Propagator::LabelFlow::Plus}});
   // lines
   Line line_B321ax_12(Vertex::v1, Vertex::v2, prop_B321ax_k23);
   Line line_B321ax_23(Vertex::v2, Vertex::v3, prop_B321ax_k3);
   std::vector<Line> lines_B321ax {line_B321ax_12, line_B321ax_23};
   // define the diagram
   DiagramTree* B321ax = new DiagramTree(lines_B321ax, _vertextypes, _kerneltypes);

   // define the tree diagrams
   _tree = {B411x, B321ax};
   _diagrams = LabelMap<Graphs_3point, DiagramBase*> {{Graphs_3point::B411x, B411x}, {Graphs_3point::B321ax, B321ax}};
}

} // namespace fnfast

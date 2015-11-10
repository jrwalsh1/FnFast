//------------------------------------------------------------------------------
/// \file DiagramSet2pointEFT.cpp
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
//    Implementation of class DiagramSet2pointEFT
//------------------------------------------------------------------------------

#include <iostream>

#include "DiagramSet2pointEFT.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramSet2pointEFT::DiagramSet2pointEFT(Order order, LabelMap<Vertex, KernelType> kerneltypes)
: DiagramSetBase(order), _kerneltypes(kerneltypes)
{
   // External momentum labels of the diagrams
   _extmomlabels = { Momentum::k1, Momentum::k2};

   // set the vertex and kernel types
   _vertextypes = LabelMap<Vertex, VertexType>({{Vertex::v1, VertexType::type1}, {Vertex::v2, VertexType::type2}});

   // set up the diagrams, starting with the tree level
   // P31x
   // propagators
   Propagator prop_P31x_k2(LabelMap<Momentum, Propagator::LabelFlow> {{Momentum::k2, Propagator::LabelFlow::kPlus}});
   // lines
   Line line_P31x_12(Vertex::v1, Vertex::v2, prop_P31x_k2);
   std::vector<Line> lines_P31x {line_P31x_12};
   // define the diagram
   DiagramTree* P31x = new DiagramTree(lines_P31x, _vertextypes, _kerneltypes);

   // define the tree diagrams
   _tree = {P31x};
   _diagrams = LabelMap<Graphs_2point, DiagramBase*> {{Graphs_2point::P31x, P31x}};
}

} // namespace fnfast

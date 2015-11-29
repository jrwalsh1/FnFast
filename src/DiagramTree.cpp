//------------------------------------------------------------------------------
/// \file DiagramTree.cpp
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
//    Implementation of class DiagramTree
//------------------------------------------------------------------------------

#include <cassert>
#include <iostream>

#include "DiagramTree.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramTree::DiagramTree(std::vector<Line> lines) : DiagramBase(lines)
{
   _order = Order::Tree;
   // check to ensure that the diagram is really tree level (no loop momentum)
   bool isLoop = false;
   for (auto line : _lines) {
      if (line.propagator.hasLabel(Momentum::q) || line.propagator.hasLabel(Momentum::q)) {
         isLoop = true;
      }
   }
   assert(!isLoop);
}

//------------------------------------------------------------------------------
DiagramTree::DiagramTree(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes) : DiagramBase(lines, vertextypes)
{
   _order = Order::Tree;
   // check to ensure that the diagram is really tree level (no loop momentum)
   bool isLoop = false;
   for (auto line : _lines) {
      if (line.propagator.hasLabel(Momentum::q) || line.propagator.hasLabel(Momentum::q)) {
         isLoop = true;
      }
   }
   assert(!isLoop);
}

//------------------------------------------------------------------------------
DiagramTree::DiagramTree(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes, LabelMap<Vertex, KernelType> kerneltypes) : DiagramBase(lines, vertextypes, kerneltypes)
{
   _order = Order::Tree;
   // check to ensure that the diagram is really tree level (no loop momentum)
   bool isLoop = false;
   for (auto line : _lines) {
      if (line.propagator.hasLabel(Momentum::q) || line.propagator.hasLabel(Momentum::q)) {
         isLoop = true;
      }
   }
   assert(!isLoop);
}

//------------------------------------------------------------------------------
double DiagramTree::value(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // the diagram value is:
   // symmetry factor * propagators * vertices
   // summed over external momentum permutations
   double value = 0;
   for (auto perm : _perms) {
      // get the momentum permutation
      LabelMap<Momentum, ThreeVector> mom_perm = mom;
      mom_perm.permute(perm);
      double diagvalue = _symfac;
      // iterate over lines
      for (auto line : _lines) {
         diagvalue *= (*PL)(line.propagator.p(mom_perm).magnitude());
      }
      // now do vertex factors
      for (auto vertex : _vertices) {
         std::vector<ThreeVector> p;
         p.reserve(_vertexmomenta[vertex].size());
         // loop over propagators attached to the vertex
         for (auto vx_prop : _vertexmomenta[vertex]) {
            p.push_back(vx_prop.p(mom_perm));
         }
         if (_kerneltypes[vertex] == KernelType::delta) {
            diagvalue *= kernels[vertex]->Fn_sym(p);
         } else {
            diagvalue *= kernels[vertex]->Gn_sym(p);
         }
      }
      value += diagvalue;
   }
   return value;
}

} // namespace fnfast

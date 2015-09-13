//------------------------------------------------------------------------------
/// \file DiagramTree.cpp
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
//    Implementation of class DiagramTree
//------------------------------------------------------------------------------

#include <cassert>

#include "DiagramTree.hpp"

//------------------------------------------------------------------------------
DiagramTree::DiagramTree(vector<Line> lines) : DiagramBase(lines)
{
   _order = Order::kTree;
   // check to ensure that the diagram is really tree level (no loop momentum)
   bool isLoop = false;
   for (auto line : _lines) {
      if (line.propagator.hasLabel(MomentumLabel::q) || line.propagator.hasLabel(MomentumLabel::q)) {
         isLoop = true;
      }
   }
   assert(!isLoop);
}

//------------------------------------------------------------------------------
double DiagramTree::value(const MomentumMap<ThreeVector>& mom, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // the diagram value is:
   // symmetry factor * propagators * vertices
   double value = _symfac;
   // iterate over lines
   for (auto line : _lines) {
      value *= (*PL)(line.propagator.p(mom).magnitude());
   }
   // now do vertex factors
   for (auto vertex : _vertices) {
      vector<ThreeVector> p;
      p.reserve(_vertexmomenta[vertex].size());
      // loop over propagators attached to the vertex
      for (auto vx_prop : _vertexmomenta[vertex]) {
         p.push_back(vx_prop.p(mom));
      }
      value *= kernels[vertex]->Fn_sym(p);
   }
   return value;
}

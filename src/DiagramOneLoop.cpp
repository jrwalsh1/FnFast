//------------------------------------------------------------------------------
/// \file DiagramOneLoop.cpp
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
//    Implementation of class DiagramOneLoop
//------------------------------------------------------------------------------

#include <limits>
#include <cassert>

#include "DiagramOneLoop.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramOneLoop::DiagramOneLoop(std::vector<Line> lines) : DiagramBase(lines), _qmax(std::numeric_limits<double>::infinity())
{
   _order = Order::kOneLoop;
   // check to ensure that the diagram is really one loop
   bool isLoop = false;
   bool is2Loop = false;
   // find the nontrivial poles
   for (auto line : _lines) {
      // check if the line has the loop momentum in it
      // if so, it has an IR pole that must be regulated if it is away from 0
      if (line.propagator.hasLabel(Momentum::q)) {
         isLoop = true;
         _order = Order::kOneLoop;
         Propagator pole = line.propagator.IRpole(Momentum::q);
         if (!pole.isNull()) {
            _IRpoles.push_back(pole);
         }
      }
      if (line.propagator.hasLabel(Momentum::q2)) {
         is2Loop = true;
      }
   }
   assert(isLoop && !is2Loop);
}

//------------------------------------------------------------------------------
DiagramOneLoop::DiagramOneLoop(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes) : DiagramBase(lines, vertextypes), _qmax(std::numeric_limits<double>::infinity())
{
   _order = Order::kOneLoop;
   // check to ensure that the diagram is really one loop
   bool isLoop = false;
   bool is2Loop = false;
   // find the nontrivial poles
   for (auto line : _lines) {
      // check if the line has the loop momentum in it
      // if so, it has an IR pole that must be regulated if it is away from 0
      if (line.propagator.hasLabel(Momentum::q)) {
         isLoop = true;
         _order = Order::kOneLoop;
         Propagator pole = line.propagator.IRpole(Momentum::q);
         if (!pole.isNull()) {
            _IRpoles.push_back(pole);
         }
      }
      if (line.propagator.hasLabel(Momentum::q2)) {
         is2Loop = true;
      }
   }
   assert(isLoop && !is2Loop);
}

//------------------------------------------------------------------------------
DiagramOneLoop::DiagramOneLoop(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes, LabelMap<Vertex, KernelType> kerneltypes) : DiagramBase(lines, vertextypes, kerneltypes), _qmax(std::numeric_limits<double>::infinity())
{
   _order = Order::kOneLoop;
   // check to ensure that the diagram is really one loop
   bool isLoop = false;
   bool is2Loop = false;
   // find the nontrivial poles
   for (auto line : _lines) {
      // check if the line has the loop momentum in it
      // if so, it has an IR pole that must be regulated if it is away from 0
      if (line.propagator.hasLabel(Momentum::q)) {
         isLoop = true;
         _order = Order::kOneLoop;
         Propagator pole = line.propagator.IRpole(Momentum::q);
         if (!pole.isNull()) {
            _IRpoles.push_back(pole);
         }
      }
      if (line.propagator.hasLabel(Momentum::q2)) {
         is2Loop = true;
      }
   }
   assert(isLoop && !is2Loop);
}

//------------------------------------------------------------------------------
double DiagramOneLoop::value_base(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // check to see if the loop momentum is above the cutoff, if so return 0
   if (mom[Momentum::q].magnitude() > _qmax) { return 0; }

   // the diagram value is:
   // symmetry factor * propagators * vertices
   double value = _symfac;
   // iterate over lines
   for (auto line : _lines) {
      value *= (*PL)(line.propagator.p(mom).magnitude());
   }
   // now do vertex factors
   for (auto vertex : _vertices) {
      std::vector<ThreeVector> p;
      p.reserve(_vertexmomenta[vertex].size());
      // loop over propagators attached to the vertex
      for (auto vx_prop : _vertexmomenta[vertex]) {
         p.push_back(vx_prop.p(mom));
      }
      if (_kerneltypes[vertex] == KernelType::delta) {
         value *= kernels[vertex]->Fn_sym(p);
      } else {
         value *= kernels[vertex]->Gn_sym(p);
      }
   }
   return value;
}

//------------------------------------------------------------------------------
double DiagramOneLoop::value_base_IRreg(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   // no IR regulation necessary if there are no poles away from q = 0
   if (_IRpoles.empty()) return value_base(mom, kernels, PL);

   // To regulate the diagram in the IR, we map each region with
   // an IR pole at q = qIR != 0 onto coordinates with the pole at q = 0
   double value = 0;
   // need to regulate only the unique IR poles
   // e.g. in the covariance limit, two IR poles can be degenerate
   // and we should treat them simultaneously
   std::vector<ThreeVector> uniqueIRpoles;
   // pole at q = 0
   uniqueIRpoles.push_back(ThreeVector(0, 0, 0));
   // loop over the nonzero poles
   for (auto& pole_prop : _IRpoles) {
      // check if pole is unique
      bool is_unique = true;
      ThreeVector pole = pole_prop.p(mom);
      for (auto unique_pole : uniqueIRpoles) {
         if (pole == unique_pole) {
            is_unique = false;
            break;
         }
      }
      if (is_unique) { uniqueIRpoles.push_back(pole); }
   }
   // now loop over all the unique IR poles
   for (size_t i = 0; i < uniqueIRpoles.size(); i++) {
      // for these poles we change variables: q -> q + pole
      // so that the pole maps to 0 and we exclude all other poles
      ThreeVector pole = uniqueIRpoles[i];
      double PSregion = 1;
      // loop over all other poles and make PS cuts for each
      for (size_t j = 0; j < uniqueIRpoles.size(); j++) {
         if (j != i) {
            ThreeVector pole_j = uniqueIRpoles[j];
            PSregion *= theta(mom[Momentum::q], mom[Momentum::q] + pole - pole_j);
         }
      }
      // copy and shift the diagram momentum for the pole
      LabelMap<Momentum, ThreeVector> mom_shift = mom;
      mom_shift[Momentum::q] = mom[Momentum::q] + pole;
      // add the diagram value for this shifted momentum, times the PS factor
      value += PSregion * value_base(mom_shift, kernels, PL);
   }

   return value;
}

//------------------------------------------------------------------------------
double DiagramOneLoop::value(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const
{
   /*
    * To return the IR regulated diagram symmetrized over external momenta,
    * we symmetrize the IR regulated diagram with the input momentum routing
    * over external momentum configurations.
    * To make the symmetrization more efficient, we compute only the
    * momentum configurations giving distinct diagram values,
    * and multiply each by the appropriate symmetry factor
    */

   double value = 0;
   // loop over external momentum permutations
   // symmetrize over q -> -q
   for (auto perm : _perms) {
      LabelMap<Momentum, ThreeVector> mom_perm = mom;
      mom_perm.permute(perm);
      value += 0.5 * value_base_IRreg(mom_perm, kernels, PL);
      ThreeVector mq = -1 * mom_perm[Momentum::q];
      mom_perm[Momentum::q] = mq;
      value += 0.5 * value_base_IRreg(mom_perm, kernels, PL);
   }

   return value;
}

} // namespace fnfast

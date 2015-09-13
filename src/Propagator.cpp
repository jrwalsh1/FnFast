//------------------------------------------------------------------------------
/// \file Propagator.cpp
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
//    Implementation of classes Propagator
//------------------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "Propagator.hpp"

//------------------------------------------------------------------------------
Propagator::Propagator(MomentumMap<LabelFlow> components)
: _components(components) {}

//------------------------------------------------------------------------------
ThreeVector Propagator::p(MomentumMap<ThreeVector> mom) const
{
   // output container
   ThreeVector pvec;
   // loop and add momentum components
   for (auto const& label : _components.labels()) {
      pvec += static_cast<int>(_components[label]) * mom[label];
   }

   return pvec;
}

//------------------------------------------------------------------------------
vector<MomentumLabel> Propagator::labels() const
{
   vector<MomentumLabel> labels;
   // add label if its coefficient is not kNull
   for (auto const& label : _components.labels()) {
      if (_components[label] != LabelFlow::kNull) {
         labels.push_back(label);
      }
   }
   return labels;
}

//------------------------------------------------------------------------------
bool Propagator::hasLabel(MomentumLabel label) const
{
   return _components.hasLabel(label);
}

//------------------------------------------------------------------------------
bool Propagator::isNull() const
{
   // return true if any label has nonzero coefficient
   for (auto const& label : _components.labels()) {
      if (_components[label] != LabelFlow::kNull) { return false; }
   }

   return true;
}

//------------------------------------------------------------------------------
Propagator Propagator::reverse() const
{
   // copy the underlying MomentumMap, reverse labels
   MomentumMap<LabelFlow> rev_comp = _components;
   for (auto label : rev_comp.labels()) {
      rev_comp[label] = reverse_flow(rev_comp[label]);
   }

   return Propagator(rev_comp);
}

//------------------------------------------------------------------------------
Propagator Propagator::IRpole(MomentumLabel label) const
{
   // set up the map for the new propagator
   unordered_map<MomentumLabel, LabelFlow> pole;
   // first check to make sure the label we're solving for is present in the propagator
   if ( !hasLabel(label) ) {
      cout << "Propagator::IRpole : no component of propagator with given label!" << endl;
      return Propagator(MomentumMap<LabelFlow>(pole));
   }

   // now the case where the label is present
   LabelFlow pole_flow = _components[label];
   // if the component has flow = 0, we also return a null propagator
   if (pole_flow == LabelFlow::kNull) {
      cout << "Propagator::IRpole : component of propagator has flow 0!" << endl;
      return Propagator(MomentumMap<LabelFlow>(pole));
   }
   // otherwise solve for the pole (effectively scale all factors by -1 / fac)
   for (auto const& proplabel : _components.labels()) {
      if (proplabel != label) {
         if (pole_flow == LabelFlow::kMinus) {
            pole[proplabel] = _components[proplabel];
         } else {
            pole[proplabel] = reverse_flow(_components[proplabel]);
         }
      }
   }

   return Propagator(MomentumMap<LabelFlow>(pole));
}

//------------------------------------------------------------------------------
Propagator::LabelFlow Propagator::reverse_flow(Propagator::LabelFlow sign) {
   return (sign == LabelFlow::kMinus) ? LabelFlow::kPlus  :
          (sign == LabelFlow::kPlus)  ? LabelFlow::kMinus :
                                        LabelFlow::kNull;
}

//------------------------------------------------------------------------------
std::ostream& Propagator::print(std::ostream& out) const
{
   stringstream ss;
   ss << "propagator: ";
   for (auto const& label : _components.labels()) {
      string sign = (_components[label] == LabelFlow::kPlus) ? " + k" :
                    (_components[label] == LabelFlow::kMinus) ? " - k" :
                                                               " 0 k";
      ss << sign << static_cast<int>(label);
   }
   return out << ss.str();
}

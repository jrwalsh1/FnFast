//------------------------------------------------------------------------------
/// \file Propagator.cpp
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
//    Implementation of class Propagator
//------------------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "Propagator.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
Propagator::Propagator(LabelMap<Momentum, LabelFlow> components)
: _components(components) {}

//------------------------------------------------------------------------------
ThreeVector Propagator::p(LabelMap<Momentum, ThreeVector> mom) const
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
std::vector<Momentum> Propagator::labels() const
{
   std::vector<Momentum> labels;
   // add label if its coefficient is not Null
   for (auto const& label : _components.labels()) {
      if (_components[label] != LabelFlow::Null) {
         labels.push_back(label);
      }
   }
   return labels;
}

//------------------------------------------------------------------------------
bool Propagator::hasLabel(Momentum label) const
{
   return _components.hasLabel(label);
}

//------------------------------------------------------------------------------
bool Propagator::isNull() const
{
   // return true if any label has nonzero coefficient
   for (auto const& label : _components.labels()) {
      if (_components[label] != LabelFlow::Null) { return false; }
   }

   return true;
}

//------------------------------------------------------------------------------
Propagator Propagator::reverse() const
{
   // copy the underlying Momentum map, reverse labels
   LabelMap<Momentum, LabelFlow> rev_comp = _components;
   for (auto label : rev_comp.labels()) {
      rev_comp[label] = reverse_flow(rev_comp[label]);
   }

   return Propagator(rev_comp);
}

//------------------------------------------------------------------------------
Propagator Propagator::IRpole(Momentum label) const
{
   // set up the map for the new propagator
   std::unordered_map<Momentum, LabelFlow> pole;
   // first check to make sure the label we're solving for is present in the propagator
   if ( !hasLabel(label) ) {
      std::cout << "Propagator::IRpole : no component of propagator with given label!" << std::endl;
      return Propagator(LabelMap<Momentum, LabelFlow>(pole));
   }

   // now the case where the label is present
   LabelFlow pole_flow = _components[label];
   // if the component has flow = 0, we also return a null propagator
   if (pole_flow == LabelFlow::Null) {
      std::cout << "Propagator::IRpole : component of propagator has flow 0!" << std::endl;
      return Propagator(LabelMap<Momentum, LabelFlow>(pole));
   }
   // otherwise solve for the pole (effectively scale all factors by -1 / fac)
   for (auto const& proplabel : _components.labels()) {
      if (proplabel != label) {
         if (pole_flow == LabelFlow::Minus) {
            pole[proplabel] = _components[proplabel];
         } else {
            pole[proplabel] = reverse_flow(_components[proplabel]);
         }
      }
   }

   return Propagator(LabelMap<Momentum, LabelFlow>(pole));
}

//------------------------------------------------------------------------------
Propagator::LabelFlow Propagator::reverse_flow(Propagator::LabelFlow sign) {
   return (sign == LabelFlow::Minus) ? LabelFlow::Plus  :
          (sign == LabelFlow::Plus)  ? LabelFlow::Minus :
                                        LabelFlow::Null;
}

//------------------------------------------------------------------------------
std::ostream& Propagator::print(std::ostream& out) const
{
   std::stringstream ss;
   ss << "propagator: ";
   for (auto const& label : _components.labels()) {
      std::string sign = (_components[label] == LabelFlow::Plus) ? " + k" :
                         (_components[label] == LabelFlow::Minus) ? " - k" :
                                                                     " 0 k";
      ss << sign << static_cast<int>(label);
   }
   return out << ss.str();
}

} // namespace fnfast

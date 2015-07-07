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

#include "Propagator.hpp"

//------------------------------------------------------------------------------
Propagator::Propagator(unordered_map<Momenta::MomentumLabel, double> components)
: _components(components)
{
   // assign any missing keys to 0 using map::insert
   // exisiting keys will not be replaced
   for (size_t c = 0; c < Momenta::momentumlabels.size(); c++) {
      _components.insert(pair<Momenta::MomentumLabel, double>(Momenta::momentumlabels[c], 0));
   }
}

//------------------------------------------------------------------------------
ThreeVector Propagator::p(DiagramMomenta mom)
{
   // output container
   ThreeVector pvec;
   // handle the case where the momentum has no components
   if (_components.size() == 0) {
      return pvec;
   }
   // loop over components, add to p
   for (size_t c = 0; c < mom.labels.size(); c++) {
      double fac = _components[mom.labels[c]];
      pvec += fac * mom[mom.labels[c]];
   }
   
   return pvec;
}

//------------------------------------------------------------------------------
bool Propagator::hasLabel(Momenta::MomentumLabel label)
{
   return (_components[label] != 0);
}

//------------------------------------------------------------------------------
bool Propagator::isNull()
{
   // return true if any label has nonzero coefficient
   for (size_t c = 0; c < Momenta::momentumlabels.size(); c++) {
      if (_components[Momenta::momentumlabels[c]] != 0) { return false; }
   }

   return true;
}

//------------------------------------------------------------------------------
Propagator Propagator::reverse()
{
   // get the map
   unordered_map<Momenta::MomentumLabel, double> comp = _components;
   // multiply all factors (values) by -1
   for (size_t c = 0; c < Momenta::momentumlabels.size(); c++) {
      comp[Momenta::momentumlabels[c]] *= -1;
   }

   // return the propagator
   Propagator prop(comp);
   return prop;
}

//------------------------------------------------------------------------------
Propagator Propagator::IRpole(Momenta::MomentumLabel label)
{
   // first check to make sure the label we're solving for is present in the propagator
   if ( !hasLabel(label) ) {
      cout << "Propagator::IRpole : no component of propagator with given label!" << endl;
      unordered_map<Momenta::MomentumLabel, double> nullcomp;
      Propagator nullprop(nullcomp);
      return nullprop;
   }
   
   // now the case where the label is present
   unordered_map<Momenta::MomentumLabel, double> comp = _components;
   // get the factor for the label, then set the label factor to 0
   int fac = _components[label];
   comp[label] = 0;
   // now scale all factors by -1 / fac
   for (size_t c = 0; c < Momenta::momentumlabels.size(); c++) {
      comp[Momenta::momentumlabels[c]] *= (-1. / fac);
   }
   Propagator prop(comp);
   
   return comp;
}

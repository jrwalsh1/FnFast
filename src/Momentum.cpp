//------------------------------------------------------------------------------
/// \file Momentum.cpp
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
//    Implementation of classes DiagramMomenta
//------------------------------------------------------------------------------

#include "Momentum.hpp"


//------------------------------------------------------------------------------
DiagramMomenta::DiagramMomenta(vector<Momenta::MomentumLabel> labelset)
: labels(labelset)
{
   // assign all keys a null momentum
   ThreeVector pnull;
   for (size_t c = 0; c < labels.size(); c++) {
      _momenta[labels[c]] = pnull;
   }
}

//------------------------------------------------------------------------------
DiagramMomenta::DiagramMomenta(unordered_map<Momenta::MomentumLabel, ThreeVector> momenta)
: _momenta(momenta)
{
   for (auto mom : momenta) {
      labels.push_back(mom.first);
   }
}

//------------------------------------------------------------------------------
void DiagramMomenta::set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> momenta)
{
   _momenta = momenta;
   vector<Momenta::MomentumLabel> inlabels;
   for (auto mom : momenta) {
      inlabels.push_back(mom.first);
   }
   labels = inlabels;
}

//------------------------------------------------------------------------------
void DiagramMomenta::permute(unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> perm)
{
   unordered_map<Momenta::MomentumLabel, ThreeVector> curr_momenta = _momenta;
   for (auto item : perm) {
      _momenta[item.first] = curr_momenta[item.second];
   }
}

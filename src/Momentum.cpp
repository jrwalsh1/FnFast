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
DiagramMomenta::DiagramMomenta(map<Momenta::MomentumLabel, ThreeVector> momenta)
: _momenta(momenta)
{
   // assign any missing keys to 0 using map::insert
   // exisiting keys will not be replaced
   ThreeVector pnull;
   for (size_t c = 0; c < Momenta::momentumlabels.size(); c++) {
      _momenta.insert(pair<Momenta::MomentumLabel, ThreeVector>(Momenta::momentumlabels[c], pnull));
   }
}

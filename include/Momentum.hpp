//------------------------------------------------------------------------------
/// \file Momentum.hpp
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
//    Interface of class DiagramMomenta
//------------------------------------------------------------------------------

#ifndef MOMENTUM_HPP
#define MOMENTUM_HPP

#include <vector>
#include <unordered_map>

#include "Labels.hpp"
#include "ThreeVector.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class DiagramMomenta
 *
 * \brief class for momenta of loop, external lines
 *
 * DiagramMomenta(map<MomentumLabel, ThreeVector> momenta)
 *
 * Container for the loop and external momenta.
 * Instantiated with a map from MomentumLabel to ThreeVector
 *
 * Provides functions for:
 * - Accessing momenta
 * - Modifying momenta
 */
//------------------------------------------------------------------------------

class DiagramMomenta
{
   public:
      vector<Momenta::MomentumLabel> labels;                               ///< momentum labels

   private:
      unordered_map<Momenta::MomentumLabel, ThreeVector> _momenta;         ///< ThreeVector for a given MomentumLabel

   public:
      /// constructors
      DiagramMomenta() {}
      DiagramMomenta(vector<Momenta::MomentumLabel> labelset);
      DiagramMomenta(unordered_map<Momenta::MomentumLabel, ThreeVector> momenta);
      /// destructor
      virtual ~DiagramMomenta() {}

      /// set the momenta
      void set_momenta(unordered_map<Momenta::MomentumLabel, ThreeVector> momenta);

      /// accessors
      unordered_map<Momenta::MomentumLabel, ThreeVector> momenta() { return _momenta; }
      ThreeVector& operator[](Momenta::MomentumLabel label) { return _momenta[label]; }

      /// permute momenta using a map
      void permute(unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> perm);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // MOMENTUM_HPP

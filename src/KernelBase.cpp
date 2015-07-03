//------------------------------------------------------------------------------
/// \file KernelBase.cpp
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
//    Implementation of non-virtual functions in KernelBase
//------------------------------------------------------------------------------

#include <algorithm>

#include "ThreeVector.hpp"
#include "KernelBase.hpp"

//------------------------------------------------------------------------------
double KernelBase::Fn_sym(vector<ThreeVector> p)
{
   double value = 0;
   int nperm = 0; // count the permutations
   // use the next_permutation algorithm together with the comparison operator
   // in ThreeVector to generate permutations
   sort(p.begin(),p.end());
   do {
      nperm++;
      value += Fn(p);
   } while (next_permutation(p.begin(), p.end()));

   return value / nperm;
}

//------------------------------------------------------------------------------
double KernelBase::Gn_sym(vector<ThreeVector> p)
{
   double value = 0;
   int nperm = 0; // count the permutations
   // use the next_permutation algorithm together with the comparison operator
   // in ThreeVector to generate permutations
   sort(p.begin(),p.end());
   do {
      nperm++;
      value += Gn(p);
   } while (next_permutation(p.begin(), p.end()));

   return value / nperm;
}

//------------------------------------------------------------------------------
/// \file SPTkernels.cpp
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
//    Implementation of functions in SPTkernels
//------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>

#include "ThreeVector.hpp"
#include "SPTkernels.hpp"

//------------------------------------------------------------------------------
double SPTkernels::cF_alpha(int n)
{
   return (2*n + 1.) / ((n - 1) * (2 * n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::cF_beta(int n)
{
   return 2. / ((n - 1) * (2*n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::cG_alpha(int n)
{
   return 3. / ((n - 1) * (2*n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::cG_beta(int n)
{
   return (2. * n) / ((n - 1) * (2*n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::alpha(ThreeVector p1, ThreeVector p2)
{
   // handle the IR limit with an explicit cutoff
   double eps = 1e-12;
   if (p1*p1 < eps) { return 0; }

   // otherwise
   return ((p1 + p2)*p1) / (p1*p1);
}

//------------------------------------------------------------------------------
double SPTkernels::beta(ThreeVector p1, ThreeVector p2)
{
   // handle the IR limit with an explicit cutoff
   double eps = 1e-12;
   if ((p1*p1 < eps) || (p2*p2 < eps)) { return 0; }
   
   // otherwise
   return ((p1 + p2)*(p1 + p2)) * (p1*p2) / (2 * (p1*p1) * (p2*p2));
}

//------------------------------------------------------------------------------
double SPTkernels::Fn(vector<ThreeVector> p)
{
   int n = p.size();
   // handle the nonsense and trivial case
   if (n == 0) { return 0; }
   if (n == 1) { return 1; }
   
   // handle the root case for the recursion
   if (n == 2) { return cF_alpha(n) * alpha(p[0], p[1]) + cF_beta(n) * beta(p[0], p[1]); }
   
   // now the nontrivial recursion cases
   double Fnval = 0;
   for (int k = 1; k < n; k++) {
      vector<ThreeVector> p_k(p.begin(), p.begin() + k);
      vector<ThreeVector> p_nk(p.begin() + k, p.end());
      // sum the subsets of momenta (start with p0 = 0)
      ThreeVector p0;
      ThreeVector pktot = accumulate(p_k.begin(), p_k.end(), p0);
      ThreeVector pnktot = accumulate(p_nk.begin(), p_nk.end(), p0);
      // add term to the sum
      Fnval += Gn(p_k) * (cF_alpha(n) * alpha(pktot, pnktot) * Fn(p_nk) + cF_beta(n) * beta(pktot, pnktot) * Gn(p_nk));
   }
   return Fnval;
}

//------------------------------------------------------------------------------
double SPTkernels::Gn(vector<ThreeVector> p)
{
   int n = p.size();
   // handle the nonsense and trivial case
   if (n == 0) { return 0; }
   if (n == 1) { return 1; }
   
   // handle the root case for the recursion
   if (n == 2) { return cG_alpha(n) * alpha(p[0], p[1]) + cG_beta(n) * beta(p[0], p[1]); }
   
   // now the nontrivial recursion cases
   double Gnval = 0;
   for (int k = 1; k < n; k++) {
      vector<ThreeVector> p_k(p.begin(), p.begin() + k);
      vector<ThreeVector> p_nk(p.begin() + k, p.end());
      // sum the subsets of momenta (start with p0 = 0)
      ThreeVector p0;
      ThreeVector pktot = accumulate(p_k.begin(), p_k.end(), p0);
      ThreeVector pnktot = accumulate(p_nk.begin(), p_nk.end(), p0);
      // add term to the sum
      Gnval += Gn(p_k) * (cG_alpha(n) * alpha(pktot, pnktot) * Fn(p_nk) + cG_beta(n) * beta(pktot, pnktot) * Gn(p_nk));
   }
   return Gnval;
}

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

namespace fnfast {

//------------------------------------------------------------------------------
SPTkernels::SPTkernels()
{
   // precompute alpha, beta coefficients up to F7: this will be all that's ever needed
   for (int c = 0; c <= 7; c++) {
      _cFalpha[c] = cF_alpha(c);
      _cFbeta[c] = cF_beta(c);
      _cGalpha[c] = cG_alpha(c);
      _cGbeta[c] = cG_beta(c);
}
   for (int i = 0; i < 7; i++) {
      for (int j = 0; j < 7; j++) {
         _alpha[i][j] = 0;
         _beta[i][j] = 0;
      }
   }
}

//------------------------------------------------------------------------------
double SPTkernels::cF_alpha(int n)
{
   if (n < 2) { return 0; }
   return (2*n + 1.) / ((n - 1) * (2 * n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::cF_beta(int n)
{
   if (n < 2) { return 0; }
   return 2. / ((n - 1) * (2*n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::cG_alpha(int n)
{
   if (n < 2) { return 0; }
   return 3. / ((n - 1) * (2*n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::cG_beta(int n)
{
   if (n < 2) { return 0; }
   return (2. * n) / ((n - 1) * (2*n + 3));
}

//------------------------------------------------------------------------------
double SPTkernels::alpha(ThreeVector& p1, ThreeVector& p2)
{
   // handle the IR limit with an explicit cutoff
   double eps = 1e-12;
   if (p1*p1 < eps) { return 0; }

   // otherwise
   return ((p1 + p2)*p1) / (p1*p1);
}

//------------------------------------------------------------------------------
double SPTkernels::beta(ThreeVector& p1, ThreeVector& p2)
{
   // handle the IR limit with an explicit cutoff
   double eps = 1e-12;
   if ((p1*p1 < eps) || (p2*p2 < eps)) { return 0; }

   // otherwise
   return ((p1 + p2)*(p1 + p2)) * (p1*p2) / (2 * (p1*p1) * (p2*p2));
}

//------------------------------------------------------------------------------
double SPTkernels::Fn_sym(std::vector<ThreeVector>& p)
{
   size_t n = p.size();
   // handle the low multiplicity cases
   if (n == 0) { return 0; }
   if (n == 1) { return 1; }
   if (n == 2) { return 0.5 * (_cFalpha[n] * (alpha(p[0], p[1]) + alpha(p[1], p[0])) + _cFbeta[n] * 2 * beta(p[0], p[1])); } // note beta symmetric

   // stash the alpha and beta function values
   // and set up the vector of index labels
   size_t np = p.size();
   std::vector<int> indices;
   for (int i = 0; i < np; i++) {
      indices.push_back(i);
      for (size_t j = 0; j < i; j++) {
         _alpha[i][j] = alpha(p[i], p[j]);
         _alpha[j][i] = alpha(p[j], p[i]);
         _beta[i][j] = beta(p[i], p[j]);
         _beta[j][i] = _beta[i][j]; // beta is symmetric!
      }
   }

   double value = 0;
   int nperm = 0; // count the permutations
   // use the next_permutation algorithm to loop through the permutations
   std::sort(indices.begin(), indices.end());
   do {
      nperm++;
      std::vector<ThreeVector> pperm;
      for (size_t c = 0; c < np; c++) {
         pperm.push_back(p[indices[c]]);
      }
      value += Fn(pperm, indices);
   } while (std::next_permutation(indices.begin(), indices.end()));

   return value / nperm;
}

//------------------------------------------------------------------------------
double SPTkernels::Gn_sym(std::vector<ThreeVector>& p)
{
   size_t n = p.size();
   // handle the low multiplicity cases
   if (n == 0) { return 0; }
   if (n == 1) { return 1; }
   if (n == 2) { return 0.5 * (_cGalpha[n] * (alpha(p[0], p[1]) + alpha(p[1], p[0])) + _cGbeta[n] * 2 * beta(p[0], p[1])); } // note beta symmetric

   // stash the alpha and beta function values
   // and set up the vector of index labels
   size_t np = p.size();
   std::vector<int> indices;
   for (int i = 0; i < np; i++) {
      indices[i] = i;
      for (size_t j = 0; j < i; j++) {
         _alpha[i][j] = alpha(p[i], p[j]);
         _alpha[j][i] = alpha(p[j], p[i]);
         _beta[i][j] = beta(p[i], p[j]);
         _beta[j][i] = _beta[i][j]; // beta is symmetric!
      }
   }

   double value = 0;
   int nperm = 0; // count the permutations
   // use the next_permutation algorithm to loop through the permutations
   std::sort(indices.begin(), indices.end());
   do {
      nperm++;
      value += Gn(p, indices);
   } while (std::next_permutation(indices.begin(), indices.end()));

   return value / nperm;
}

//------------------------------------------------------------------------------
double SPTkernels::Fn(std::vector<ThreeVector>& p, std::vector<int>& indices)
{
   size_t n = indices.size();
   // handle the nonsense and trivial case
   if (n == 0) { return 0; }
   if (n == 1) { return 1; }

   // handle the root case for the recursion
   if (n == 2) { return _cFalpha[n] * _alpha[indices[0]][indices[1]] + _cFbeta[n] * _beta[indices[0]][indices[1]]; }

   // now the nontrivial recursion cases
   double Fnval = 0;
   for (int k = 1; k < n; k++) {
      // split the vector list for the recursion
      std::vector<ThreeVector> p_k(p.begin(), p.begin() + k);
      std::vector<ThreeVector> p_nk(p.begin() + k, p.end());
      // split the index list for the recursion
      std::vector<int> indices_k(indices.begin(), indices.begin() + k);
      std::vector<int> indices_nk(indices.begin() + k, indices.end());
      // sum the subsets of momenta (start with p0 = 0)
      ThreeVector p0;
      ThreeVector pktot = std::accumulate(p_k.begin(), p_k.end(), p0);
      ThreeVector pnktot = std::accumulate(p_nk.begin(), p_nk.end(), p0);
      // add term to the sum
      Fnval += Gn(p_k, indices_k) * (_cFalpha[n] * alpha(pktot, pnktot) * Fn(p_nk, indices_nk) + _cFbeta[n] * beta(pktot, pnktot) * Gn(p_nk, indices_nk));
   }
   return Fnval;
}


//------------------------------------------------------------------------------
double SPTkernels::Gn(std::vector<ThreeVector>& p, std::vector<int>& indices)
{
   size_t n = indices.size();
   // handle the nonsense and trivial case
   if (n == 0) { return 0; }
   if (n == 1) { return 1; }

   // handle the root case for the recursion
   if (n == 2) { return _cGalpha[n] * _alpha[indices[0]][indices[1]] + _cGbeta[n] * _beta[indices[0]][indices[1]]; }

   // now the nontrivial recursion cases
   double Gnval = 0;
   for (int k = 1; k < n; k++) {
      // split the vector list for the recursion
      std::vector<ThreeVector> p_k(p.begin(), p.begin() + k);
      std::vector<ThreeVector> p_nk(p.begin() + k, p.end());
      // split the index list for the recursion
      std::vector<int> indices_k(indices.begin(), indices.begin() + k);
      std::vector<int> indices_nk(indices.begin() + k, indices.end());
      // sum the subsets of momenta (start with p0 = 0)
      ThreeVector p0;
      ThreeVector pktot = std::accumulate(p_k.begin(), p_k.end(), p0);
      ThreeVector pnktot = std::accumulate(p_nk.begin(), p_nk.end(), p0);
      // add term to the sum
      Gnval += Gn(p_k, indices_k) * (_cGalpha[n] * alpha(pktot, pnktot) * Fn(p_nk, indices_nk) + _cGbeta[n] * beta(pktot, pnktot) * Gn(p_nk, indices_nk));
   }
   return Gnval;
}

} // namespace fnfast

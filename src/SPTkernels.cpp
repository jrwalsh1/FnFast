//------------------------------------------------------------------------------
/// \file SPTkernels.cpp
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
//    Implementation of functions in SPTkernels
//------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <utility>

#include "ThreeVector.hpp"
#include "SPTkernels.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
SPTkernels::SPTkernels()
{
   // precompute useful quantities for fast evaluation up to n = 7
   // higher n should never be needed, but easy to add if desired
   // ----------------------------------------
   // precompute alpha, beta coefficients
   for (int c = 0; c <= 7; c++) {
      _cFalpha[c] = cF_alpha(c);
      _cFbeta[c] = cF_beta(c);
      _cGalpha[c] = cG_alpha(c);
      _cGbeta[c] = cG_beta(c);
   }

   // compute binomial coefficients
   _binom = std::vector<std::vector<int> > (8);
   for (int i = 0; i < 8; i++) {
      // default fill all entries with 0
      _binom[i] = std::vector<int> (i+1, 0);
      for (int j = 0; j <= i; j++) {
         // fill with actual values
         _binom[i][j] = fact(i) / (fact(i - j) * fact(j));
      }
   }

   // compute permutation subsets, cache
   _permset = std::vector<std::vector<std::vector<std::vector<int> > > > (8);
   // i counts the total number of indices
   for (int i = 1; i < 8; i++) {
      // default fill each top level
      _permset[i] = std::vector<std::vector<std::vector<int> > > (i+1);
      for (int j = 1; j <= i; j++) {
         // fill with actual values
         _permset[i][j] = generate_permset(j, i);
      }
   }

   // compute subset pairs, cache
   _subsetpairs = std::vector<std::vector<std::vector<std::vector<SubsetPair> > > > (8);
   // i indexes the number of total indices
   for (int i = 1; i < 8; i++) {
      // default fill the top level
      _subsetpairs[i] = std::vector<std::vector<std::vector<SubsetPair> > > (i+1);
      // j counts the number of indices in the subset
      for (int j = 1; j <= i; j++) {
         // container for all subsets of the given subset from j, i
         std::vector<std::vector<SubsetPair> > subsetpairs_ij;
         // generate all subsets of j indices built from i indices
         std::vector<std::vector<int> > perms_ij = generate_permset(j, i);
         // for each index subset, find all subset pairs of that subset
         // remember we need the hashing function from i indices
         for (auto perm : perms_ij) {
            std::vector<SubsetPair> perm_pairs = generate_pairedsubsets(perm, i);
            subsetpairs_ij.push_back(perm_pairs);
         }
         _subsetpairs[i][j] = subsetpairs_ij;
      }
   }

   // fill the Fn_sym and Gn_sym arrays
   _Fn_sym = std::vector<std::vector<std::vector<double> > > (8);
   // n is the total number of indices (n for the kernel we're calculating)
   for (int n = 1; n < 8; n++) {
      // default fill the top level
      _Fn_sym[n] = std::vector<std::vector<double> > (n+1);
      // k is the number of indices in the daughter kernel calculations
      // the ones in the recursion relation
      for (int k = 1; k <= n; k++) {
         _Fn_sym[n][k] = std::vector<double> (_permset[n][k].size(), 0);
      }
   }
   _Gn_sym = _Fn_sym;
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
double SPTkernels::alpha(const ThreeVector& p1, const ThreeVector& p2)
{
   // handle the IR limit with an explicit cutoff
   double eps = 1e-12;
   if (p1*p1 < eps) { return 0; }

   // otherwise
   return 1 + (p2*p1) / (p1*p1);
}

//------------------------------------------------------------------------------
double SPTkernels::beta(const ThreeVector& p1, const ThreeVector& p2)
{
   // handle the IR limit with an explicit cutoff
   double eps = 1e-12;
   if ((p1*p1 < eps) || (p2*p2 < eps)) { return 0; }

   // otherwise
   return ((p1 + p2)*(p1 + p2)) * (p1*p2) / (2 * (p1*p1) * (p2*p2));
}

//------------------------------------------------------------------------------
double SPTkernels::Fn_sym(const std::vector<ThreeVector>& p)
{
   // calculates Fn_sym by calculating all lower multiplicity cases first
   // and using those results to build the requested case
   int n = p.size();

   // need to do this in a specific order so that the higher multiplicity
   // functions can use the lower multiplicity function results
   // k represents the number of momenta in the set,
   // e.g. for k = 3 we calculate F_3^sym (p1, p2, p3) + F_3^sym (p1, p2, p4) + ...
   for (int k = 1; k < n; k++) {
      // for every momentum permutation of k elements from the n,
      // calculate the symmetrized Fn and Gn
      for (int j = 0; j < _permset[n][k].size(); j++) {
         _Fn_sym[n][k][j] = Fn_sym_build(p, _permset[n][k][j], j);
         _Gn_sym[n][k][j] = Gn_sym_build(p, _permset[n][k][j], j);
      }
   }
   // now calculate the main result
   double result = Fn_sym_build(p, _permset[n][n][0], 0);
   _Fn_sym[n][n][0] = result;

   return result;
}

//------------------------------------------------------------------------------
double SPTkernels::Gn_sym(const std::vector<ThreeVector>& p)
{
   // calculates Gn_sym by calculating all lower multiplicity cases first
   // and using those results to build the requested case
   int n = p.size();

   // need to do this in a specific order so that the higher multiplicity
   // functions can use the lower multiplicity function results
   // k represents the number of momenta in the set,
   // e.g. for k = 3 we calculate F_3^sym (p1, p2, p3) + F_3^sym (p1, p2, p4) + ...
   for (int k = 1; k < n; k++) {
      // for every momentum permutation of k elements from the n,
      // calculate the symmetrized Fn and Gn
      for (int j = 0; j < _permset[n][k].size(); j++) {
         _Fn_sym[n][k][j] = Fn_sym_build(p, _permset[n][k][j], j);
         _Gn_sym[n][k][j] = Gn_sym_build(p, _permset[n][k][j], j);
      }
   }
   // now calculate the main result
   double result = Gn_sym_build(p, _permset[n][n][0], 0);
   _Gn_sym[n][n][0] = result;

   return result;
}

//------------------------------------------------------------------------------
double SPTkernels::Fn_sym_build(const std::vector<ThreeVector>& p, const std::vector<int> indices, int hashvalue)
{
   // calculates Fn_sym using the results of lower multiplicity calculations
   int n = p.size();
   int k = indices.size();

   // handle the base cases
   if (k == 1) { return 1; }
   if (k == 2) { return 0.5 * (_cFalpha[2] * (alpha(p[indices[0] - 1], p[indices[1] - 1])
                                             + alpha(p[indices[1] - 1], p[indices[0] - 1]))
                              + _cFbeta[2] * 2 * beta(p[indices[0] - 1], p[indices[1] - 1]));
               } // note beta symmetric

   // now do the recursion case
   double result = 0;
   // need to sum over all subsets of the given index set
   for (const auto& subsetpair : _subsetpairs[n][k][hashvalue]) {
      int nA = subsetpair.subsetA.size();
      int nB = subsetpair.subsetB.size();
      ThreeVector pA;
      for (const auto& index : subsetpair.subsetA ) { pA += p[index - 1]; }
      ThreeVector pB;
      for (const auto& index : subsetpair.subsetB ) { pB += p[index - 1]; }
      double combfac = 1. / _binom[k][nA];
      // atomic quantities
      double FnA = _Fn_sym[n][nA][subsetpair.hashA];
      double FnB = _Fn_sym[n][nB][subsetpair.hashB];
      double GnA = _Gn_sym[n][nA][subsetpair.hashA];
      double GnB = _Gn_sym[n][nB][subsetpair.hashB];
      double alphaAB = alpha(pA, pB);
      double alphaBA = alpha(pB, pA);
      double betaval = beta(pA, pB); // note beta symmetric
      // add subset result
      result += combfac * GnA * (_cFalpha[k] * alphaAB * FnB + _cFbeta[k] * betaval * GnB);
      result += combfac * GnB * (_cFalpha[k] * alphaBA * FnA + _cFbeta[k] * betaval * GnA);
   }

   return result;
}

//------------------------------------------------------------------------------
double SPTkernels::Gn_sym_build(const std::vector<ThreeVector>& p, const std::vector<int> indices, int hashvalue)
{
   // calculates Fn_sym using the results of lower multiplicity calculations
   int n = p.size();
   int k = indices.size();

   // handle the base cases
   if (k == 1) { return 1; }
   if (k == 2) { return 0.5 * (_cGalpha[2] * (alpha(p[indices[0] - 1], p[indices[1] - 1])
                                             + alpha(p[indices[1] - 1], p[indices[0] - 1]))
                              + _cGbeta[2] * 2 * beta(p[indices[0] - 1], p[indices[1] - 1]));
               } // note beta symmetric

   // now do the recursion case
   double result = 0;
   // need to sum over all subsets of the given index set
   for (const auto& subsetpair : _subsetpairs[n][k][hashvalue]) {
      int nA = subsetpair.subsetA.size();
      int nB = subsetpair.subsetB.size();
      ThreeVector pA;
      for (const auto& index : subsetpair.subsetA ) { pA += p[index - 1]; }
      ThreeVector pB;
      for (const auto& index : subsetpair.subsetB ) { pB += p[index - 1]; }
      double combfac = 1. / _binom[k][nA];
      // atomic quantities
      double FnA = _Fn_sym[n][nA][subsetpair.hashA];
      double FnB = _Fn_sym[n][nB][subsetpair.hashB];
      double GnA = _Gn_sym[n][nA][subsetpair.hashA];
      double GnB = _Gn_sym[n][nB][subsetpair.hashB];
      double alphaAB = alpha(pA, pB);
      double alphaBA = alpha(pB, pA);
      double betaval = beta(pA, pB); // note beta symmetric
      // add subset result
      result += combfac * GnA * (_cGalpha[k] * alphaAB * FnB + _cGbeta[k] * betaval * GnB);
      result += combfac * GnB * (_cGalpha[k] * alphaBA * FnA + _cGbeta[k] * betaval * GnA);
   }

   return result;
}

//------------------------------------------------------------------------------
std::vector<std::vector<int> > SPTkernels::generate_permset(int k, int n)
{
   std::vector<std::vector<int> > permset;
   int index = 1;
   int maxindex = n - k + 1;
   // seed the permutation set with the first index
   for (int i = index; i <= maxindex; i++) {
      permset.push_back(std::vector<int> {i});
   }
   index++;
   maxindex++;
   // loop over spots in the permset up to k
   while (index <= k) {
      std::vector<std::vector<int> > permset_iter;
      // loop over existing vectors
      for (auto elem : permset) {
         int minindex = elem.back() + 1;
         for (int i = minindex; i <= maxindex; i++) {
            permset_iter.push_back(elem);
            permset_iter.back().push_back(i);
         }
      }
      // move onto the next spot
      permset = std::move(permset_iter);
      index++;
      maxindex++;
   }
   return permset;
}

//------------------------------------------------------------------------------
std::vector<SPTkernels::SubsetPair> SPTkernels::generate_pairedsubsets(const std::vector<int>& indices, int nmax)
{
   int n = indices.size();
   // output container
   std::vector<SubsetPair> subsets;
   // there are 2^{n-1} - 1 total pairs of subsets
   // since each index is in one set or the other (and we ignore the full/empty case)
   int totpairs = int_pow(2, n-1) - 1;
   subsets.reserve(totpairs);
   // construct the subsets explicity using the binary representation of
   for (int i = 1; i <= totpairs; i++) {
      std::vector<int> subsetA;
      std::vector<int> subsetB;
      for (int j = 0; j < n; j++) {
         // look at the j^th bit of the pair number
         // 1: subset A, 0: subset B
         if ((i & (1 << j)) != 0) {
            subsetA.push_back(indices[j]);
         } else {
            subsetB.push_back(indices[j]);
         }
      }
      // hash the subsets
      int hashA = hash_perm(subsetA, nmax);
      int hashB = hash_perm(subsetB, nmax);
      SubsetPair subsetpair(subsetA, hashA, subsetB, hashB);
      subsets.push_back(subsetpair);
   }
   return subsets;
}

//------------------------------------------------------------------------------
int SPTkernels::hash_perm(const std::vector<int>& indices, int nmax)
{
   int nslots = indices.size();
   // trivial cases
   if (nslots == 0) { return 0; }
   if (nslots == 1) { return indices[0] - 1; }
   // function to find the number of the index subset in a lexicographic ordering
   // serves as a monotone minimal perfect hash
   // example: 2-index subset out of 4 indices
   // order: [1,2], [1,3], [1,4], [2,3], [2,4], [3,4]
   int navail = nmax;
   int minindex = 1;
   int pos = 0;
   int value = 0;
   // step through the index list
   while (pos < indices.size()) {
      // compute the number of skips
      int nskip = indices[pos] - minindex;
      value += _binom[navail][nslots] - _binom[navail - nskip][nslots];
      // step forward in the vector
      minindex = indices[pos] + 1;
      navail = nmax - minindex + 1;
      nslots--;
      pos++;
   }
   return value;
}


} // namespace fnfast

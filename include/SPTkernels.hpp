//------------------------------------------------------------------------------
/// \file SPTkernels.hpp
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
//    Definition of SPTkernels
//------------------------------------------------------------------------------

#ifndef SPT_KERNELS_HPP
#define SPT_KERNELS_HPP

#include <numeric>

#include "KernelBase.hpp"

namespace fnfast {

class ThreeVector;

//------------------------------------------------------------------------------
/**
 * \namespace SPTkernels
 *
 * \brief Defines the SPT kernels with fast evaluation.
 *
 * Defines the base functions and recursion relations for the SPT kernels
 * Uses fast evaluation methods
 */
//------------------------------------------------------------------------------
class SPTkernels : public KernelBase
{
   private:
      struct SubsetPair {
         std::vector<int> subsetA;
         int hashA;
         std::vector<int> subsetB;
         int hashB;

         SubsetPair(std::vector<int> sA, int hA, std::vector<int> sB, int hB)
         : subsetA(sA), hashA(hA), subsetB(sB), hashB(hB) {}
      };

      double _cFalpha[8];                                   ///< pre-computed constants for Fn coefficient of alpha term
      double _cFbeta[8];                                    ///< pre-computed constants for Fn coefficient of beta term
      double _cGalpha[8];                                   ///< pre-computed constants for Gn coefficient of alpha term
      double _cGbeta[8];                                    ///< pre-computed constants for Gn coefficient of beta term
      std::vector<std::vector<int> > _binom;                ///< binomial coefficients
      std::vector<std::vector<std::vector<std::vector<int> > > > _permset;    ///< set of all perms of j indices from 1..i, for 1 <= j <= i <= 7
      std::vector<std::vector<std::vector<std::vector<SubsetPair> > > > _subsetpairs;   ///< all subset pairs for any subsets of 1..7
      std::vector<std::vector<std::vector<double> > > _Fn_sym;    ///< container for the coefficients
      std::vector<std::vector<std::vector<double> > > _Gn_sym;    ///< container for the coefficients

   public:
      /// constructor
      SPTkernels();
      /// destructor
      ~SPTkernels() {}

      double cF_alpha(int n);   ///< constant for Fn coefficient of alpha term
      double cF_beta(int n);    ///< constant for Fn coefficient of beta term
      double cG_alpha(int n);   ///< constant for Gn coefficient of alpha term
      double cG_beta(int n);    ///< constant for Gn coefficient of beta term

      double alpha(const ThreeVector& p1, const ThreeVector& p2);       ///< kernel function alpha
      double beta(const ThreeVector& p1, const ThreeVector& p2);        ///< kernel function alpha

      double Fn_sym(const std::vector<ThreeVector>& p);    ///< symmetrized SPT kernel Fn (q1, ..., qn)
      double Gn_sym(const std::vector<ThreeVector>& p);    ///< symmetrized SPT kernel Gn (q1, ..., qn)

   private:
      double Fn_sym_build(const std::vector<ThreeVector>& p, const std::vector<int> indices, int hashvalue);    ///< symmetrized SPT kernel Fn (q1, ..., qn), uses precomputed results to calculate
      double Gn_sym_build(const std::vector<ThreeVector>& p, const std::vector<int> indices, int hashvalue);    ///< symmetrized SPT kernel Gn (q1, ..., qn), uses precomputed results to calculate

      std::vector<std::vector<int> > generate_permset(int k, int n);     ///< function to generate all ordered k-element subsets of indices 1..n
      std::vector<SubsetPair> generate_pairedsubsets(const std::vector<int>& indices, int nmax);   ///< function to generate all paired subsets of a set of indices
      int hash_perm(const std::vector<int>& indices, int nmax);     ///< hashing function for

      /// factorial
      static int fact(int n) { return (n == 0 || n == 1) ? 1 : n * fact(n-1); }
      static int int_pow(int base, int exp);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline int SPTkernels::int_pow(int base, int exp)
{
   if (exp == 0) { return 1; }
   if (exp == 1) { return base; }
   if (base == 1) { return 1; }
   if (exp < 0) { return 0; }

   int result = 1;
   while (exp)
   {
      if (exp & 1)
           result *= base;
      exp >>= 1;
      base *= base;
   }

   return result;
}


} // namespace fnfast

#endif // SPT_KERNELS_HPP

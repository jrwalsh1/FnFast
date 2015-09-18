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
 * \brief Defines the SPT kernels.
 *
 * Defines the base functions and recursion relations for the SPT kernels
 * NOTE: the containers in this class are only set up to hold data for up to n = 7
 * But you should never need to compute higher than this!
 */
//------------------------------------------------------------------------------
class SPTkernels : public KernelBase
{
   private:
      double _cFalpha[8];                                   ///< pre-computed constants for Fn coefficient of alpha term
      double _cFbeta[8];                                    ///< pre-computed constants for Fn coefficient of beta term
      double _cGalpha[8];                                   ///< pre-computed constants for Gn coefficient of alpha term
      double _cGbeta[8];                                    ///< pre-computed constants for Gn coefficient of beta term
      double _alpha[7][7];                                  ///< stashed coefficients alpha for base recursion cases
      double _beta[7][7];                                   ///< stashed coefficients beta for base recursion cases

   public:
      /// constructor
      SPTkernels();
      /// destructor
      ~SPTkernels() {}

      double cF_alpha(int n);                               ///< constant for Fn coefficient of alpha term
      double cF_beta(int n);                                ///< constant for Fn coefficient of beta term
      double cG_alpha(int n);                               ///< constant for Gn coefficient of alpha term
      double cG_beta(int n);                                ///< constant for Gn coefficient of beta term

      double alpha(ThreeVector& p1, ThreeVector& p2);       ///< kernel function alpha
      double beta(ThreeVector& p1, ThreeVector& p2);        ///< kernel function alpha

      double Fn_sym(std::vector<ThreeVector>& p);                ///< symmetrized SPT kernel Fn (q1, ..., qn)
      double Gn_sym(std::vector<ThreeVector>& p);                ///< symmetrized SPT kernel Gn (q1, ..., qn)

   private:
      double Fn(std::vector<ThreeVector>& p, std::vector<int>& indices);       ///< SPT kernel Fn, in terms of indices
      double Gn(std::vector<ThreeVector>& p, std::vector<int>& indices);       ///< SPT kernel Gn, in terms of indices
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // SPT_KERNELS_HPP

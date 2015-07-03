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

using namespace std;

class ThreeVector;

//------------------------------------------------------------------------------
/**
 * \namespace SPTkernels
 *
 * \brief Defines the SPT kernels.
 *
 * Defines the base functions and recursion relations for the SPT kernels
 */
//------------------------------------------------------------------------------
class SPTkernels : public KernelBase
{
   public:
      double cF_alpha(int n);                            ///< constant for Fn coefficient of alpha term
      double cF_beta(int n);                             ///< constant for Fn coefficient of beta term
      double cG_alpha(int n);                            ///< constant for Gn coefficient of alpha term
      double cG_beta(int n);                             ///< constant for Gn coefficient of beta term

      double alpha(ThreeVector p1, ThreeVector p2);      ///< kernel function alpha
      double beta(ThreeVector p1, ThreeVector p2);       ///< kernel function alpha

      double Fn(vector<ThreeVector> p);                  ///< SPT kernel Fn (q1, ..., qn)
      double Gn(vector<ThreeVector> p);                  ///< SPT kernel Gn (q1, ..., qn)
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // SPT_KERNELS_HPP

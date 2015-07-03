//------------------------------------------------------------------------------
/// \file KernelBase.hpp
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
//    Definition of base class KernelBase
//------------------------------------------------------------------------------

#ifndef KERNEL_BASE_HPP
#define KERNEL_BASE_HPP

#include <vector>

#include "ThreeVector.hpp"

using namespace std;

//class ThreeVector;

//------------------------------------------------------------------------------
/**
 * \class KernelBase
 *
 * \brief Defines the base class for recusive kernels (SPT, EFTofLSS, etc.).
 *
 * Defines the interface functions shared between kernels.
 * The symmetrization over kernels is defined in the base class.
 */
//------------------------------------------------------------------------------
class KernelBase
{
   public:
      virtual double alpha(ThreeVector p1, ThreeVector p2) = 0;      ///< kernel function alpha
      virtual double beta(ThreeVector p1, ThreeVector p2) = 0;       ///< kernel function alpha

      virtual double Fn(vector<ThreeVector> p) = 0;                  ///< kernel Fn (q1, ..., qn)
      virtual double Gn(vector<ThreeVector> p) = 0;                  ///< kernel Gn (q1, ..., qn)
      
      double Fn_sym(vector<ThreeVector> p);                          ///< symmetrized kernel Fn
      double Gn_sym(vector<ThreeVector> p);                          ///< symmetrized kernel Gn
};

#endif // KERNEL_BASE_HPP

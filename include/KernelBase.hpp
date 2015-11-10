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

namespace fnfast {

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
      virtual double Fn_sym(const std::vector<ThreeVector>& p) = 0;    ///< symmetrized kernel Fn
      virtual double Gn_sym(const std::vector<ThreeVector>& p) = 0;    ///< symmetrized kernel Gn
};

} // namespace fnfast

#endif // KERNEL_BASE_HPP

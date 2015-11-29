//------------------------------------------------------------------------------
/// \file WindowFunctionBase.hpp
//
// Author(s):
//    Daniele Bertolini
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
//    Interface of class WindowFunctionBase
//------------------------------------------------------------------------------

#ifndef WINDOW_FUNCTION_BASE_HPP
#define WINDOW_FUNCTION_BASE_HPP

#include "ThreeVector.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
/**
 * \class WindowFunctionBase
 *
 * \brief Base class for Window function
 *
 * Provides virtual functions:
 * - to evaluate the window function
 */
//------------------------------------------------------------------------------

class WindowFunctionBase
{
   public:
      /// returns the linear power spectrum
      virtual double operator()(ThreeVector x) = 0;
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace std

#endif // WINDOW_FUNCTION_BASE_HPP

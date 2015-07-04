//------------------------------------------------------------------------------
/// \file WindowFunctionTopHat.hpp
//
// Author(s):
//    Daniele Bertolini
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
//    Interface of class WindowFunctionTopHat
//------------------------------------------------------------------------------

#ifndef WINDOW_FUNCTION_TOPHAT_HPP
#define WINDOW_FUNCTION_TOPHAT_HPP

#include "ThreeVector.hpp"
#include <cmath>

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class WindowFunctionTopHat
 *
 * \brief Base class for Window function
 *
 * Provides functions:
 * - to evaluate the window function top hat in Fourier space.
 * In real space W(x) = 1 if x is inside the survey volume and W(x) = 0 otherwise.
 *
 */
//------------------------------------------------------------------------------

class WindowFunctionTopHat: public WindowFunctionBase
{
   private:
      /// survey linear size and volume
      double _L, _V;

   public:
      /// constructor
      WindowFunctionTopHat(double L): _L(L), _V(L * L * L) {}

      /// returns the window function value
      double operator()(ThreeVector x) { return _V * _Wi(x.p1()) * _Wi(x.p2()) * _Wi(x.p3()); }

   private:
      /// helper function
      /// Fourier transform of a single component
      double _Wi(double xi) {
         if (xi == 0) return 1;
         else return sin(xi * _L/2) / (xi * _L/2);
      }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // WINDOW_FUNCTION_TOPHAT_HPP

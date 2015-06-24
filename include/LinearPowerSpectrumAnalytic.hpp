//------------------------------------------------------------------------------
/// \file LinearPowerSpectrumAnalytic.hpp
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
//    Interface of class LinearPowerSpectrumAnalytic
//------------------------------------------------------------------------------

#ifndef LINEAR_POWER_SPECTRUM_ANALYTIC_HPP
#define LINEAR_POWER_SPECTRUM_ANALYTIC_HPP

#include <cmath>

#include "LinearPowerSpectrumBase.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class LinearPowerSpectrumAnalytic
 *
 * \brief class for linear power spectra using an analytic power law formula
 *
 * Provides functions:
 * - to evaluate the power spectrum
 */
//------------------------------------------------------------------------------

class LinearPowerSpectrumAnalytic : public LinearPowerSpectrumBase
{
   private:
      int _n;           ///< exponent in power law

   public:
      /// constructor
      LinearPowerSpectrumAnalytic(int n) : _n(n) {}
      /// destructor
      virtual ~LinearPowerSpectrumAnalytic() {}

      /// returns the linear power spectrum
      double operator()(double x) { return pow(x, _n); }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // LINEAR_POWER_SPECTRUM_ANALYTIC_HPP

//------------------------------------------------------------------------------
/// \file LinearPowerSpectrumBase.hpp
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
//    Interface of class LinearPowerSpectrumBase
//------------------------------------------------------------------------------

#ifndef LINEAR_POWER_SPECTRUM_BASE_HPP
#define LINEAR_POWER_SPECTRUM_BASE_HPP

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class LinearPowerSpectrumBase
 *
 * \brief Base class for linear power spectra
 *
 * Provides virtual functions:
 * - to evaluate the power spectrum
 */
//------------------------------------------------------------------------------

class LinearPowerSpectrumBase
{
   public:
      /// returns the linear power spectrum
      virtual double operator()(double x) = 0;
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // LINEAR_POWER_SPECTRUM_BASE_HPP

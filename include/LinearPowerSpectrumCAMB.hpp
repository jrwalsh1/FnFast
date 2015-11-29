//------------------------------------------------------------------------------
/// \file LinearPowerSpectrumCAMB.hpp
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
//    Interface of class LinearPowerSpectrumCAMB
//------------------------------------------------------------------------------

#ifndef LINEAR_POWER_SPECTRUM_CAMB_HPP
#define LINEAR_POWER_SPECTRUM_CAMB_HPP

#include <cmath>
#include <string>
#include <vector>

#include "LinearPowerSpectrumBase.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fit.h>

namespace fnfast {

//------------------------------------------------------------------------------
/**
 * \class LinearPowerSpectrumCAMB
 *
 * \brief class for linear power spectra using input from CAMB
 *
 * Provides functions:
 * - to evaluate the power spectrum
 */
//------------------------------------------------------------------------------

class LinearPowerSpectrumCAMB : public LinearPowerSpectrumBase
{
   private:
      std::string _input_file;                       ///< input file
      std::vector<double> _kvec, _kvec_patches;      ///< data storage vectors
      std::vector<double> _Pvec, _Pvec_patches;      ///< data storage vectors
      gsl_interp_accel* _accel_ptr;                  ///< interpolation objects in gsl
      gsl_spline* _spline_ptr;                       ///< interpolation objects in gsl
      double _c0_low, _c1_low, _c0_high, _c1_high;   ///< fit parameters to define high and low k patches

   public:
      /// constructors
      LinearPowerSpectrumCAMB(const std::string& input_file);
      /// destructor
      virtual ~LinearPowerSpectrumCAMB()
      {
          gsl_spline_free(_spline_ptr);
          gsl_interp_accel_free(_accel_ptr);
      }

      /// returns the linear power spectrum
      double operator()(double x);

   private:
      /// Helper function to generate points equally spaced in log
      std::vector<double> _log_gen(double xmin, double xmax, int n);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // LINEAR_POWER_SPECTRUM_CAMB_HPP

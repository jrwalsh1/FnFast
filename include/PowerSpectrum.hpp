//------------------------------------------------------------------------------
/// \file PowerSpectrum.hpp
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
//    Interface of class PowerSpectrum
//------------------------------------------------------------------------------

#ifndef POWER_SPECTRUM_HPP
#define POWER_SPECTRUM_HPP

#include <vector>
#include <map>

#include "Diagram.hpp"
#include "SPTkernels.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class PowerSpectrum
 *
 * \brief class to calculate the power spectrum
 *
 * PowerSpectrum(Order)
 *
 * Contains the power spectrum at various levels:
 * - tree level:
 *    - differential in k
 *    - integrated over k
 * - one loop (only available if called with Order oneLoop)
 *    - differential in k, q
 *    - integrated over q, differential in k
 *    - integrated over q, k
 *
 * Provides functions for access to the power spectrum at these levels
 */
//------------------------------------------------------------------------------

class PowerSpectrum
{
   private:
      Order _order;                    ///< order of the calculation
      LinearPowerSpectrumBase* _PL;    ///< the linear power spectrum used in the calculation
      SPTkernels* _SPTkernels;         ///< SPT kernels instance
      vector<Diagram> _tree;           ///< tree level diagrams
      vector<Diagram> _loop;           ///< loop diagrams
      vector<Diagram> _cterms;         ///< counterterms

   public:
      /// constructor
      PowerSpectrum(Order order, LinearPowerSpectrumBase* PL);
      /// destructor
      virtual ~PowerSpectrum() {}
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // POWER_SPECTRUM_HPP

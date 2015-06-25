//------------------------------------------------------------------------------
/// \file Trispectrum.hpp
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
//    Interface of class Trispectrum
//------------------------------------------------------------------------------

#ifndef TRISPECTRUM_HPP
#define TRISPECTRUM_HPP

#include <vector>
#include <map>

#include "Diagram.hpp"
#include "SPTkernels.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class Trispectrum
 *
 * \brief class to calculate the trispectrum
 *
 * Trispectrum(Order)
 *
 * Contains the trispectrum at various levels:
 * - tree level:
 *    - differential in k_i
 *    - integrated over k_i
 * - one loop (only available if called with Order oneLoop)
 *    - differential in k_i, q
 *    - integrated over q, differential in k_i
 *    - integrated over q, k_i
 *
 * Provides functions for access to the trispectrum at these levels
 */
//------------------------------------------------------------------------------

class Trispectrum
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
      Trispectrum(Order order, LinearPowerSpectrumBase* PL);
      /// destructor
      virtual ~Trispectrum() {}
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // TRISPECTRUM_HPP

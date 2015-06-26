//------------------------------------------------------------------------------
/// \file Bispectrum.hpp
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
//    Interface of class Bispectrum
//------------------------------------------------------------------------------

#ifndef BISPECTRUM_HPP
#define BISPECTRUM_HPP

#include <vector>
#include <map>

#include "Diagram.hpp"
#include "SPTkernels.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class Bispectrum
 *
 * \brief class to calculate the bispectrum
 *
 * Bispectrum(Order)
 *
 * Contains the bispectrum at various levels:
 * - tree level:
 *    - differential in k_i
 *    - integrated over k_i
 * - one loop (only available if called with Order oneLoop)
 *    - differential in k_i, q
 *    - integrated over q, differential in k_i
 *    - integrated over q, k_i
 *
 * Provides functions for access to the bispectrum at these levels
 */
//------------------------------------------------------------------------------

class Bispectrum
{
   public:
      enum Graphs {
         B211,
         B411,
         B321a,
         B321b,
         B222
      };

   private:
      Order _order;                          ///< order of the calculation
      LinearPowerSpectrumBase* _PL;          ///< the linear power spectrum used in the calculation
      SPTkernels* _SPTkernels;               ///< SPT kernels instance
      vector<Diagram*> _tree;                ///< tree level diagrams
      vector<Diagram*> _loop;                ///< loop diagrams
      vector<Diagram*> _cterms;              ///< counterterms
      map<Graphs, Diagram*> _diagrams;       ///< container for diagrams

   public:
      /// constructor
      Bispectrum(Order order, LinearPowerSpectrumBase* PL);
      /// destructor
      virtual ~Bispectrum() {}

      /// access diagrams
      Diagram* operator[](Graphs graph) { return _diagrams[graph]; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // BISPECTRUM_HPP

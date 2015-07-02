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

#include "Momentum.hpp"
#include "Diagram.hpp"
#include "SPTkernels.hpp"
#include "EFTkernels.hpp"
#include "WindowFunctionBase.hpp"

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
         B222,
         B411x,
         B321bx
      };

   private:
      Order _order;                                   ///< order of the calculation
      LinearPowerSpectrumBase* _PL;                   ///< the linear power spectrum used in the calculation
      SPTkernels* _SPTkernels;                        ///< SPT kernels instance
      EFTkernels* _EFTkernels;                        ///< EFT kernels instance
      vector<Diagram*> _tree;                         ///< tree level diagrams
      vector<Diagram*> _loop;                         ///< loop diagrams
      vector<Diagram*> _cterms;                       ///< counterterms
      map<Graphs, Diagram*> _diagrams;                ///< container for diagrams
      EFTcoefficients* _eftcoefficients;              ///< EFT coefficients
      vector<Momenta::MomentumLabel> _labels;         ///< external momenta labels
      DiagramMomenta _momenta;                        ///< diagram momenta
      double _UVcutoff;                               ///< UV cutoff for loop integrations
      double _kBin;                                   ///< size of k bins
      WindowFunctionBase* _W;                         ///< Window function

   public:
      /// constructor
      Bispectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients);
      /// destructor
      virtual ~Bispectrum() {}

      /// access diagrams
      Diagram* operator[](Graphs graph) { return _diagrams[graph]; }
    
      /// set size of k bins
      void set_kBinSize(double kBin) { _kBin = kBin;}
      /// set window function
      void set_windowFunction(WindowFunctionBase* W) { _W = W; }

      /// access results
      /// Differential in k
      /// tree level
      double treeLevel_value(ThreeVector k2, ThreeVector k3);
      /// one loop differential in q and integrated in q
      double oneLoopSPT_value(ThreeVector k2, ThreeVector k3, ThreeVector q);
      double oneLoopSPT_value(ThreeVector k2, ThreeVector k3);
      /// one loop counterterms
      double oneLoopCterms_value(ThreeVector k2, ThreeVector k3);

      /// Averaged over k bins
      /// tree level
      double treeLevel_value(double k2, double k3);
      /// one loop integrated in q
      double oneLoopSPT_value(double k2, double k3);
      /// one loop counterterms
      double oneLoopCterms_value(double k2, double k3);

      /// Averaged over k bins + convolution with window function
      /// tree level
      double treeLevel_value_win(double k2, double k3);
      /// one loop integrated in q
      double oneLoopSPT_value_win(double k2, double k3);
      /// one loop counterterms
      double oneLoopCterms_value_win(double k2, double k3);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // BISPECTRUM_HPP

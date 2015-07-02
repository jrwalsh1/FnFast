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

#include <map>

#include "Momentum.hpp"
#include "Diagram.hpp"
#include "SPTkernels.hpp"
#include "EFTkernels.hpp"
#include "WindowFunctionBase.hpp"

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
   public:
      enum Graphs {
         P11,
         P31,
         P22,
         P31x
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
      PowerSpectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients);
      /// destructor
      virtual ~PowerSpectrum() {}

      /// access diagrams
      Diagram* operator[](Graphs graph) { return _diagrams[graph]; }

      /// set size of k bins
      void set_kBinSize(double kBin) { _kBin = kBin; }
      /// set window function
      void set_windowFunction(WindowFunctionBase* W) { _W = W; }

      /// access results
      /// Differential in k
      /// tree level
      double treeLevel_value(ThreeVector k);
      /// one loop differential in q and integrated in q
      double oneLoopSPT_value(ThreeVector k, ThreeVector q);
      double oneLoopSPT_value(ThreeVector k);
      /// one loop counterterms
      double oneLoopCterms_value(ThreeVector k);
    
      /// Averaged over k bins
      /// tree level
      double treeLevel_value(double k);
      /// one loop integrated in q
      double oneLoopSPT_value(double k);
      /// one loop counterterms
      double oneLoopCterms_value(double k);
    
      /// Averaged over k bins + convolution with window function
      /// tree level
      double treeLevel_value_win(double k);
      /// one loop integrated in q
      double oneLoopSPT_value_win(double k);
      /// one loop counterterms
      double oneLoopCterms_value_win(double k);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // POWER_SPECTRUM_HPP

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
         B321ax
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
      double _UVcutoff;                               ///< UV cutoff for loop integrations
      double _kBin;                                   ///< size of k bins
      WindowFunctionBase* _W;                         ///< Window function
      int _seed;                                      ///< random number seed for VEGAS

      /// phase space for the loop momentum
      class LoopPhaseSpace {
         private:
            double _qmax;           ///< upper limit on q integral
            double _jacobian;       ///< jacobian for the phase space point
            ThreeVector _q;         ///< loop momentum value

         public:
            static constexpr double pi = 3.14159265358979;      ///< pi

         public:
            /// constructor
            LoopPhaseSpace(double qmax) : _qmax(qmax) {}
            /// destructor
            virtual ~LoopPhaseSpace() {}

            /// set the loop phase space, returns the jacobian
            double setPS(double qpts[3]);

            /// returns the loop momentum
            ThreeVector q() { return _q; }
      };

      /// container for the integration options
      struct LoopIntegrationOptions {
         double k1;
         double k2;
         double costheta12;
         Bispectrum* bispectrum;
         LoopPhaseSpace* loopPS;
      };

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

      /// set the loop momentum cutoff
      void set_qmax(double qmax);

      /// set the random number seed
      void set_seed(int seed) { _seed = seed; }

      /// get results
      /// Differential in k
      /// tree level
      double tree(double k1, double k2, double costheta12);
      /// one loop differential in q and integrated in q
      double loopSPT_excl(ThreeVector k1, ThreeVector k2, ThreeVector q);
      double loopSPT(double k1, double k2, double costheta12);
      /// one loop counterterms
      double ctermsEFT(double k1, double k2, double costheta12);

      /// Averaged over k bins
      /// tree level
      double tree_kbin(double k1, double k2);
      /// one loop integrated in q
      double loopSPT_kbin(double k1, double k2);
      /// one loop counterterms
      double ctermsEFT_kbin(double k1, double k2);

      /// Averaged over k bins + convolution with window function
      /// tree level
      double tree_kbin_win(double k1, double k2);
      /// one loop integrated in q
      double loopSPT_kbin_win(double k1, double k2);
      /// one loop counterterms
      double ctermsEFT_kbin_win(double k1, double k2);

   private:
      /// loop integrand function
      static int loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline void Bispectrum::set_qmax(double qmax) {
   for (size_t c = 0; c < _loop.size(); c++) {
      _loop[c]->set_qmax(qmax);
   }
}

#endif // BISPECTRUM_HPP

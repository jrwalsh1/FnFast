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

#include "Momentum.hpp"
#include "Diagram.hpp"
#include "SPTkernels.hpp"
#include "EFTkernels.hpp"
#include "WindowFunctionBase.hpp"

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
   public:
      enum Graphs {
         T3111,
         T2211,
         T5111,
         T4211a,
         T4211b,
         T3311a,
         T3311b,
         T3221a,
         T3221b,
         T3221c,
         T2222,
         T5111x,
         T4211ax,
         T3311ax,
         T3221ax
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
            double set_loopPS(double qpts[3]);

            /// returns the loop momentum
            ThreeVector q() { return _q; }
      };

      /// container for the integration options
      struct LoopIntegrationOptions {
         double k;
         double kp;
         double costheta;
         Trispectrum* trispectrum;
         LoopPhaseSpace* loopPS;
      };

      /// container for the integration options
      struct AngularIntegrationOptions {
         double k;
         double kp;
         Trispectrum* trispectrum;
      };

   public:
      /// constructor
      Trispectrum(Order order, LinearPowerSpectrumBase* PL, EFTcoefficients* eftcoefficients);
      /// destructor
      virtual ~Trispectrum() {}

      /// access diagrams
      Diagram* operator[](Graphs graph) { return _diagrams[graph]; }
   
      /// set size of k bins
      void set_kBinSize(double kBin) { _kBin = kBin;}
      /// set window function
      void set_windowFunction(WindowFunctionBase* W) { _W = W; }

      /// get results
      /// covariance limit, differential in k
      /// tree level
      double cov_tree(double k, double kp, double costheta);
      /// one loop integrated in q
      double cov_loopSPT(double k, double kp, double costheta);
      /// one loop counterterms
      double cov_ctermsEFT(double k, double kp, double costheta);

      /// covariance limit integrated over angles
      /// tree level
      double cov_tree(double k, double kp);
      /// full one loop integrated in q
      double cov_loop(double k, double kp);

      /// Averaged over k bins
      /// tree level
      double cov_tree_kbin(double k2, double k3, double k4);
      /// one loop integrated in q
      double cov_loopSPT_kbin(double k2, double k3, double k4);
      /// one loop counterterms
      double cov_ctermsEFT_kbin(double k2, double k3, double k4);

      /// Averaged over k bins + convolution with window function
      /// tree level
      double cov_tree_kbin_win(double k2, double k3, double k4);
      /// one loop integrated in q
      double cov_loopSPT_kbin_win(double k2, double k3, double k4);
      /// one loop counterterms
      double cov_ctermsEFT_kbin_win(double k2, double k3, double k4);

      /// full trispectrum, differential in k
      /// tree level
      double tree(ThreeVector k1, ThreeVector k2, ThreeVector k3);
      /// one loop differential in q and integrated in q
      double loopSPT_excl(ThreeVector k1, ThreeVector k2, ThreeVector k3, ThreeVector q);
      double loopSPT(ThreeVector k1, ThreeVector k2, ThreeVector k3);
      /// one loop counterterms
      double ctermsEFT(ThreeVector k1, ThreeVector k2, ThreeVector k3);

   private:
      /// tree level: polar angle integrand function, covariance limit
      static int tree_angular_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
      /// loop integrand function, covariance limit
      static int loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
      /// polar angle integrand function, covariance limit
      static int angular_loop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // TRISPECTRUM_HPP

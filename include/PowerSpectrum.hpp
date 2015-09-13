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

#include "MomentumMap.hpp"
#include "VertexMap.hpp"
#include "DiagramSet2point.hpp"
#include "KernelBase.hpp"
#include "Integration.hpp"

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
 * - one loop
 *    - differential in k, q
 *    - integrated over q, differential in k
 * - two loop
 *    - differential in k, q
 *    - integrated over q, differential in k
 *
 * Provides functions for access to the power spectrum at these levels
 */
//------------------------------------------------------------------------------

class PowerSpectrum
{
   private:
      Order _order;                       ///< order of the calculation
      DiagramSet2point _diagrams;         ///< 2-point diagrams
      double _UVcutoff;                   ///< UV cutoff for loop integrations
      int _seed;                          ///< random number seed for VEGAS

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
            double setPS(double qpts[2]);

            /// returns the loop momentum
            ThreeVector q() { return _q; }
      };

      /// container for the integration options
      struct LoopIntegrationOptions {
         double k;
         const PowerSpectrum* powerspectrum;
         LoopPhaseSpace* loopPS;
         const VertexMap<KernelBase*>* kernels;
         LinearPowerSpectrumBase* PL;
      };

   public:
      /// constructor
      PowerSpectrum(Order order);
      /// destructor
      virtual ~PowerSpectrum() {}

      /// access diagrams
      DiagramBase* operator[](DiagramSet2point::Graphs_2point graph) { return _diagrams[graph]; }

      /// set the loop momentum cutoff
      void set_qmax(double qmax) { _diagrams.set_qmax(qmax); }

      /// set the random number seed
      void set_seed(int seed) { _seed = seed; }

      /// get results differential in k
      /// tree level
      double tree(double k, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
      /// one loop differential in q and integrated over q
      double oneLoop_excl(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const;
      IntegralResult oneLoop(double k, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
      /// two loop differential in q, q2 and integrated over q, q2
      double twoLoop_excl(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const;
      IntegralResult twoLoop(double k, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

   private:
      /// one loop integrand
      static int oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
      /// two loop integrand
      static int twoLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // POWER_SPECTRUM_HPP

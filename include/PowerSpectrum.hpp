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

#include "DiagramSet2pointSPT.hpp"
#include "KernelBase.hpp"
#include "Integration.hpp"

namespace fnfast {

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
      DiagramSet2pointSPT _diagrams;      ///< 2-point diagrams
      double _UVcutoff;                   ///< UV cutoff for loop integrations
      int _seed;                          ///< random number seed for VEGAS

      /// container for the integration options
      struct LoopPhaseSpace
      {
         int ndim;
         double k;
         LabelMap<Momentum, ThreeVector> momenta;
         double qmax;
         const LabelMap<Vertex, KernelBase*>* kernels;
         LinearPowerSpectrumBase* PL;
         const PowerSpectrum* powerspectrum;
         static constexpr double pi = 3.14159265358979;

         /// constructor
         LoopPhaseSpace(double k, double qlim, const LabelMap<Vertex, KernelBase*>* kern, LinearPowerSpectrumBase* linPS, const PowerSpectrum* powerspec);

         /// sample phase space; return the point along with the jacobian
         std::pair<double, LabelMap<Momentum, ThreeVector>* const> generate_point_oneLoop(std::vector<double> xpts);
         std::pair<double, LabelMap<Momentum, ThreeVector>* const> generate_point_twoLoop(std::vector<double> xpts);
      };

   public:
      /// constructor
      PowerSpectrum(Order order);
      /// destructor
      virtual ~PowerSpectrum() {}

      /// access diagrams
      const DiagramSetBase* diagrams() const { return &_diagrams; }
      DiagramBase* operator[](Graphs_2point graph) { return _diagrams[graph]; }

      /// set the loop momentum cutoff
      void set_qmax(double qmax) { _diagrams.set_qmax(qmax); }

      /// set the random number seed
      void set_seed(int seed) { _seed = seed; }

      /// get results differential in k
      /// tree level
      double tree(double k, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
      /// one loop integrated over q
      IntegralResult oneLoop(double k, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
      /// two loop integrated over q, q2
      IntegralResult twoLoop(double k, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

   private:
      /// one loop integrand
      static int oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
      /// two loop integrand
      static int twoLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline PowerSpectrum::LoopPhaseSpace::LoopPhaseSpace(double kmag, double qlim, const LabelMap<Vertex, KernelBase*>* kern, LinearPowerSpectrumBase* linPS, const PowerSpectrum* powerspec)
: ndim(2), k(kmag), momenta(LabelMap<Momentum, ThreeVector> {{Momentum::k1, ThreeVector(0, 0, -k)}, {Momentum::k2, ThreeVector(0, 0, -k)}, {Momentum::q, ThreeVector()}}), qmax(qlim), kernels(kern), PL(linPS), powerspectrum(powerspec)
{}

} // namespace fnfast

#endif // POWER_SPECTRUM_HPP

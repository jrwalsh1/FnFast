//------------------------------------------------------------------------------
/// \file Bispectrum.hpp
//
// Author(s):
//    Jon Walsh
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
//    Interface of class Bispectrum
//------------------------------------------------------------------------------

#ifndef BISPECTRUM_HPP
#define BISPECTRUM_HPP

#include "DiagramSet3pointSPT.hpp"
#include "KernelBase.hpp"
#include "Integration.hpp"

namespace fnfast {

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
 *    - differential in k1, k2 (magnitudes) and theta12
 * - one loop
 *    - differential in k1, k2 (magnitudes) and theta12, q
 *    - integrated over q, differential in k1, k2 (magnitudes) and theta12
 *
 * Provides functions for access to the bispectrum at these levels
 */
//------------------------------------------------------------------------------

class Bispectrum
{
   private:
      Order _order;                       ///< order of the calculation
      DiagramSet3pointSPT _diagrams;      ///< 3-point diagrams
      double _UVcutoff;                   ///< UV cutoff for loop integrations
      int _seed;                          ///< random number seed for VEGAS

      /// container for the integration options
      struct LoopPhaseSpace
      {
         int ndim;
         double k1;
         double k2;
         double theta12;
         LabelMap<Momentum, ThreeVector> momenta;
         double qmax;
         const LabelMap<Vertex, KernelBase*>* kernels;
         LinearPowerSpectrumBase* PL;
         const Bispectrum* bispectrum;
         static constexpr double pi = 3.14159265358979;

         /// constructor
         LoopPhaseSpace(double k1, double k2, double theta12, double qlim, const LabelMap<Vertex, KernelBase*>* kern, LinearPowerSpectrumBase* linPS, const Bispectrum* bispec);

         /// sample phase space; return the point along with the jacobian
         std::pair<double, LabelMap<Momentum, ThreeVector>* const> generate_point_oneLoop(std::vector<double> xpts);
      };

   public:
      /// constructor
      Bispectrum(Order order);
      /// destructor
      virtual ~Bispectrum() {}

      /// access diagrams
      const DiagramSetBase* diagrams() const { return &_diagrams; }
      DiagramBase* operator[](Graphs_3point graph) { return _diagrams[graph]; }

      /// set the loop momentum cutoff
      void set_qmax(double qmax) { _diagrams.set_qmax(qmax); }

      /// set the random number seed
      void set_seed(int seed) { _seed = seed; }

      /// get results differential in k
      /// tree level
      double tree(double k1, double k2, double theta12, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
      /// one loop integrated over q
      IntegralResult oneLoop(double k1, double k2, double theta12, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

   private:
      /// one loop integrand
      static int oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline Bispectrum::LoopPhaseSpace::LoopPhaseSpace(double k1mag, double k2mag, double theta12val, double qlim, const LabelMap<Vertex, KernelBase*>* kern, LinearPowerSpectrumBase* linPS, const Bispectrum* bispec)
: ndim(3), k1(k1mag), k2(k2mag), theta12(theta12val), momenta(LabelMap<Momentum, ThreeVector> {{Momentum::k1, ThreeVector()}, {Momentum::k2, ThreeVector()}, {Momentum::k3, ThreeVector()}, {Momentum::q, ThreeVector()}}), qmax(qlim), kernels(kern), PL(linPS), bispectrum(bispec)
{
   // set the external momenta
   momenta[Momentum::k1] = ThreeVector(0, 0, k1);
   momenta[Momentum::k2] = ThreeVector(k2 * sin(theta12val), 0, k2 * cos(theta12val));
   momenta[Momentum::k3] = -momenta[Momentum::k1] - momenta[Momentum::k2];
}

} // namespace fnfast

#endif // BISPECTRUM_HPP

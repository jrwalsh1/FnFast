//------------------------------------------------------------------------------
/// \file Covariance.hpp
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
//    Interface of class Covariance
//------------------------------------------------------------------------------

#ifndef COVARIANCE_HPP
#define COVARIANCE_HPP

#include "DiagramSet4pointSPT.hpp"
#include "KernelBase.hpp"
#include "Integration.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
/**
 * \class Covariance
 *
 * \brief class to calculate the covariance
 *
 * Covariance(Order)
 *
 * Contains the covariance at various levels:
 * - tree level:
 *    - differential in k, k' (magnitudes) and theta
 * - one loop
 *    - differential in k, k' (magnitudes) and theta, q
 *    - integrated over q, differential in k, k' (magnitudes) and theta
 *
 * Provides functions for access to the bispectrum at these levels
 */
//------------------------------------------------------------------------------

class Covariance
{
   private:
      Order _order;                       ///< order of the calculation
      DiagramSet4pointSPT _diagrams;      ///< 4-point diagrams
      double _UVcutoff;                   ///< UV cutoff for loop integrations
      int _seed;                          ///< random number seed for VEGAS

      /// container for the integration options
      struct PhaseSpace
      {
         int ndim;
         double k;
         double kprime;
         LabelMap<Momentum, ThreeVector> momenta;
         double qmax;
         const LabelMap<Vertex, KernelBase*>* kernels;
         LinearPowerSpectrumBase* PL;
         const Covariance* covariance;
         static constexpr double pi = 3.14159265358979;

         /// constructor
         PhaseSpace(double kmag, double kprimemag, double qlim, const LabelMap<Vertex, KernelBase*>* kern, LinearPowerSpectrumBase* linPS, const Covariance* cov);

         /// sample phase space; return the point along with the jacobian
         std::pair<double, LabelMap<Momentum, ThreeVector>* const> generate_point_tree(std::vector<double> xpts);
         std::pair<double, LabelMap<Momentum, ThreeVector>* const> generate_point_oneLoop(std::vector<double> xpts);
      };

   public:
      /// constructor
      Covariance(Order order);
      /// destructor
      virtual ~Covariance() {}

      /// access diagrams
      const DiagramSetBase* diagrams() const { return &_diagrams; }
      DiagramBase* operator[](Graphs_4point graph) { return _diagrams[graph]; }

      /// set the loop momentum cutoff
      void set_qmax(double qmax) { _diagrams.set_qmax(qmax); }

      /// set the random number seed
      void set_seed(int seed) { _seed = seed; }

      /// get results differential in k
      /// tree level
      IntegralResult tree(double k, double kprime, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
      /// one loop integrated over q, theta
      IntegralResult oneLoop(double k, double kprime, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

   private:
      /// one loop integrand
      static int tree_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
      /// one loop integrand
      static int oneLoop_integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline Covariance::PhaseSpace::PhaseSpace(double kmag, double kprimemag, double qlim, const LabelMap<Vertex, KernelBase*>* kern, LinearPowerSpectrumBase* linPS, const Covariance* cov)
: k(kmag), kprime(kprimemag), momenta(LabelMap<Momentum, ThreeVector> {{Momentum::k1, ThreeVector()}, {Momentum::k2, ThreeVector()}, {Momentum::k3, ThreeVector()}, {Momentum::k4, ThreeVector()}, {Momentum::q, ThreeVector()}}), qmax(qlim), kernels(kern), PL(linPS), covariance(cov)
{}

} // namespace fnfast

#endif // COVARIANCE_HPP

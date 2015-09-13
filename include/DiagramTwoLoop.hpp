//------------------------------------------------------------------------------
/// \file DiagramTwoLoop.hpp
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
//    Interface of class DiagramTwoLoop
//------------------------------------------------------------------------------

#ifndef DIAGRAM_TWO_LOOP_HPP
#define DIAGRAM_TWO_LOOP_HPP

#include "DiagramBase.hpp"

using namespace std;

class DiagramTwoLoop : public DiagramBase
{
   private:
      Order _order;                   ///< order of the calculation
      vector<Propagator> _IRpoles;    ///< IR poles
      double _qmax;                   ///< upper limit on the magnitude of the loop momentum (default is infinity)

   public:
      /// constructor
      DiagramTwoLoop(vector<Line> lines);
      /// destructor
      virtual ~DiagramTwoLoop() {}

      /// returns the diagram value with the input momentum routing
      double value_base(const MomentumMap<ThreeVector>& mom, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

      /// returns the IR regulated diagram value with the input momentum routing
      double value_base_IRreg(const MomentumMap<ThreeVector>& mom, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

      /// returns the IR regulated diagram value, symmetrized over external momenta
      double value(const MomentumMap<ThreeVector>& mom, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

      /// set a cutoff on the magnitude of the loop momentum
      void set_qmax(double qmax) { _qmax = qmax; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // DIAGRAM_TWO_LOOP_HPP

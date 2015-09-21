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

namespace fnfast {

class DiagramTwoLoop : public DiagramBase
{
   private:
      Order _order;                        ///< order of the calculation
      std::vector<Propagator> _IRpoles;    ///< IR poles
      double _qmax;                        ///< upper limit on the magnitude of the loop momentum (default is infinity)

   public:
      /// base constructor, assumes all vertices are the same type and we're computing delta correlators
      DiagramTwoLoop(std::vector<Line> lines);
      /// constructor with vertex types specified, assumes we're computing delta correlators
      DiagramTwoLoop(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes);
      /// constructor specifying vertex and kernel types
      DiagramTwoLoop(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes, LabelMap<Vertex, KernelType> kerneltypes);
      /// destructor
      virtual ~DiagramTwoLoop() {}

      /// returns the diagram value with the input momentum routing
      double value_base(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

      /// returns the IR regulated diagram value with the input momentum routing
      double value_base_IRreg(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

      /// returns the IR regulated diagram value, symmetrized over external momenta
      double value(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;

      /// set a cutoff on the magnitude of the loop momentum
      void set_qmax(double qmax) { _qmax = qmax; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // DIAGRAM_TWO_LOOP_HPP

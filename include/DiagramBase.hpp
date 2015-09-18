//------------------------------------------------------------------------------
/// \file DiagramBase.hpp
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
//    Definition of base class DiagramBase
//------------------------------------------------------------------------------

#ifndef DIAGRAM_BASE_HPP
#define DIAGRAM_BASE_HPP

#include <string>
#include <map>
#include <vector>
#include <unordered_map>

#include "KernelBase.hpp"
#include "LinearPowerSpectrumBase.hpp"
#include "Line.hpp"
#include "LabelMap.hpp"

namespace fnfast {

class DiagramBase
{
   protected:
      Order _order;                                               ///< order of the calculation
      std::vector<Vertex> _vertices;                              ///< list of vertices in the diagram
      double _symfac;                                             ///< symmetry factor
      std::vector<Line> _lines;                                   ///< lines in the graph
      LabelMap<Vertex, std::vector<Propagator> > _vertexmomenta;  ///< vertex factors
      std::vector<VertexPair> _vertexpairs;                       ///< container for endpoint vertices of lines
      std::vector<Momentum> _extmomlabels;                        ///< momentum labels in the graph
      std::vector<LabelMap<Momentum, Momentum> > _perms;          ///< permutations of external momenta for the graph

   public:
      /// constructor
      DiagramBase(std::vector<Line> lines);
      /// destructor
      virtual ~DiagramBase() {}

      /// return the symmetry factor
      double symmetry_factor() const { return _symfac; }

      /// return the number of permutations
      int nperms() const { return _perms.size(); }

      /// get the external momentum permutations used in the diagram calculation
      std::vector<LabelMap<Momentum, Momentum> > get_perms() const { return _perms; }

      /// set the external momentum permutations to be used in the diagram calculation
      void set_perms(std::vector<LabelMap<Momentum, Momentum> > perms) { _perms = perms; }

      /// returns the diagram value with the input momentum routing
      virtual double value(const LabelMap<Momentum, ThreeVector>& mom, const LabelMap<Vertex, KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const = 0;

   protected:
      /// function theta(|p1| < |p2|)
      static double theta(ThreeVector p1, ThreeVector p2);

      /// computes symmetry factor
      double calc_symmetry_factor();

      /// computes permutations of external momenta
      std::vector<LabelMap<Momentum, Momentum> > calc_permutations();

      /// factorial
      static int factorial(int n) { return (n == 0 || n == 1) ? 1 : n * factorial(n-1); }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline double DiagramBase::theta(ThreeVector p1, ThreeVector p2)
{
   if (p1.square() < p2.square()) { return 1; }
   else { return 0; }
}

} // namespace fnfast

#endif // DIAGRAM_BASE_HPP

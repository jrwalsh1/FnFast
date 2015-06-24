//------------------------------------------------------------------------------
/// \file Diagram.hpp
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
//    Interface of class Diagram
//------------------------------------------------------------------------------

#ifndef DIAGRAM_HPP
#define DIAGRAM_HPP

#include <string>
#include <map>
#include <vector>

#include "KernelBase.hpp"
#include "LinearPowerSpectrumBase.hpp"
#include "Line.hpp"

using namespace std;

class Diagram
{
   private:
      Order _order;                                                        ///< order of the calculation
      double _symfac;                                                      ///< symmetry factor
      vector<Line> _lines;                                                 ///< lines in the graph
      map<Vertices::VertexLabel, vector<Propagator> > _vertexmomenta;      ///< vertex factors
      vector<Propagator> _IRpoles;                                         ///< IR poles
      map<Vertices::VertexLabel, KernelBase*> _kernels;                    ///< pointer to a kernel instance
      LinearPowerSpectrumBase* _PL;                                        ///< pointer to a linear power spectrum instance

   public:
      /// constructor
      Diagram(vector<Line> lines, map<Vertices::VertexLabel, KernelBase*> kernels, LinearPowerSpectrumBase* PL);
      /// destructor
      virtual ~Diagram() {}

      /// set the linear power spectrum
      void setLinearPowerSpectrum(LinearPowerSpectrumBase* PL);

      /// set the kernels
      void setKernels(map<Vertices::VertexLabel, KernelBase*> kernels);

      /// returns the diagram value with the input momentum routing
      double value_base(DiagramMomenta mom);

      /// returns the IR regulated diagram value with the input momentum routing
      double value_base_IRreg(DiagramMomenta mom);

      /// returns the IR regulated diagram value, symmetrized over external momenta
      double value_IRreg(DiagramMomenta mom);

   private:
      /// function theta(|p1| < |p2|)
      static double theta(ThreeVector p1, ThreeVector p2);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline void Diagram::setLinearPowerSpectrum(LinearPowerSpectrumBase* PL)
{
   _PL = PL;
}

//------------------------------------------------------------------------------
inline void Diagram::setKernels(map<Vertices::VertexLabel, KernelBase*> kernels)
{
   _kernels = kernels;
}

//------------------------------------------------------------------------------
inline double Diagram::theta(ThreeVector p1, ThreeVector p2)
{
   if (p1.square() < p2.square()) { return 1; }
   else { return 0; }
}

#endif // DIAGRAM_HPP

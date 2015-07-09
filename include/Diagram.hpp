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
#include <unordered_map>
#include <limits>

#include "KernelBase.hpp"
#include "LinearPowerSpectrumBase.hpp"
#include "Line.hpp"

using namespace std;

class Diagram
{
   private:
      Order _order;                                                                    ///< order of the calculation
      vector<Vertices::VertexLabel> _vertices;                                         ///< list of vertices in the diagram
      double _symfac;                                                                  ///< symmetry factor
      vector<Line> _lines;                                                             ///< lines in the graph
      unordered_map<Vertices::VertexLabel, vector<Propagator> > _vertexmomenta;        ///< vertex factors
      vector<Propagator> _IRpoles;                                                     ///< IR poles
      unordered_map<Vertices::VertexLabel, KernelBase*> _kernels;                      ///< pointer to a kernel instance
      LinearPowerSpectrumBase* _PL;                                                    ///< pointer to a linear power spectrum instance
      vector<VertexPair> _vertexpairs;                                                 ///< container for endpoint vertices of lines
      vector<Momenta::MomentumLabel> _extmomlabels;                                    ///< momentum labels in the graph
      vector<unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> > _perms;   ///< permutations of external momenta for the graph
      double _qmax;                                                                    ///< upper limit on the magnitude of the loop momentum (default is infinity)

   public:
      /// constructor
      Diagram(vector<Line> lines, unordered_map<Vertices::VertexLabel, KernelBase*> kernels, LinearPowerSpectrumBase* PL);
      /// destructor
      virtual ~Diagram() {}

      /// return the symmetry factor
      double symmetry_factor() { return _symfac; }

      /// return the number of permutations
      int nperms() { return _perms.size(); }

      /// get the external momentum permutations used in the diagram calculation
      vector<unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> > get_perms() { return _perms; }

      /// set the external momentum permutations to be used in the diagram calculation
      void set_perms(vector<unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> > perms) { _perms = perms; }

      /// set the linear power spectrum
      void setLinearPowerSpectrum(LinearPowerSpectrumBase* PL);

      /// set the kernels
      void setKernels(unordered_map<Vertices::VertexLabel, KernelBase*> kernels);

      /// returns the diagram value with the input momentum routing
      double value_base(DiagramMomenta mom);

      /// returns the IR regulated diagram value with the input momentum routing
      double value_base_IRreg(DiagramMomenta mom);

      /// returns the IR regulated diagram value, symmetrized over external momenta
      double value_IRreg(DiagramMomenta mom);

      /// set a cutoff on the magnitude of the loop momentum
      void set_qmax(double qmax) { _qmax = qmax; }

   private:
      /// function theta(|p1| < |p2|)
      static double theta(ThreeVector p1, ThreeVector p2);

      /// computes symmetry factor
      double calc_symmetry_factor();

      /// computes permutations of external momenta
      vector<unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> > calc_permutations();

      /// factorial
      static int factorial(int n) { return (n == 0 || n == 1) ? 1 : n * factorial(n-1); }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

/*
//------------------------------------------------------------------------------
inline Diagram::~Diagram()
{
   delete _PL;
   for (auto kern : _kernels) {
      delete kern.second;
   }
}
*/

//------------------------------------------------------------------------------
inline void Diagram::setLinearPowerSpectrum(LinearPowerSpectrumBase* PL)
{
   _PL = PL;
}

//------------------------------------------------------------------------------
inline void Diagram::setKernels(unordered_map<Vertices::VertexLabel, KernelBase*> kernels)
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

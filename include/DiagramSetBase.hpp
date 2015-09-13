//------------------------------------------------------------------------------
/// \file DiagramSetBase.hpp
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
//    Definition of base class DiagramSetBase
//------------------------------------------------------------------------------

#ifndef DIAGRAM_SET_BASE_HPP
#define DIAGRAM_SET_BASE_HPP

#include "DiagramTree.hpp"
#include "DiagramOneLoop.hpp"
#include "DiagramTwoLoop.hpp"

using namespace std;

class DiagramSetBase
{
   protected:
      Order _order;                            ///< order of the calculation to work up to
      vector<DiagramTree*> _tree;              ///< tree diagrams
      vector<DiagramOneLoop*> _oneLoop;        ///< one loop diagrams
      vector<DiagramTwoLoop*> _twoLoop;        ///< two loop diagrams
      vector<MomentumLabel> _extmomlabels;     ///< external momentum labels in the graph

   public:
      /// constructor
      DiagramSetBase(Order order);
      /// destructor
      virtual ~DiagramSetBase() {}

      /// return the highest order in the set of diagrams
      Order order() const { return _order; }

      /// return the external momentum labels
      vector<MomentumLabel> external_labels() const { return _extmomlabels; }

      /// get the tree level diagrams
      vector<DiagramTree*> tree() const { return _tree; }

      /// get the one loop diagrams
      vector<DiagramOneLoop*> oneLoop() const { return _oneLoop; }

      /// get the two loop diagrams
      vector<DiagramTwoLoop*> twoLoop() const { return _twoLoop; }

      /// get the value of the tree level diagrams
      double value_tree(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const;

      /// get the value of the one loop diagrams
      double value_oneLoop(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const;

      /// get the value of the two loop diagrams
      double value_twoLoop(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const;

      /// set the loop momentum restriction for all loop diagrams
      void set_qmax(double qmax);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
inline DiagramSetBase::DiagramSetBase(Order order)
: _order(order) {}

//------------------------------------------------------------------------------
inline double DiagramSetBase::value_tree(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const
{
   double value = 0;
   for (auto diagram : _tree) {
      value += diagram->value(mom, kernels, PL);
   }
   return value;
}

//------------------------------------------------------------------------------
inline double DiagramSetBase::value_oneLoop(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const
{
   double value = 0;
   for (auto diagram : _oneLoop) {
      value += diagram->value(mom, kernels, PL);
   }
   return value;
}

//------------------------------------------------------------------------------
inline double DiagramSetBase::value_twoLoop(const MomentumMap<ThreeVector>& mom, VertexMap<KernelBase*> kernels, LinearPowerSpectrumBase* PL) const
{
   double value = 0;
   for (auto diagram : _twoLoop) {
      value += diagram->value(mom, kernels, PL);
   }
   return value;
}

//------------------------------------------------------------------------------
inline void DiagramSetBase::set_qmax(double qmax)
{
   for (auto diagram : _oneLoop) {
      diagram->set_qmax(qmax);
   }
   for (auto diagram : _twoLoop) {
      diagram->set_qmax(qmax);
   }
}


#endif // DIAGRAM_SET_BASE_HPP

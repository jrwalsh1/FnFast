//------------------------------------------------------------------------------
/// \file DiagramTree.hpp
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
//    Interface of class DiagramTree
//------------------------------------------------------------------------------

#ifndef DIAGRAM_TREE_HPP
#define DIAGRAM_TREE_HPP

#include "DiagramBase.hpp"

using namespace std;

class DiagramTree : public DiagramBase
{
   private:
      Order _order;     ///< order of the calculation

   public:
      /// constructor
      DiagramTree(vector<Line> lines);
      /// destructor
      virtual ~DiagramTree() {}

      /// returns the diagram value with the input momentum routing
      double value(const MomentumMap<ThreeVector>& mom, const VertexMap<KernelBase*>& kernels, LinearPowerSpectrumBase* PL) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // DIAGRAM_TREE_HPP
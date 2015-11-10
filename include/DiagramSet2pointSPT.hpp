//------------------------------------------------------------------------------
/// \file DiagramSet2pointSPT.hpp
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
//    Definition of base class DiagramSet2pointSPT
//------------------------------------------------------------------------------

#ifndef DIAGRAM_SET_2_POINT_SPT_HPP
#define DIAGRAM_SET_2_POINT_SPT_HPP

#include "DiagramSetBase.hpp"

namespace fnfast {

class DiagramSet2pointSPT : public DiagramSetBase
{
   private:
      LabelMap<Graphs_2point, DiagramBase*> _diagrams;      ///< container for diagrams
      LabelMap<Vertex, VertexType> _vertextypes;            ///< vertex types
      LabelMap<Vertex, KernelType> _kerneltypes;            ///< kernel types

   public:
      /// constructor
      DiagramSet2pointSPT(Order order, LabelMap<Vertex, KernelType> kerneltypes = {{Vertex::v1, KernelType::delta}, {Vertex::v2, KernelType::delta}});
      /// destructor
      virtual ~DiagramSet2pointSPT() {}

      /// access diagrams
      DiagramBase* operator[](const Graphs_2point& graph) { return _diagrams[graph]; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // DIAGRAM_SET_2_POINT_SPT_HPP

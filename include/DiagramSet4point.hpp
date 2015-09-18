//------------------------------------------------------------------------------
/// \file DiagramSet4point.hpp
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
//    Definition of base class DiagramSet4point
//------------------------------------------------------------------------------

#ifndef DIAGRAM_SET_4_POINT_HPP
#define DIAGRAM_SET_4_POINT_HPP

#include "DiagramSetBase.hpp"

namespace fnfast {

class DiagramSet4point : public DiagramSetBase
{
   public:
      /// graph labels
      enum class Graphs_4point : int {
         // tree
         T3111,
         T2211,
         // one loop
         T5111,
         T4211a,
         T4211b,
         T3311a,
         T3311b,
         T3221a,
         T3221b,
         T3221c,
         T2222
      };

   private:
      std::map<Graphs_4point, DiagramBase*> _diagrams;      ///< container for diagrams

   public:
      /// constructor
      DiagramSet4point(Order order);
      /// destructor
      virtual ~DiagramSet4point() {}

      /// access diagrams
      DiagramBase* operator[](const Graphs_4point& graph) { return _diagrams[graph]; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // DIAGRAM_SET_4_POINT_HPP

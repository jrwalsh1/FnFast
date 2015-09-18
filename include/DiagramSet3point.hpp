//------------------------------------------------------------------------------
/// \file DiagramSet3point.hpp
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
//    Definition of base class DiagramSet3point
//------------------------------------------------------------------------------

#ifndef DIAGRAM_SET_3_POINT_HPP
#define DIAGRAM_SET_3_POINT_HPP

#include "DiagramSetBase.hpp"

namespace fnfast {

class DiagramSet3point : public DiagramSetBase
{
   public:
      /// graph labels
      enum class Graphs_3point : int {
         // tree
         B211,
         // one loop
         B411,
         B321a,
         B321b,
         B222
      };

   private:
      std::map<Graphs_3point, DiagramBase*> _diagrams;      ///< container for diagrams

   public:
      /// constructor
      DiagramSet3point(Order order);
      /// destructor
      virtual ~DiagramSet3point() {}

      /// access diagrams
      DiagramBase* operator[](const Graphs_3point& graph) { return _diagrams[graph]; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // DIAGRAM_SET_3_POINT_HPP

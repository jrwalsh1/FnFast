//------------------------------------------------------------------------------
/// \file DiagramSet2point.hpp
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
//    Definition of base class DiagramSet2point
//------------------------------------------------------------------------------

#ifndef DIAGRAM_SET_2_POINT_HPP
#define DIAGRAM_SET_2_POINT_HPP

#include "DiagramSetBase.hpp"

using namespace std;

class DiagramSet2point : public DiagramSetBase
{
   public:
      /// graph labels
      enum class Graphs_2point : int {
         P11,
         P31,
         P22,
         P51,
         P42,
         P33a,
         P33b
      };

   private:
      map<Graphs_2point, DiagramBase*> _diagrams;      ///< container for diagrams

   public:
      /// constructor
      DiagramSet2point(Order order);
      /// destructor
      virtual ~DiagramSet2point() {}

      /// access diagrams
      DiagramBase* operator[](const Graphs_2point& graph) const { return NULL; } //_diagrams[graph]; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // DIAGRAM_SET_2_POINT_HPP

//------------------------------------------------------------------------------
/// \file Line.hpp
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
//    Interface of class Line
//------------------------------------------------------------------------------

#ifndef LINE_HPP
#define LINE_HPP

#include "Propagator.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
/**
 * \class Line
 *
 * \brief class for lines in diagram
 *
 * Line(int index_start, int index_end, Propagator prop)
 *
 * Provides public data for the following quantities:
 * - index of the starting and ending vertices
 * - propagator
 */
//------------------------------------------------------------------------------

class Line
{
   public:
      Vertex start;           ///< label of the starting vertex
      Vertex end;             ///< label of the ending vertex
      Propagator propagator;  ///< propagator

      /// constructor
      Line(Vertex start, Vertex end, Propagator prop);
      /// destructor
      virtual ~Line() {}
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

inline Line::Line(Vertex vstart, Vertex vend, Propagator prop) :
start(vstart), end(vend), propagator(prop) {}

} // namespace fnfast

#endif // LINE_HPP

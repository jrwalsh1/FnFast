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

//------------------------------------------------------------------------------
/**
 * \class Line
 *
 * \brief class for lines in diagram
 *
 * Line(int index_start, int index_end, Propagator prop)
 *
 * Provides accessors to the following quantities:
 * - index of the starting and ending vertices
 * - propagator
 */
//------------------------------------------------------------------------------

class Line
{
   private:
      Vertices::VertexLabel _start;    ///< label of the starting vertex
      Vertices::VertexLabel _end;      ///< label of the ending vertex
      Propagator _prop;                ///< propagator

   public:
      /// constructor
      Line(Vertices::VertexLabel start, Vertices::VertexLabel end, Propagator prop);
      /// destructor
      virtual ~Line() {}
   
      /// accessors
      Vertices::VertexLabel start() { return _start; }
      Vertices::VertexLabel end() { return _end; }
      Propagator propagator() { return _prop; }
      
      /// set the starting index
      void set_vertex_start(Vertices::VertexLabel start) { _start = start; }

      /// set the ending index
      void set_vertex_end(Vertices::VertexLabel end) { _end = end; }

      /// set the momentum
      void set_propagator(Propagator prop) { _prop = prop; }
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

inline Line::Line(Vertices::VertexLabel start, Vertices::VertexLabel end, Propagator prop) :
_start(start), _end(end), _prop(prop) {}

#endif // LINE_HPP

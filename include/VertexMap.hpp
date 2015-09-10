//------------------------------------------------------------------------------
/// \file VertexMap.hpp
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
//    Interface of class VertexMap
//------------------------------------------------------------------------------

#ifndef VERTEX_MAP_HPP
#define VERTEX_MAP_HPP

#include <vector>
#include <unordered_map>
#include <initializer_list>

#include "Labels.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class VertexMap
 *
 * \brief templated class acting as a map from VertexLabels to arbitrary containers.
 *
 * VertexMap<T>(map<VertexLabel, T>) : builds map given input argument
 *    (also accepts initializer list)
 *
 * Generic container for objects associated with vertex labels
 * Tracks active labels
 *
 * Provides functions (via labels) for:
 * - Accessing objects
 * - Modifying objects
 * - Retreival of objects
 */
//------------------------------------------------------------------------------

template <typename T>
class VertexMap
{
   private:
      unordered_map<VertexLabel, T> _map;        ///< the map from labels to objects
      vector<VertexLabel> _labels;               ///< active labels

   public:
      /// constructors
      VertexMap();
      VertexMap(unordered_map<VertexLabel, T> vx_map);
      VertexMap(initializer_list<pair<VertexLabel, T> > vx_map);
      /// destructor
      virtual ~VertexMap() {}

      /// size of the map
      size_t size() const;

      /// get allowed labels
      vector<VertexLabel> labels() const;

      /// test if a given label is present
      bool hasLabel(const VertexLabel& label) const;

      /// reset the underlying map, including allowed labels
      void set_map(unordered_map<VertexLabel, T> vx_map);

      /// the subscript operator
      const T& operator[](const VertexLabel& label) const;
      T& operator[](VertexLabel& label);

      /// permute objects using a map on the labels
      void permute(unordered_map<VertexLabel, VertexLabel> perm);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // VERTEX_MAP_HPP

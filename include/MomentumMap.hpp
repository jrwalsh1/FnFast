//------------------------------------------------------------------------------
/// \file MomentumMap.hpp
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
//    Interface of class MomentumMap
//------------------------------------------------------------------------------

#ifndef MOMENTUM_MAP_HPP
#define MOMENTUM_MAP_HPP

#include <vector>
#include <unordered_map>
#include <initializer_list>

#include "Labels.hpp"

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class MomentumMap
 *
 * \brief templated class acting as a map from MomentumLabels to arbitrary containers.
 *
 * MomentumMap<T>(map<MomentumLabel, T>) : builds map given input argument
 *    (also accepts initializer list)
 *
 * Generic container for objects associated with momentum labels
 * Tracks active labels
 *
 * Provides functions (via labels) for:
 * - Accessing objects
 * - Modifying objects
 * - Retreival of objects
 */
//------------------------------------------------------------------------------

template <typename T>
class MomentumMap
{
   private:
      unordered_map<MomentumLabel, T> _map;        ///< the map from labels to objects
      vector<MomentumLabel> _labels;               ///< active labels

   public:
      /// constructors
      MomentumMap();
      MomentumMap(unordered_map<MomentumLabel, T> mom_map);
      MomentumMap(initializer_list<pair<MomentumLabel, T> > mom_map);
      /// destructor
      virtual ~MomentumMap() {}

      /// size of the map
      size_t size() const;

      /// get allowed labels
      vector<MomentumLabel> labels() const;

      /// test if a given label is present
      bool hasLabel(const MomentumLabel& label) const;

      /// reset the underlying map, including allowed labels
      void set_map(unordered_map<MomentumLabel, T> mom_map);

      /// the subscript operator
      const T& operator[](const MomentumLabel& label) const;
      T& operator[](const MomentumLabel& label);

      /// permute objects using a map on the labels
      void permute(MomentumMap<MomentumLabel> perm);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // MOMENTUM_MAP_HPP

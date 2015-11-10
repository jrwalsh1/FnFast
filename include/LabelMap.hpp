//------------------------------------------------------------------------------
/// \file LabelMap.hpp
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
//    Interface of class LabelMap
//------------------------------------------------------------------------------

#ifndef LABEL_MAP_HPP
#define LABEL_MAP_HPP

#include <vector>
#include <unordered_map>
#include <initializer_list>
#include <cassert>

#include "Labels.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
/**
 * \class LabelMap
 *
 * \brief templated class acting as a map from arbitrary labels to arbitrary containers.
 *
 * LabelMap<S, T>(map<S, T>) : builds map given input argument
 *    (also accepts initializer list)
 *
 * Generic container for objects associated with labels
 * Labels must be hashable
 * Tracks active labels
 *
 * Provides functions (via labels) for:
 * - Accessing objects
 * - Modifying objects
 * - Retreival of objects
 */
//------------------------------------------------------------------------------

template <typename S, typename T>
class LabelMap
{
   private:
      std::unordered_map<S, T> _map;        ///< the map from labels to objects
      std::vector<S> _labels;               ///< active labels

   public:
      /// constructors
      LabelMap();
      LabelMap(std::unordered_map<S, T> label_map);
      LabelMap(std::initializer_list<std::pair<S, T> > label_map);
      /// destructor
      virtual ~LabelMap() {}

      /// size of the map
      size_t size() const;

      /// get allowed labels
      std::vector<S> labels() const;

      /// test if a given label is present
      bool hasLabel(const S& label) const;

      /// reset the underlying map, including allowed labels
      void set_map(std::unordered_map<S, T> label_map);

      /// the subscript operator
      const T& operator[](const S& label) const;
      T& operator[](const S& label);

      /// permute objects using a map on the labels
      void permute(LabelMap<S, S> perm);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
template <typename S, typename T>
inline LabelMap<S, T>::LabelMap() {}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline LabelMap<S, T>::LabelMap(std::unordered_map<S, T> label_map)
: _map(label_map) {
   // save the keys from the map
   for (auto item : _map) {
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline LabelMap<S, T>::LabelMap(std::initializer_list<std::pair<S, T> > label_map)
{
   // loop over the elements, add them to the map, save the labels
   for (auto item : label_map) {
      _map[item.first] = item.second;
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline size_t LabelMap<S, T>::size() const
{
   return _map.size();
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline std::vector<S> LabelMap<S, T>::labels() const
{
   return _labels;
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline bool LabelMap<S, T>::hasLabel(const S& label) const
{
   return (_map.find(label) != _map.end());
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline void LabelMap<S, T>::set_map(std::unordered_map<S, T> label_map)
{
   _map = label_map;
   // save the keys from the map
   _labels.clear();
   for (auto item : _map) {
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline const T& LabelMap<S, T>::operator[](const S& label) const
{
   // throw an exception if there is no key found
   if (_map.find(label) == _map.end())
      assert(false && "LabelMap::[] : key not found in map! Quitting.");
   return _map.at(label);
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline T& LabelMap<S, T>::operator[](const S& label)
{
   // throw an exception if there is no key found
   if (_map.find(label) == _map.end())
      assert(false && "LabelMap::[] : key not found in map! Quitting.");
   return _map[label];
}

//------------------------------------------------------------------------------
template <typename S, typename T>
inline void LabelMap<S, T>::permute(LabelMap<S, S> perm)
{
   // copy the underlying map, permute based on that
   std::unordered_map<S, T> curr_map = _map;
   for (auto label : perm.labels()) {
      _map[label] = curr_map[perm[label]];
   }
}

} // namespace fnfast

#endif // LABEL_MAP_HPP

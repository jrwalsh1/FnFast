//------------------------------------------------------------------------------
/// \file MomentumMap.cpp
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
//    Implementation of class MomentumMap
//------------------------------------------------------------------------------

#include <cassert>
#include <iostream>

#include "MomentumMap.hpp"
#include "ThreeVector.hpp"
#include "Propagator.hpp"

//------------------------------------------------------------------------------
template <typename T>
MomentumMap<T>::MomentumMap() {}

//------------------------------------------------------------------------------
template <typename T>
MomentumMap<T>::MomentumMap(unordered_map<MomentumLabel, T> mom_map)
: _map(mom_map) {
   // save the keys from the map
   for (auto item : _map) {
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename T>
MomentumMap<T>::MomentumMap(initializer_list<pair<MomentumLabel, T> > mom_map)
{
   // loop over the elements, add them to the map, save the labels
   for (auto item : mom_map) {
      _map[item.first] = item.second;
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename T>
size_t MomentumMap<T>::size() const
{
   return _map.size();
}

//------------------------------------------------------------------------------
template <typename T>
vector<MomentumLabel> MomentumMap<T>::labels() const
{
   return _labels;
}

//------------------------------------------------------------------------------
template <typename T>
bool MomentumMap<T>::hasLabel(const MomentumLabel& label) const
{
   return (_map.find(label) != _map.end());
}

//------------------------------------------------------------------------------
template <typename T>
void MomentumMap<T>::set_map(unordered_map<MomentumLabel, T> mom_map)
{
   _map = mom_map;
}

//------------------------------------------------------------------------------
template <typename T>
const T& MomentumMap<T>::operator[](const MomentumLabel& label) const
{
   // throw an exception if there is no key found
   if (_map.find(label) == _map.end()) {
      cout << "MomentumMap::[] : key not found in map! Quitting." << endl;
      assert(false);
   }
   return _map.at(label);
}

//------------------------------------------------------------------------------
template <typename T>
T& MomentumMap<T>::operator[](MomentumLabel& label)
{
   // throw an exception if there is no key found
   if (_map.find(label) == _map.end()) {
      cout << "MomentumMap::[] : key not found in map! Quitting." << endl;
      assert(false);
   }
   return _map[label];
}

//------------------------------------------------------------------------------
template <typename T>
void MomentumMap<T>::permute(MomentumMap<MomentumLabel> perm)
{
   // copy the underlying map, permute based on that
   unordered_map<MomentumLabel, T> curr_map = _map;
   for (auto label : perm.labels()) {
      _map[label] = curr_map[perm[label]];
   }
}

template class MomentumMap<MomentumLabel>;
template class MomentumMap<ThreeVector>;
template class MomentumMap<Propagator::LabelFlow>;

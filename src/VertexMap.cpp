//------------------------------------------------------------------------------
/// \file VertexMap.cpp
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
//    Implementation of class VertexMap
//------------------------------------------------------------------------------

#include <cassert>
#include <iostream>

#include "VertexMap.hpp"
#include "Propagator.hpp"
#include "KernelBase.hpp"

//------------------------------------------------------------------------------
template <typename T>
VertexMap<T>::VertexMap() {}

//------------------------------------------------------------------------------
template <typename T>
VertexMap<T>::VertexMap(unordered_map<VertexLabel, T> vx_map)
: _map(vx_map) {
   // save the keys from the map
   for (auto item : _map) {
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename T>
VertexMap<T>::VertexMap(initializer_list<pair<VertexLabel, T> > vx_map)
{
   // loop over the elements, add them to the map, save the labels
   for (auto item : vx_map) {
      _map[item.first] = item.second;
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename T>
size_t VertexMap<T>::size() const
{
   return _map.size();
}

//------------------------------------------------------------------------------
template <typename T>
vector<VertexLabel> VertexMap<T>::labels() const
{
   return _labels;
}

//------------------------------------------------------------------------------
template <typename T>
bool VertexMap<T>::hasLabel(const VertexLabel& label) const
{
   return (_map.find(label) != _map.end());
}

//------------------------------------------------------------------------------
template <typename T>
void VertexMap<T>::set_map(unordered_map<VertexLabel, T> vx_map)
{
   _map = vx_map;
   // save the keys from the map
   _labels.clear();
   for (auto item : _map) {
      _labels.push_back(item.first);
   }
}

//------------------------------------------------------------------------------
template <typename T>
const T& VertexMap<T>::operator[](const VertexLabel& label) const
{
   // throw an exception if there is no key found
   assert(_map.find(label) != _map.end());
   return _map.at(label);
}

//------------------------------------------------------------------------------
template <typename T>
T& VertexMap<T>::operator[](VertexLabel& label)
{
   // throw an exception if there is no key found
   assert(_map.find(label) != _map.end());
   return _map[label];
}

//------------------------------------------------------------------------------
template <typename T>
void VertexMap<T>::permute(unordered_map<VertexLabel, VertexLabel> perm)
{
   // copy the underlying map, permute based on that
   unordered_map<VertexLabel, T> curr_map = _map;
   for (auto item : perm) {
      _map[item.first] = curr_map[item.second];
   }
}

template class VertexMap<vector<Propagator> >;
template class VertexMap<KernelBase*>;
template class VertexMap<int>;

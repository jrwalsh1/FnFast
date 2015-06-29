//------------------------------------------------------------------------------
/// \file Labels.hpp
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
//    Definition of Order, MomentumLabel constants
//------------------------------------------------------------------------------

#ifndef LABELS_HPP
#define LABELS_HPP

#include <vector>
#include <functional>

using namespace std;

//------------------------------------------------------------------------------
/**
 * \enum Order
 *
 * \brief Defines constants to label the order of a diagram / calculation.
 *
 * One of kTree or kOneLoop.
 */
//------------------------------------------------------------------------------
enum Order
{
   kTree  = 0,
   kOneLoop = 1
};

//------------------------------------------------------------------------------
/**
 * \struct Vertices
 *
 * \brief Defines constants to label the vertices and a vector of all values
 *
 * Current labels are v1, v2, v3, v4.
 * If desired, more slots for vertices may be added.
 */
//------------------------------------------------------------------------------
struct Vertices
{
   enum VertexLabel
   {
      v1 = 1,
      v2 = 2,
      v3 = 3,
      v4 = 4
   };

   static const vector<VertexLabel> vertexlabels;
};

//------------------------------------------------------------------------------
/**
 * \struct VertexPair
 *
 * \brief A container for a pair of VertexLabels
 *
 * Defines equality operator for use in maps, comparisons
 */
//------------------------------------------------------------------------------
struct VertexPair
{
   Vertices::VertexLabel vA;     ///< label A
   Vertices::VertexLabel vB;     ///< label B

   /// constructor
   VertexPair(Vertices::VertexLabel vxA, Vertices::VertexLabel vxB) : vA(vxA), vB(vxB) {}

   /// equality operator
   bool operator==(const VertexPair& rhs) const {
      if (((vA == rhs.vA) && (vB == rhs.vB)) || ((vA == rhs.vB) && (vB == rhs.vA))) { return true; }
      else { return false; }
   }

   /// comparison operator
   bool operator<(const VertexPair& rhs) const {
      if (min(vA, vB) < min(rhs.vA, rhs.vB)) {
         return true;
      } else if ((min(vA, vB) == min(rhs.vA, rhs.vB)) && (max(vA, vB) < max(rhs.vA, rhs.vB))) {
         return true;
      }
      return false;
   }
};

//------------------------------------------------------------------------------
/**
 * \struct MomentumLabels
 *
 * \brief Defines constants to label the momenta and a vector of all values
 *
 * Current labels are q (loop momentum), k1, k2, k3, k4 (external momenta).
 * If desired, more slots for loop and external momenta may be added.
 */
//------------------------------------------------------------------------------
struct Momenta{
   enum MomentumLabel
   {
      q  = 0,
      k1 = 1,
      k2 = 2,
      k3 = 3,
      k4 = 4
   };

   static const vector<MomentumLabel> momentumlabels;
};

/// hash functions for label objects
namespace std
{
   /// hash function for VertexLabel
   template <>
   struct hash<Vertices::VertexLabel>
   {
      size_t operator()(const Vertices::VertexLabel &v) const
      {
         // Compute individual hash values for VertexLabel
         return ((hash<int>()(v)) >> 1);
      }
   };

   /// hash function for VertexPair
   template <>
   struct hash<VertexPair>
   {
      size_t operator()(const VertexPair& vp) const
      {
         // Compute individual hash values for two data members and combine them using XOR and bit shifting
         return ((hash<int>()(vp.vA) ^ (hash<int>()(vp.vB) << 1)) >> 1);
      }
   };

   /// hash function for MomentumLabel
   template <>
   struct hash<Momenta::MomentumLabel>
   {
      size_t operator()(const Momenta::MomentumLabel& k) const
      {
         // Compute individual hash values for MomentumLabel
         return ((hash<int>()(k)) >> 1);
      }
   };
}


#endif // LABELS_HPP

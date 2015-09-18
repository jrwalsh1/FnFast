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
//    Definition of Order, Momentum, Vertex constants
//------------------------------------------------------------------------------

#ifndef LABELS_HPP
#define LABELS_HPP

#include <vector>
#include <functional>

namespace fnfast {

//------------------------------------------------------------------------------
/**
 * \enum class Order
 *
 * \brief Defines constants to label the order of a diagram / calculation.
 *
 * One of kTree, kOneLoop, kTwoLoop.
 */
//------------------------------------------------------------------------------
enum class Order
{
   kTree,
   kOneLoop,
   kTwoLoop
};

//------------------------------------------------------------------------------
/**
 * \enum class Vertex
 *
 * \brief Defines constants to label the vertices
 *
 * Current labels are v1, v2, v3, v4.
 * If desired, more slots for vertices may be added.
 */
//------------------------------------------------------------------------------
enum class Vertex : int
{
   v1 = 1,
   v2 = 2,
   v3 = 3,
   v4 = 4
};

//------------------------------------------------------------------------------
/**
 * \struct VertexPair
 *
 * \brief A container for a pair of Vertex labels
 *
 * Defines equality operator for use in maps, comparisons
 */
//------------------------------------------------------------------------------
struct VertexPair
{
   Vertex vA;     ///< label A
   Vertex vB;     ///< label B

   /// constructor
   VertexPair(Vertex vxA, Vertex vxB) : vA(vxA), vB(vxB) {}

   /// equality operator
   bool operator==(const VertexPair& rhs) const {
      if (((vA == rhs.vA) && (vB == rhs.vB)) || ((vA == rhs.vB) && (vB == rhs.vA))) { return true; }
      else { return false; }
   }

   /// comparison operator
   bool operator<(const VertexPair& rhs) const {
      if (std::min(vA, vB) < std::min(rhs.vA, rhs.vB)) {
         return true;
      } else if ((std::min(vA, vB) == std::min(rhs.vA, rhs.vB)) && (std::max(vA, vB) < std::max(rhs.vA, rhs.vB))) {
         return true;
      }
      return false;
   }
};

//------------------------------------------------------------------------------
/**
 * \enum class Momentum
 *
 * \brief Defines constants to label the momenta
 *
 * Current labels are q, q2 (loop momenta), k1, k2, k3, k4 (external momenta).
 * q2 should only be used in two-loop diagrams (as the 2nd loop momentum),
 * q should be used as the loop momentum in one-loop diagrams.
 * If desired, more slots for loop and external momenta may be added.
 */
//------------------------------------------------------------------------------
enum class Momentum : int
{
   q2 = -1,
   q = 0,
   k1 = 1,
   k2 = 2,
   k3 = 3,
   k4 = 4
};

} // namespace fnfast

/// hash functions for label objects
//------------------------------------------------------------------------------
/// hash function for Vertex
template <>
struct std::hash<fnfast::Vertex>
{
   size_t operator()(const fnfast::Vertex &v) const
   {
      // Compute individual hash values for Vertex
      return ((std::hash<int>()(static_cast<int>(v))) >> 1);
   }
};

/// hash function for VertexPair
template <>
struct std::hash<fnfast::VertexPair>
{
   size_t operator()(const fnfast::VertexPair& vp) const
   {
      // Compute individual hash values for two data members and combine them using XOR and bit shifting
      return ((std::hash<int>()(static_cast<int>(vp.vA)) ^ (std::hash<int>()(static_cast<int>(vp.vB)) << 1)) >> 1);
   }
};

/// hash function for Momentum
template <>
struct std::hash<fnfast::Momentum>
{
   size_t operator()(const fnfast::Momentum& k) const
   {
      // Compute individual hash values for Momentum
      return ((std::hash<int>()(static_cast<int>(k))) >> 1);
   }
};

#endif // LABELS_HPP

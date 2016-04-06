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
 * \enum class VertexType
 *
 * \brief Defines dummy types for vertices
 *
 * Current labels are type1, type2, type3, type4.
 * If desired, more slots for vertex types may be added.
 */
//------------------------------------------------------------------------------
enum class VertexType : int
{
   type1 = 1,
   type2 = 2,
   type3 = 3,
   type4 = 4
};

//------------------------------------------------------------------------------
/**
 * \enum class KernelType
 *
 * \brief Defines types for kernels in N-point functions, delta or theta
 *
 * labels are delta, theta
 */
//------------------------------------------------------------------------------
enum class KernelType : int
{
   delta,
   theta
};

//------------------------------------------------------------------------------
/**
 * \struct VertexObjectPair
 *
 * \brief A container for a pair of complete vertex objects
 *
 * Includes the vertex label, vertex type, and kernel type
 * for both vertices in the pair
 */
//------------------------------------------------------------------------------
struct VertexObjectPair
{
   Vertex vertexA;            ///< vertex label for A
   Vertex vertexB;            ///< vertex label for B
   VertexType vertexAtype;    ///< vertex type for A
   VertexType vertexBtype;    ///< vertex type for B
   KernelType kernelAtype;    ///< kernel type for A
   KernelType kernelBtype;    ///< kernel type for B

   /// constructor
   VertexObjectPair(Vertex vxA, Vertex vxB, VertexType vxAtype, VertexType vxBtype, KernelType kAtype, KernelType kBtype)
   : vertexA(vxA), vertexB(vxB), vertexAtype(vxAtype), vertexBtype(vxBtype), kernelAtype(kAtype), kernelBtype(kBtype) {}

   /// equality operator
   bool operator==(const VertexObjectPair& rhs) const {
      if ((((vertexA == rhs.vertexA) && (vertexB == rhs.vertexB)) || ((vertexA == rhs.vertexB) && (vertexB == rhs.vertexA)))
            && (vertexAtype == rhs.vertexAtype) && (vertexBtype == rhs.vertexBtype)
            && (kernelAtype == rhs.kernelAtype) && (kernelBtype == rhs.kernelBtype))
               { return true; }
      else { return false; }
   }

   /// comparison operator: duplicate operator< for the underlying VertexPair
   bool operator<(const VertexObjectPair& rhs) const {
      VertexPair vxpair(vertexA, vertexB);
      VertexPair vxpairRHS(rhs.vertexA, rhs.vertexB);
      return (vxpair < vxpairRHS);
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

//------------------------------------------------------------------------------
/**
 * \enum class Graphs_2point
 *
 * \brief Defines constants to label 2-point diagrams
 *
 * Defines constants to label 2-point diagrams
 */
//------------------------------------------------------------------------------
enum class Graphs_2point : int {
   // ---------- SPT graph labels ----------
   // tree
   P11,
   // one loop
   P31,
   P22,
   // two loop
   P51,
   P42,
   P33a,
   P33b,
   // ---------- EFT graph labels ----------
   // one loop counterterms
   P31x,
   // two loop counterterms
   P51x,
   P42x,
   P33ax
};

//------------------------------------------------------------------------------
/**
 * \enum class Graphs_3point
 *
 * \brief Defines constants to label 3-point diagrams
 *
 * Defines constants to label 3-point diagrams
 */
//------------------------------------------------------------------------------
enum class Graphs_3point : int {
   // ---------- SPT graph labels ----------
   // tree
   B211,
   // one loop
   B411,
   B321a,
   B321b,
   B222,
   // ---------- EFT graph labels ----------
   B411x,
   B321ax
};

//------------------------------------------------------------------------------
/**
 * \enum class Graphs_4point
 *
 * \brief Defines constants to label 4-point diagrams
 *
 * Defines constants to label 4-point diagrams
 */
//------------------------------------------------------------------------------
enum class Graphs_4point : int {
   // ---------- SPT graph labels ----------
   // tree
   T3111,
   T2211,
   // one loop
   T5111,
   T4211a,
   T4211b,
   T3311a,
   T3311b,
   T3221a,
   T3221b,
   T3221c,
   T2222,
   // ---------- EFT graph labels ----------
   T5111x,
   T4211ax,
   T3311ax,
   T3221ax
};

} // namespace fnfast

namespace std {

/// hash functions for label objects
//------------------------------------------------------------------------------
/// hash function for Vertex
template <>
struct hash<fnfast::Vertex>
{
   size_t operator()(const fnfast::Vertex &v) const
   {
      // Compute individual hash values for Vertex
      return ((hash<int>()(static_cast<int>(v))) >> 1);
   }
};

/// hash function for VertexPair
template <>
struct hash<fnfast::VertexPair>
{
   size_t operator()(const fnfast::VertexPair& vp) const
   {
      // Compute individual hash values for two data members and combine them using XOR and bit shifting
      return ((hash<int>()(static_cast<int>(vp.vA)) ^ (hash<int>()(static_cast<int>(vp.vB)) << 1)) >> 1);
   }
};

/// hash function for Momentum
template <>
struct hash<fnfast::Momentum>
{
   size_t operator()(const fnfast::Momentum& k) const
   {
      // Compute individual hash values for Momentum
      return ((hash<int>()(static_cast<int>(k))) >> 1);
   }
};

/// hash function for Graphs_2point
template <>
struct hash<fnfast::Graphs_2point>
{
   size_t operator()(const fnfast::Graphs_2point& k) const
   {
      // Compute individual hash values for Graphs_2point
      return ((hash<int>()(static_cast<int>(k))) >> 1);
   }
};

/// hash function for Graphs_3point
template <>
struct hash<fnfast::Graphs_3point>
{
   size_t operator()(const fnfast::Graphs_3point& k) const
   {
      // Compute individual hash values for Graphs_2point
      return ((hash<int>()(static_cast<int>(k))) >> 1);
   }
};

/// hash function for Graphs_4point
template <>
struct hash<fnfast::Graphs_4point>
{
   size_t operator()(const fnfast::Graphs_4point& k) const
   {
      // Compute individual hash values for Graphs_2point
      return ((hash<int>()(static_cast<int>(k))) >> 1);
   }
};

} // namespace std

#endif // LABELS_HPP

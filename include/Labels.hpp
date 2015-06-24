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

#endif // LABELS_HPP

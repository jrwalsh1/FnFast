//------------------------------------------------------------------------------
/// \file Labels.cpp
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
//    Implementation of objects in Vertices and Momenta structs
//------------------------------------------------------------------------------

#include "Labels.hpp"

//------------------------------------------------------------------------------
const vector<Vertices::VertexLabel> Vertices::vertexlabels = {v1, v2, v3, v4};

//------------------------------------------------------------------------------
const vector<Momenta::MomentumLabel> Momenta::momentumlabels = {q, k1, k2, k3, k4};

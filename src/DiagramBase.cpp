//------------------------------------------------------------------------------
/// \file DiagramBase.cpp
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
//    Implementation of base class DiagramBase
//------------------------------------------------------------------------------

#include <iostream>
#include <limits>
#include <algorithm>

#include "DiagramBase.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
DiagramBase::DiagramBase(std::vector<Line> lines) : _lines(lines)
{
   // construct the vertex momenta map
   std::unordered_map<Vertex, std::vector<Propagator> > vx_momenta;
   // iterate over lines, fill in other diagram objects
   for (auto line : _lines) {
      // for the start and end vertices, add the momentum to the vertex factor list
      vx_momenta[line.start].push_back(line.propagator);
      vx_momenta[line.end].push_back(line.propagator.reverse());
      // store the vertex pairs for each line
      _vertexpairs.push_back(VertexPair(line.start, line.end));
   }
   _vertexmomenta.set_map(vx_momenta);

   // assumes the momenta and vertices are canonically ordered,
   // e.g. {k1, k2, k3} and {v1, v2, v3} for a 3-point graph
   // and not {k2, k3, k4} or {v2, v3, v4}
   size_t nvertices = _vertexmomenta.labels().size();
   for (size_t i = 1; i <= nvertices; i++) {
      _vertices.push_back(static_cast<Vertex>(i));
      _extmomlabels.push_back(static_cast<Momentum>(i));
   }

   // create a map of constant vertex type (since it was not passed as an argument)
   std::unordered_map<Vertex, VertexType> vxtypes;
   for (auto vertex : _vertices) {
      vxtypes[vertex] = VertexType::type1;
   }
   _vertextypes.set_map(vxtypes);

   // create a map of constant kernel type
   std::unordered_map<Vertex, KernelType> ktypes;
   for (auto vertex : _vertices) {
      ktypes[vertex] = KernelType::delta;
   }
   _kerneltypes.set_map(ktypes);

   // calculate the symmetry factor
   _symfac = calc_symmetry_factor();

   // calculate the permutations of the external momenta
   _perms = calc_permutations();
}

//------------------------------------------------------------------------------
DiagramBase::DiagramBase(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes) : _lines(lines), _vertextypes(vertextypes)
{
   // construct the vertex momenta map
   std::unordered_map<Vertex, std::vector<Propagator> > vx_momenta;
   // iterate over lines, fill in other diagram objects
   for (auto line : _lines) {
      // for the start and end vertices, add the momentum to the vertex factor list
      vx_momenta[line.start].push_back(line.propagator);
      vx_momenta[line.end].push_back(line.propagator.reverse());
      // store the vertex pairs for each line
      _vertexpairs.push_back(VertexPair(line.start, line.end));
   }
   _vertexmomenta.set_map(vx_momenta);

   // assumes the momenta and vertices are canonically ordered,
   // e.g. {k1, k2, k3} and {v1, v2, v3} for a 3-point graph
   // and not {k2, k3, k4} or {v2, v3, v4}
   size_t nvertices = _vertexmomenta.labels().size();
   for (size_t i = 1; i <= nvertices; i++) {
      _vertices.push_back(static_cast<Vertex>(i));
      _extmomlabels.push_back(static_cast<Momentum>(i));
   }

   // create a map of constant kernel type
   std::unordered_map<Vertex, KernelType> ktypes;
   for (auto vertex : _vertices) {
      ktypes[vertex] = KernelType::delta;
   }
   _kerneltypes.set_map(ktypes);

   // calculate the symmetry factor
   _symfac = calc_symmetry_factor();

   // calculate the permutations of the external momenta
   _perms = calc_permutations();
}

//------------------------------------------------------------------------------
DiagramBase::DiagramBase(std::vector<Line> lines, LabelMap<Vertex, VertexType> vertextypes, LabelMap<Vertex, KernelType> kerneltypes) : _lines(lines), _vertextypes(vertextypes), _kerneltypes(kerneltypes)
{
   // construct the vertex momenta map
   std::unordered_map<Vertex, std::vector<Propagator> > vx_momenta;
   // iterate over lines, fill in other diagram objects
   for (auto line : _lines) {
      // for the start and end vertices, add the momentum to the vertex factor list
      vx_momenta[line.start].push_back(line.propagator);
      vx_momenta[line.end].push_back(line.propagator.reverse());
      // store the vertex pairs for each line
      _vertexpairs.push_back(VertexPair(line.start, line.end));
   }
   _vertexmomenta.set_map(vx_momenta);

   // assumes the momenta and vertices are canonically ordered,
   // e.g. {k1, k2, k3} and {v1, v2, v3} for a 3-point graph
   // and not {k2, k3, k4} or {v2, v3, v4}
   size_t nvertices = _vertexmomenta.labels().size();
   for (size_t i = 1; i <= nvertices; i++) {
      _vertices.push_back(static_cast<Vertex>(i));
      _extmomlabels.push_back(static_cast<Momentum>(i));
   }

   // calculate the symmetry factor
   _symfac = calc_symmetry_factor();

   // calculate the permutations of the external momenta
   _perms = calc_permutations();
}

//------------------------------------------------------------------------------
double DiagramBase::calc_symmetry_factor()
{
   /*
    * calculation of the (internal) symmetry factor
    * this factor is defined as the number of diagrams with an equivalent topology
    * that give the same value for any momentum routing.
    * This function uses the Solon formula, prod_i N_i! / prod_{i,j} P_ij!
    * where N_i is the number of lines from vertex i,
    * and P_ij is the number of lines between vertices i and j
    */
   // count the number of lines for each vertex
   std::unordered_map<Vertex, int> vertexcounts;
   for (auto vertex : _vertices) { vertexcounts[vertex] = 0; }
   for (auto vx_pair : _vertexpairs) {
      vertexcounts[vx_pair.vA]++;
      vertexcounts[vx_pair.vB]++;
   }
   // calculate the numerator
   int numerator = 1;
   for (auto vertex : _vertices) {
      numerator *= factorial(vertexcounts[vertex]);
   }
   // calculate the number of lines between each vertex pair
   std::unordered_map<VertexPair, int> linecounts;
   int denominator = 1;
   // loop over all possible pairs of vertices
   for (size_t i = 0; i < _vertices.size(); i++) {
      for (size_t j = 0; j <= i; j++) {
         VertexPair vxpair(_vertices[i], _vertices[j]);
         linecounts[vxpair] = 0;
         // see how many vertex pairs in the diagram match the given one
         for (auto vx_pair : _vertexpairs) {
            if (vx_pair == vxpair) {
               linecounts[vxpair]++;
               // need to double-add in the case that
               // a line connects a vertex to itself
               if (i == j) { linecounts[vxpair]++; }
            }
         }
         denominator *= factorial(linecounts[vxpair]);
      }
   }
   // the symmetry factor
   double symfac = numerator * 1. / denominator;

   return symfac;
}

//------------------------------------------------------------------------------
std::vector<LabelMap<Momentum, Momentum> > DiagramBase::calc_permutations()
{
   /*
    * Calculation of the external momentum permutations.
    * Returns the set of permutations needed to correctly sum over all
    * external momentum routings of the graph.
    * Eliminates degenerate configurations.
    * The permutations are stored as a map of momentum label mappings.
    */
   std::vector<LabelMap<Momentum, Momentum> > perms;
   int nperm = 0;

   // create a list of VertexPair objects, where each represents
   // a distinct relabeling of the vertices (and hence a permutation of external momenta).
   // We will iterate over vertex labeling permutations and store only those
   // which will give a unique diagram value.
   std::vector<std::vector<VertexObjectPair> > vertexconnections;

   // Loop over vertex permutations.
   // Since we will be applying the same permutation to momentum labels and vertex labels,
   // we index permutations with a simple list of integers
   std::vector<int> indices;
   for (size_t i = 0; i < _vertices.size(); i++) { indices.push_back(i); }
   // make sure it's sorted...
   std::sort(indices.begin(), indices.end());
   // main permutation loop
   do {
      // create a map from the canonical vertex labels to the permuted ones
      std::unordered_map<Vertex, Vertex> vertexmap;
      for (size_t i = 0; i < _vertices.size(); i++) {
         vertexmap[_vertices[i]] = _vertices[indices[i]];
      }

      // use this map to create a vector of VertexObjectPair objects for each line
      // that are the diagram lines under the vertex permutation
      std::vector<VertexObjectPair> vertexpairs;
      for (auto line : _lines) {
         // store the vertex pairs for each line
         Vertex vxstart = vertexmap[line.start];
         Vertex vxend = vertexmap[line.end];
         vertexpairs.push_back(VertexObjectPair(vxstart, vxend, _vertextypes[vxstart], _vertextypes[vxend], _kerneltypes[vxstart], _kerneltypes[vxend]));
      }
      // now sort this object and compare it to existing ones
      // to find unique permutations
      std::sort(vertexpairs.begin(), vertexpairs.end());

      bool newitem = true;
      for (auto vxconn : vertexconnections) {
         if (vertexpairs == vxconn) { newitem = false; break; }
      }

      // if we have encountered a new ordering, add it to the list
      // and save the permutation in terms of the momentum labels map
      if (newitem) {
         nperm++;
         vertexconnections.push_back(vertexpairs);
         // create the map for the momentum labels (as we did for the vertices)
         std::unordered_map<Momentum, Momentum> extmommap;
         for (size_t i = 0; i < _vertices.size(); i++) {
            extmommap[_extmomlabels[i]] = _extmomlabels[indices[i]];
         }
         perms.push_back(LabelMap<Momentum, Momentum>(extmommap));
      }
   } while (std::next_permutation(indices.begin(), indices.end()));

   return perms;
}

} // namespace fnfast

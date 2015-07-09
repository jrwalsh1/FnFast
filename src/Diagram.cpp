//------------------------------------------------------------------------------
/// \file Diagram.cpp
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
//    Implementation of class Diagram
//------------------------------------------------------------------------------

#include <iostream>
#include <limits>

#include "Diagram.hpp"
#include "SPTkernels.hpp"

/*
 * The implementation of a SPT diagram
 * in terms of the lines in the graph
 */

//------------------------------------------------------------------------------
Diagram::Diagram(vector<Line> lines, unordered_map<Vertices::VertexLabel, KernelBase*> kernels, LinearPowerSpectrumBase* PL) : _lines(lines), _kernels(kernels), _PL(PL), _qmax(numeric_limits<double>::infinity())
{
   // iterate over lines, fill in other diagram objects
   _order = kTree;
   for (size_t i = 0; i < _lines.size(); i++) {
      // check if the line has the loop momentum in it
      // if so, it has an IR pole that must be regulated if it is away from 0
      if (_lines[i].propagator().hasLabel(Momenta::q)) {
         _order = kOneLoop;
         Propagator pole = _lines[i].propagator().IRpole(Momenta::q);
         if (!pole.isNull()) {
            _IRpoles.push_back(pole);
         }
      }
      // for the start and end vertices, add the momentum to the vertex factor list
      _vertexmomenta[_lines[i].start()].push_back(_lines[i].propagator());
      _vertexmomenta[_lines[i].end()].push_back(_lines[i].propagator().reverse());
      // store the vertex pairs for each line
      _vertexpairs.push_back(VertexPair(_lines[i].start(), _lines[i].end()));
   }

   // calculate the symmetry factor
   _symfac = calc_symmetry_factor();

   // store the vertex and momentum label list
   // assumes the vertices are canonically ordered,
   // e.g. {v1, v2, v3} for a bispectrum graph and not {v1, v2, v4}
   // same for momenta: {k1, k2, k3} for a bispectrum graph and not {k2, k3, k4}
   size_t nvertices = kernels.size();
   _vertices = vector<Vertices::VertexLabel>(Vertices::vertexlabels.begin(), Vertices::vertexlabels.begin() + nvertices);
   _extmomlabels = vector<Momenta::MomentumLabel>(Momenta::momentumlabels.begin() + 1, Momenta::momentumlabels.begin() + nvertices + 1);

   // calculate the permutations of the external momenta
   _perms = calc_permutations();
}

//------------------------------------------------------------------------------
double Diagram::value_base(DiagramMomenta mom)
{
   // check to see if the loop momentum is above the cutoff, if so return 0
   // first check whether the diagram has a loop or not
   bool has_loop = false;
   for (size_t c = 0; c < mom.labels.size(); c++) {
      if (mom.labels[c] == Momenta::q) {
         has_loop = true;
         break;
      }
   }
   if (has_loop) {
      if (mom[Momenta::q].magnitude() > _qmax) { return 0; }
   }

   // the diagram value is:
   // symmetry factor * propagators * vertices
   double value = _symfac;
   // iterate over lines
   for (size_t i = 0; i < _lines.size(); i++) {
      value *= (*_PL)(_lines[i].propagator().p(mom).magnitude());
   }
   // now do vertex factors
   for (size_t i = 0; i < _vertices.size(); i++) {
      Vertices::VertexLabel vertex = Vertices::vertexlabels[i];
      vector<ThreeVector> p;
      for (size_t j = 0; j < _vertexmomenta[vertex].size(); j++) {
         p.push_back(_vertexmomenta[vertex][j].p(mom));
      }
      value *= _kernels[vertex]->Fn_sym(p);
   }
   return value;
}

//------------------------------------------------------------------------------
double Diagram::value_base_IRreg(DiagramMomenta mom)
{
   // no IR regulation necessary if there are no poles away from q = 0
   if (_IRpoles.empty()) return value_base(mom);

   // To regulate the diagram in the IR, we map each region with
   // an IR pole at q = qIR != 0 onto coordinates with the pole at q = 0
   double value = 0;
   // need to regulate only the unique IR poles
   // e.g. in the covariance limit, two IR poles can be degenerate
   // and we should treat them simultaneously
   vector<ThreeVector> uniqueIRpoles;
   // pole at q = 0
   uniqueIRpoles.push_back(ThreeVector(0, 0, 0));
   // loop over the nonzero poles
   for (size_t i = 0; i < _IRpoles.size(); i++) {
      // check if pole is unique
      bool is_unique = true;
      ThreeVector pole = _IRpoles[i].p(mom);
      for (size_t j = 0; j < uniqueIRpoles.size(); j++) {
         if (pole == uniqueIRpoles[j]) {
            is_unique = false;
            break;
         }
      }
      if (is_unique) { uniqueIRpoles.push_back(pole); }
   }
   // now loop over all the unique IR poles
   for (size_t i = 0; i < uniqueIRpoles.size(); i++) {
      // for these poles we change variables: q -> q + pole
      // so that the pole maps to 0 and we exclude all other poles
      ThreeVector pole = uniqueIRpoles[i];
      double PSregion = 1;
      // loop over all other poles and make PS cuts for each
      for (size_t j = 0; j < uniqueIRpoles.size(); j++) {
         if (j != i) {
            ThreeVector pole_j = uniqueIRpoles[j];
            PSregion *= theta(mom[Momenta::q], mom[Momenta::q] + pole - pole_j);
         }
      }
      // copy and shift the diagram momentum for the pole
      DiagramMomenta mom_shift = mom;
      mom_shift[Momenta::q] = mom[Momenta::q] + pole;
      // add the diagram value for this shifted momentum, times the PS factor
      value += PSregion * value_base(mom_shift);
   }

   return value;
}

//------------------------------------------------------------------------------
double Diagram::value_IRreg(DiagramMomenta mom)
{
   /* 
    * To return the IR regulated diagram symmetrized over external momenta,
    * we symmetrize the IR regulated diagram with the input momentum routing
    * over external momentum configurations.
    * To make the symmetrization more efficient, we compute only the
    * momentum configurations giving distinct diagram values,
    * and multiply each by the appropriate symmetry factor
    */

   double value = 0;
   // loop over external momentum permutations
   // symmetrize over q -> -q
   for (size_t i = 0; i < _perms.size(); i++) {
      DiagramMomenta mom_perm = mom;
      mom_perm.permute(_perms[i]);
      value += 0.5 * value_base_IRreg(mom_perm);
      ThreeVector mq = -1 * mom_perm[Momenta::q];
      mom_perm[Momenta::q] = mq;
      value += 0.5 * value_base_IRreg(mom_perm);
   }

   return value;
}

//------------------------------------------------------------------------------
double Diagram::calc_symmetry_factor()
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
   unordered_map<Vertices::VertexLabel, int> vertexcounts;
   for (size_t i = 0; i < Vertices::vertexlabels.size(); i++) { vertexcounts[Vertices::vertexlabels[i]] = 0; }
   for (size_t i = 0; i < _vertexpairs.size(); i++) {
      vertexcounts[_vertexpairs[i].vA]++;
      vertexcounts[_vertexpairs[i].vB]++;
   }
   // calculate the numerator
   int numerator = 1;
   for (size_t i = 0; i < Vertices::vertexlabels.size(); i++) {
      numerator *= factorial(vertexcounts[Vertices::vertexlabels[i]]);
   }

   // calculate the number of lines between each vertex pair
   unordered_map<VertexPair, int> linecounts;
   int denominator = 1;
   for (size_t i = 0; i < Vertices::vertexlabels.size(); i++) {
      for (size_t j = 0; j <= i; j++) {
         VertexPair vxpair(Vertices::vertexlabels[i], Vertices::vertexlabels[j]);
         linecounts[vxpair] = 0;
         for (size_t k = 0; k < _vertexpairs.size(); k++) {
            if (_vertexpairs[k] == vxpair) {
               linecounts[vxpair]++;
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
vector<unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> > Diagram::calc_permutations()
{
   /*
    * Calculation of the external momentum permutations.
    * Returns the set of permutations needed to correctly sum over all
    * external momentum routings of the graph.
    * Eliminates degenerate configurations.
    * The permutations are stored as a map of momentum label mappings.
    */
   vector<unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> > perms;
   int nperm = 0;

   // create a list of VertexPair objects, where each represents
   // a distinct relabeling of the vertices (and hence a permutation of external momenta).
   // We will iterate over vertex labeling permutations and store only those
   // which will give a unique diagram value.
   vector<vector<VertexPair> > vertexconnections;

   // Loop over vertex permutations.
   // Since we will be applying the same permutation to momentum labels and vertex labels,
   // we index permutations with a simple list of integers
   vector<int> indices;
   for (size_t i = 0; i < _vertices.size(); i++) { indices.push_back(i); }
   // make sure it's sorted...
   sort(indices.begin(), indices.end());
   do {
      // create a map from the canonical vertex labels to the permuted ones
      unordered_map<Vertices::VertexLabel, Vertices::VertexLabel> vertexmap;
      for (size_t i = 0; i < _vertices.size(); i++) {
         vertexmap[_vertices[i]] = _vertices[indices[i]];
      }

      // use this map to create a vector of VertexPair objects for each line
      vector<VertexPair> vertexpairs;
      for (size_t i = 0; i < _lines.size(); i++) {
         // store the vertex pairs for each line
         vertexpairs.push_back(VertexPair(vertexmap[_lines[i].start()], vertexmap[_lines[i].end()]));
      }
      // now sort this object and compare it to existing ones
      sort(vertexpairs.begin(), vertexpairs.end());
//      cout << "-----permutation-----" << endl;
//      for (size_t i = 0; i < vertexpairs.size(); i++) {
//         cout << "vertex pair: " << vertexpairs[i].vA << " , " << vertexpairs[i].vB << endl;
//      }
      bool newitem = true;
      for (size_t i = 0; i < vertexconnections.size(); i++) {
         if (vertexpairs == vertexconnections[i]) { newitem = false; break; }
      }

      // if we have encountered a new ordering, add it to the list
      // and save the permutation in terms of the momentum labels map
      if (newitem) {
         nperm++;
         vertexconnections.push_back(vertexpairs);
         // create the map for the momentum labels (as we did for the vertices)
         unordered_map<Momenta::MomentumLabel, Momenta::MomentumLabel> extmommap;
         for (size_t i = 0; i < _vertices.size(); i++) {
            extmommap[_extmomlabels[i]] = _extmomlabels[indices[i]];
         }
         perms.push_back(extmommap);
      }
   } while (next_permutation(indices.begin(), indices.end()));

   return perms;
}

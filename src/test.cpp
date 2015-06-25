//------------------------------------------------------------------------------
// test bed for EFTofLSS library
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "SPTkernels.hpp"
#include "ThreeVector.hpp"
#include "Random.hpp"
#include "Diagram.hpp"
#include "LinearPowerSpectrumAnalytic.hpp"

using namespace std;

// main routine
int main()
{
   ThreeVector k1(1.2, 0.7, -0.4);
   ThreeVector k2(-1.0, 0.4, 0.1);
   ThreeVector k3(-0.5, -0.9, 0.6);
   ThreeVector k4(0.4, 0.1, 0.8);
   ThreeVector k5(-0.4, -0.1, -0.8);

   SPTkernels kernels;
   
   vector<ThreeVector> mom;
   mom.push_back(k1);
   mom.push_back(k2);
//   double F2 = kernels.Fn(mom);
//   double F2sym = kernels.Fn_sym(mom);
   mom.push_back(k3);
//   double F3 = kernels.Fn(mom);
//   double F3sym = kernels.Fn_sym(mom);
   mom.push_back(k4);
//   double F4 = kernels.Fn(mom);
//   double F4sym = kernels.Fn_sym(mom);
   mom.push_back(k5);
   double F5 = kernels.Fn(mom);
   double F5sym = kernels.Fn_sym(mom);

//   cout << "F2, F2sym = " << F2 << ", " << F2sym << endl;
//   cout << "F3, F3sym = " << F3 << ", " << F3sym << endl;
//   cout << "F4, F4sym = " << F4 << ", " << F4sym << endl;
   cout << "F5, F5sym = " << F5 << ", " << F5sym << endl;

   unordered_map<Momenta::MomentumLabel, ThreeVector> ext_mom {{Momenta::k1, k1}, {Momenta::k2, k2}};
   DiagramMomenta momenta(ext_mom);
   momenta[Momenta::q] = k3;

   unordered_map<Momenta::MomentumLabel, double> mom1 {{Momenta::q, 1}};
   unordered_map<Momenta::MomentumLabel, double> mom2 {{Momenta::q, -1}, {Momenta::k2, 1}};
   Line line1(Vertices::v1, Vertices::v2, mom1);
   Line line2(Vertices::v1, Vertices::v2, mom2);
   vector<Line> P22_lines = {line1, line2};
   unordered_map<Vertices::VertexLabel, KernelBase*> vertexkernels {{Vertices::v1, &kernels}, {Vertices::v2, &kernels}};

   LinearPowerSpectrumAnalytic PL(1);

   Diagram P22(P22_lines, vertexkernels, &PL);

   return 0;
} // end of main program

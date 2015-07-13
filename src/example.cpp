//------------------------------------------------------------------------------
// example functionality of EFTofLSS library
// program exhibits some things you can do with the code
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
#include "PowerSpectrum.hpp"
#include "Bispectrum.hpp"
#include "Trispectrum.hpp"
#include "LinearPowerSpectrumAnalytic.hpp"
#include "LinearPowerSpectrumCAMB.hpp"

using namespace std;

// main routine
int main()
{
   cout << "===== example code for FnFast =====" << endl;

   // create some momenta
   ThreeVector q(-0.1, 0.2, 0.5);
   ThreeVector k1(0.5, 0, 0.1);
   ThreeVector k2(0.2, -0.4, -0.3);
   ThreeVector k3(-0.3, 0.3, -0.6);
   ThreeVector k4 = -k1 + k2;
   cout << "momentum: " << endl;
   cout << "3 vector: " << k4.p1() << ", " << k4.p2() << ", " << k4.p3() << endl;
   cout << "magnitude: " << k4.magnitude() << endl << endl;

   // create an object that will calculate SPT kernels
   SPTkernels kernels;
   // compute an example: F3^s
   vector<ThreeVector> momenta {k1, k2, k3};
   double F3s = kernels.Fn_sym(momenta);
   cout << "F3^s: " << F3s << endl;

   // create linear power spectra
   // power law (in this case P_L(k) = k)
   LinearPowerSpectrumAnalytic PL(1);
   // interpolation of CAMB data file
   LinearPowerSpectrumCAMB PLcamb("data/LIdata.txt");
   // example call
   cout << "PL from CAMB data file: " << PLcamb(k3.magnitude()) << endl << endl;

   // create EFT counterterm coefficients object
   EFTcoefficients coeffs;
   // set a coefficient's value
   coeffs[EFTcoefficients::cs] = 1.;

   // create a SPT diagram, T5111
   // propagators
   Propagator prop_T5111_q(unordered_map<Momenta::MomentumLabel, double> {{Momenta::q, 1}});
   Propagator prop_T5111_k2(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k2, 1}});
   Propagator prop_T5111_k3(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k3, 1}});
   Propagator prop_T5111_k4(unordered_map<Momenta::MomentumLabel, double> {{Momenta::k4, 1}});
   // lines
   Line line_T5111_11(Vertices::v1, Vertices::v1, prop_T5111_q);
   Line line_T5111_12(Vertices::v1, Vertices::v2, prop_T5111_k2);
   Line line_T5111_13(Vertices::v1, Vertices::v3, prop_T5111_k3);
   Line line_T5111_14(Vertices::v1, Vertices::v4, prop_T5111_k4);
   vector<Line> lines_T5111 {line_T5111_11, line_T5111_12, line_T5111_13, line_T5111_14};
   // base kernels object
   KernelBase* kernels_SPT = new SPTkernels;
   // vertex kernels
   unordered_map<Vertices::VertexLabel, KernelBase *> vertex_kernels = {{Vertices::v1, kernels_SPT}, {Vertices::v2, kernels_SPT}, {Vertices::v3, kernels_SPT}, {Vertices::v4, kernels_SPT}};
   // define the diagram
   Diagram T5111(lines_T5111, vertex_kernels, &PLcamb);
   // display the symmetry factor and number of external momentum permutations (including the base one)
   cout << "T5111 symmetry factor, # of external momentum permutations = " << T5111.symmetry_factor() << ", " << T5111.nperms() << endl;
   // compute the value for the base momentum routing and for the fully IR regulated and symmetrized graph, for one point in phase space
   unordered_map<Momenta::MomentumLabel, ThreeVector> mom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::k4, k4}};
   DiagramMomenta diagmom(mom);
   cout << "T5111 base value, IR regulated and symmetrized value = " << T5111.value_base(diagmom) << ", " << T5111.value_IRreg(diagmom) << endl << endl;

   // create a trispectrum object at 1 loop
   Trispectrum TS(kOneLoop, &PL, &coeffs);
   // set the maximum loop momentum magnitude
   TS.set_qmax(2.);
   // set the random number seed for VEGAS
   TS.set_seed(37);
   // compute the 1-loop SPT trispectrum, integrated over the loop momentum
   double k1mag = k1.magnitude();
   double k2mag = k2.magnitude();
   double costheta12 = (k1*k2) / (k1mag*k2mag);
   IntegralResult TSresult = TS.cov_loopSPT(k1mag, k2mag, costheta12);
   cout << "1-loop SPT trispectrum result = " << TSresult.result << endl << endl;

   // create a power spectrum object at 1 loop
   PowerSpectrum PS(kOneLoop, &PL, &coeffs);
   // create a bispectrum object at 1 loop
   Bispectrum BS(kOneLoop, &PL, &coeffs);

   return 0;
} // end of main program

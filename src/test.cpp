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
#include "PowerSpectrum.hpp"
#include "Bispectrum.hpp"
#include "Trispectrum.hpp"
#include "LinearPowerSpectrumAnalytic.hpp"
#include "LinearPowerSpectrumCAMB.hpp"

using namespace std;

// main routine
int main()
{
   LinearPowerSpectrumAnalytic PL(1);
   LinearPowerSpectrumCAMB PLcamb("data/LIdata.txt");
   EFTcoefficients coeffs;

   vector<Momenta::MomentumLabel> labels = {Momenta::k1, Momenta::k2, Momenta::q};
   double kmag = 0.2;
   ThreeVector ktest(0, 0, kmag);
   ThreeVector qtest(0, kmag, 0);
   unordered_map<Momenta::MomentumLabel, ThreeVector> momentamap = {{Momenta::k1, -ktest}, {Momenta::k2, ktest}, {Momenta::q, qtest}};
   DiagramMomenta pdiag(labels);
   pdiag.set_momenta(momentamap);

   PowerSpectrum PS(kOneLoop, &PLcamb, &coeffs);

   SPTkernels kernels;
   vector<ThreeVector> pvec {qtest, -qtest, ktest};
   vector<ThreeVector> pvec2 {-qtest, qtest, ktest};
   cout << "F3: " << kernels.Fn_sym(pvec) << endl;
   cout << "F3: " << kernels.Fn_sym(pvec2) << endl;

   double P31value = PS[PowerSpectrum::P31]->value_IRreg(pdiag);
   cout << "P31 value = " << P31value << endl;
   double P22value = PS[PowerSpectrum::P22]->value_IRreg(pdiag);
   cout << "P22 value = " << P22value << endl;

   double result = PS.oneLoopSPT_value(kmag);
   cout << "1 loop SPT PS result = " << result << endl;

   return 0;
} // end of main program

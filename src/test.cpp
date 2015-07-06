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

   PowerSpectrum PS(kOneLoop, &PL, &coeffs);
   Bispectrum BS(kOneLoop, &PL, &coeffs);
   Trispectrum TS(kOneLoop, &PL, &coeffs);

   SPTkernels kernels;
   vector<ThreeVector> pvec {qtest, -qtest, ktest};
   vector<ThreeVector> pvec2 {-qtest, qtest, ktest};
   cout << "F3: " << kernels.Fn_sym(pvec) << endl;
   cout << "F3: " << kernels.Fn_sym(pvec2) << endl;

   double P31value = PS[PowerSpectrum::P31]->value_IRreg(pdiag);
   cout << "P31 value = " << P31value << endl;
   double P22value = PS[PowerSpectrum::P22]->value_IRreg(pdiag);
   cout << "P22 value = " << P22value << endl;

   double qmag = 1.1;
   double qcostheta = 0.7;
   double qphi = 2.9;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);
   ThreeVector k1(0.5, 0, 0.1);
   ThreeVector k2(-0.3, 0.3, 0.4);
   ThreeVector k3 = -k1 - k2;
   double k1mag = k1.magnitude();
   double k2mag = k2.magnitude();
   double theta12 = acos((k1*k2) / (k1mag*k2mag));
   unordered_map<Momenta::MomentumLabel, ThreeVector> BSmom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3}};


   double B411_value = BS[Bispectrum::B411]->value_base(BSmom);
   double B411_valuefull = BS[Bispectrum::B411]->value_IRreg(BSmom);
   cout << "B411 value = " << B411_value << endl;
   cout << "B411 total value = " << B411_valuefull << endl;

   double B321a_value = BS[Bispectrum::B321a]->value_base(BSmom);
   double B321a_valuefull = BS[Bispectrum::B321a]->value_IRreg(BSmom);
   cout << "B321a value = " << B321a_value << endl;
   cout << "B321a total value = " << B321a_valuefull << endl;

   double B321b_value = BS[Bispectrum::B321b]->value_base(BSmom);
   double B321b_valuefull = BS[Bispectrum::B321b]->value_IRreg(BSmom);
   cout << "B321b value = " << B321b_value << endl;
   cout << "B321b total value = " << B321b_valuefull << endl;

   double B222_value = BS[Bispectrum::B222]->value_base(BSmom);
   double B222_valuefull = BS[Bispectrum::B222]->value_IRreg(BSmom);
   cout << "B222 value = " << B222_value << endl;
   cout << "B222 total value = " << B222_valuefull << endl;


//   double PSresult = PS.loopSPT(kmag);
//   cout << "1 loop SPT PS result = " << PSresult << endl;
   
   double BSresult = BS.loopSPT(k1mag, k2mag, theta12);
   cout << "1 loop SPT BS result = " << BSresult << endl;

//   cout << "k1, k3, theta13 = " << k1mag << ", " << k3mag << ", " << theta13 << endl;
//   double TSresult = TS.cov_loopSPT(k1mag, k3mag, theta13);
//   cout << "1 loop SPT TS result = " << TSresult << endl;

   return 0;
} // end of main program

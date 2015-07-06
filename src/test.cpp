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

   double qmag = 0.4;
   double qcostheta = 0.4;
   double qphi = 0.5;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);
   ThreeVector k1(0, 0, 0.1);
   ThreeVector k2 = -k1;
   ThreeVector k3(0.2 * sqrt(1 - 0.3*0.3), 0, 0.2 * 0.3);
   ThreeVector k4 = -k3;
   double k1mag = k1.magnitude();
   double k3mag = k3.magnitude();
   double theta13 = acos((k1*k3) / (k1mag*k3mag));
   unordered_map<Momenta::MomentumLabel, ThreeVector> TSmom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::k4, k4}};

   double T5111_value = TS[Trispectrum::T5111]->value_base(TSmom);
   double T5111_valuefull = TS[Trispectrum::T5111]->value_IRreg(TSmom);
   cout << "T5111 value = " << T5111_value << endl;
   cout << "T5111 total value = " << T5111_valuefull << endl;

   double T4211a_value = TS[Trispectrum::T4211a]->value_base(TSmom);
   double T4211a_valuefull = TS[Trispectrum::T4211a]->value_IRreg(TSmom);
   cout << "T4211a value = " << T4211a_value << endl;
   cout << "T4211a total value = " << T4211a_valuefull << endl;

   double T4211b_value = TS[Trispectrum::T4211b]->value_base(TSmom);
   double T4211b_valuefull = TS[Trispectrum::T4211b]->value_IRreg(TSmom);
   cout << "T4211b value = " << T4211b_value << endl;
   cout << "T4211b total value = " << T4211b_valuefull << endl;

   double T3311a_value = TS[Trispectrum::T3311a]->value_base(TSmom);
   double T3311a_valuefull = TS[Trispectrum::T3311a]->value_IRreg(TSmom);
   cout << "T3311a value = " << T3311a_value << endl;
   cout << "T3311a total value = " << T3311a_valuefull << endl;

   double T3311b_value = TS[Trispectrum::T3311b]->value_base(TSmom);
   double T3311b_valuefull = TS[Trispectrum::T3311b]->value_IRreg(TSmom);
   cout << "T3311b value = " << T3311b_value << endl;
   cout << "T3311b total value = " << T3311b_valuefull << endl;

   double T3221a_value = TS[Trispectrum::T3221a]->value_base(TSmom);
   double T3221a_valuefull = TS[Trispectrum::T3221a]->value_IRreg(TSmom);
   cout << "T3221a value = " << T3221a_value << endl;
   cout << "T3221a total value = " << T3221a_valuefull << endl;

   double T3221b_value = TS[Trispectrum::T3221b]->value_base(TSmom);
   double T3221b_valuefull = TS[Trispectrum::T3221b]->value_IRreg(TSmom);
   cout << "T3221b value = " << T3221b_value << endl;
   cout << "T3221b total value = " << T3221b_valuefull << endl;

   double T3221c_value = TS[Trispectrum::T3221c]->value_base(TSmom);
   double T3221c_valuefull = TS[Trispectrum::T3221c]->value_IRreg(TSmom);
   cout << "T3221c value = " << T3221c_value << endl;
   cout << "T3221c total value = " << T3221c_valuefull << endl;

   double T2222_value = TS[Trispectrum::T2222]->value_base(TSmom);
   double T2222_valuefull = TS[Trispectrum::T2222]->value_IRreg(TSmom);
   cout << "T2222 value = " << T2222_value << endl;
   cout << "T2222 total value = " << T2222_valuefull << endl;

   cout << "=========================================" << endl;

   cout << "T5111 value = " << T5111_value << endl;
   cout << "T5111 total value = " << T5111_valuefull << endl;

   cout << "T4211a value = " << T4211a_value << endl;
   cout << "T4211a total value = " << T4211a_valuefull << endl;

   cout << "T4211b value = " << T4211b_value << endl;
   cout << "T4211b total value = " << T4211b_valuefull << endl;

   cout << "T3311a value = " << T3311a_value << endl;
   cout << "T3311a total value = " << T3311a_valuefull << endl;

   cout << "T3311b value = " << T3311b_value << endl;
   cout << "T3311b total value = " << T3311b_valuefull << endl;

   cout << "T3221a value = " << T3221a_value << endl;
   cout << "T3221a total value = " << T3221a_valuefull << endl;

   cout << "T3221b value = " << T3221b_value << endl;
   cout << "T3221b total value = " << T3221b_valuefull << endl;

   cout << "T3221c value = " << T3221c_value << endl;
   cout << "T3221c total value = " << T3221c_valuefull << endl;

   cout << "T2222 value = " << T2222_value << endl;
   cout << "T2222 total value = " << T2222_valuefull << endl;


   /*
   double PSresult = PS.loopSPT(kmag);
   cout << "1 loop SPT PS result = " << PSresult << endl;
   
   double BSresult = BS.loopSPT(kmag, kmag, 0.5);
   cout << "1 loop SPT BS result = " << BSresult << endl;
   */

   cout << "k1, k3, theta13 = " << k1mag << ", " << k3mag << ", " << theta13 << endl;
   double TSresult = TS.cov_loopSPT(k1mag, k3mag, theta13);
   cout << "1 loop SPT TS result = " << TSresult << endl;

   return 0;
} // end of main program

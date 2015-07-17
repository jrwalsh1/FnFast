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

   PowerSpectrum PS(kOneLoop, &PL, &coeffs);
   Bispectrum BS(kOneLoop, &PL, &coeffs);
   Trispectrum TS(kOneLoop, &PLcamb, &coeffs);
   TS.set_qmax(12.);
   TS.set_seed(37);

   double qmag = 0.01;
   double qcostheta = 0.9;
   double qphi = 0;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);

   ThreeVector k1(0.5, 0, 0.1);
   ThreeVector k2(0, 0, 0.4);
   ThreeVector k3 = -k1;
   ThreeVector k4 = -k2;
   double k1mag = k1.magnitude();
   double k2mag = k2.magnitude();
   double costheta12 = (k1*k2) / (k1mag*k2mag);

   cout << "---------- computing 1-loop trispectrum ----------" << endl;
   unordered_map<Momenta::MomentumLabel, ThreeVector> TSmom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3}, {Momenta::k4, k4}};
//   cout << "T3311b poles: " << endl;
//   TS[Trispectrum::T3311b]->value_IRreg(TSmom);


   IntegralResult TSresult = TS.cov_loopSPT(k1mag, k2mag, costheta12);
   cout << "1 loop SPT TS result = " << TSresult.result << endl;

   
/*
   cout << "---------- computing 1-loop power spectrum ----------" << endl;
   IntegralResult PSresult = PS.loopSPT(k1mag);
   cout << "1 loop SPT PS result = " << PSresult.result << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing 1-loop bispectrum ----------" << endl;
   IntegralResult BSresult = BS.loopSPT(k1mag, k2mag, costheta12);
   cout << "1 loop SPT BS result = " << BSresult.result << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing 1-loop trispectrum ----------" << endl;
   IntegralResult TSresult = TS.cov_loopSPT(k1mag, k2mag, costheta12);
   cout << "1 loop SPT TS result = " << TSresult.result << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing tree level trispectrum ----------" << endl;
   IntegralResult TSresult_tree = TS.cov_tree(k1mag, k2mag);
   cout << "tree-level TS result = " << TSresult_tree.result << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing 1-loop trispectrum, angular integrated ----------" << endl;
   IntegralResult TSresult_loop = TS.cov_loopSPT(k1mag, k2mag);
   cout << "1-loop TS result = " << TSresult_loop.result << endl;
   cout << "press any key to continue" << endl;
   cin.get();
*/

   return 0;
} // end of main program

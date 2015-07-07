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
   Trispectrum TS(kOneLoop, &PL, &coeffs);

   double qmag = 0.01;
   double qcostheta = 0.9;
   double qphi = 0;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);

   ThreeVector k1(0.5, 0, 0.1);
   ThreeVector k2(-0.3, 0.3, 0.4);
   double k1mag = k1.magnitude();
   double k2mag = k2.magnitude();
   double costheta12 = (k1*k2) / (k1mag*k2mag);

   cout << "---------- computing 1-loop power spectrum ----------" << endl;
   double PSresult = PS.loopSPT(k1mag);
   cout << "1 loop SPT PS result = " << PSresult << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing 1-loop bispectrum ----------" << endl;
   double BSresult = BS.loopSPT(k1mag, k2mag, costheta12);
   cout << "1 loop SPT BS result = " << BSresult << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing 1-loop trispectrum ----------" << endl;
   double TSresult = TS.cov_loopSPT(k1mag, k2mag, costheta12);
   cout << "1 loop SPT TS result = " << TSresult << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing tree level trispectrum ----------" << endl;
   double TSresult_tree = TS.cov_tree(k1mag, k2mag);
   cout << "tree-level TS result = " << TSresult_tree << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   cout << "---------- computing 1-loop trispectrum, angular integrated ----------" << endl;
   double TSresult_loop = TS.cov_loop(k1mag, k2mag);
   cout << "1-loop TS result = " << TSresult_loop << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   return 0;
} // end of main program

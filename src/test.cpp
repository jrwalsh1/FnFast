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

//   PowerSpectrum PS(kOneLoop, &PL, &coeffs);
   Bispectrum BS(kOneLoop, &PL, &coeffs);
//   Trispectrum TS(kOneLoop, &PL, &coeffs);

   double qmag = 1.1;
   double qcostheta = 0.7;
   double qphi = 2.9;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);
   ThreeVector k1(0, 0, 0.5);
   ThreeVector k2(0, 0.4, -0.3);
   ThreeVector k3 = -k1 - k2;
   double k1mag = k1.magnitude();
   double k2mag = k2.magnitude();
   double theta12 = acos((k1*k2) / (k1mag*k2mag));
   unordered_map<Momenta::MomentumLabel, ThreeVector> BSmom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, k2}, {Momenta::k3, k3}};

//   double PSresult = PS.loopSPT(kmag);
//   cout << "1 loop SPT PS result = " << PSresult << endl;

   /*
   unordered_map<Momenta::MomentumLabel, ThreeVector> PSmom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, -k1}};
   DiagramMomenta PSdiagmom(PSmom);
   double PSresult = PS[PowerSpectrum::P22]->value_base_IRreg(PSdiagmom);
   cout << "PS value = " << PSresult << endl;
   */

/*
   for (int c = 0; c < 100; c++) {
      qmag = pow(10, -4 + 0.04*c);
      for (int i = 0; i < 200; i++) {
         qphi = 2. * 3.14159265358979 * 0.005*i;
         q = ThreeVector(0, qmag * sin(qphi), qmag * cos(qphi));
         double BSintegrand = BS.loopSPT_excl(k1, k2, q);
         cout << qmag << "   " << qphi << "   " << qmag*qmag * BSintegrand << endl;
      }
   }
*/
/*
   for (int c = 0; c <= 100; c++) {
      qmag = pow(10, -4 + 0.04*c);
      for (int i = 0; i <= 200; i++) {
         qcostheta = 0.01*i - 1.;
         q = ThreeVector(0, qmag * sqrt(1 - qcostheta*qcostheta), qmag * qcostheta);
         unordered_map<Momenta::MomentumLabel, ThreeVector> PSmom = {{Momenta::q, q}, {Momenta::k1, k1}, {Momenta::k2, -k1}};
         DiagramMomenta PSdiagmom(PSmom);
//         double PSintegrand = PS[PowerSpectrum::P22]->value_base(PSdiagmom);
         double PSintegrand = PS.loopSPT_excl(k1, q);
         cout << qmag << "   " << qcostheta << "   " << qmag*qmag * PSintegrand << endl;
      }
   }
*/


   double BSresult = BS.loopSPT(k1mag, k2mag, theta12);
   cout << "1 loop SPT BS result = " << BSresult << endl;

//   cout << "k1, k3, theta13 = " << k1mag << ", " << k3mag << ", " << theta13 << endl;
//   double TSresult = TS.cov_loopSPT(k1mag, k3mag, theta13);
//   cout << "1 loop SPT TS result = " << TSresult << endl;

   return 0;
} // end of main program

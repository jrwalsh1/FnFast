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
#include "LinearPowerSpectrumAnalytic.hpp"
#include "LinearPowerSpectrumCAMB.hpp"

using namespace std;

// main routine
int main()
{
   LinearPowerSpectrumAnalytic PL(1);
   LinearPowerSpectrumCAMB PLcamb("data/LIdata.txt");

   PowerSpectrum PS(Order::kOneLoop);
   PS.set_qmax(12.);
   PS.set_seed(37);

   double qmag = 0.01;
   double qcostheta = 0.9;
   double qphi = 0;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);

   ThreeVector k1(0.5, 0, 0.1);
   double k1mag = k1.magnitude();
   cout << "k1: " << k1 << endl;
   cout << "q: " << q << endl;

   cout << "---------- computing 1-loop power spectrum ----------" << endl;
   SPTkernels* kernelsSPT = new SPTkernels();
   VertexMap<KernelBase*> kernels_SPT {{VertexLabel::v1, kernelsSPT}, {VertexLabel::v2, kernelsSPT}};
   IntegralResult PSresult = PS.oneLoop(k1mag, kernels_SPT, &PL);
   cout << "1 loop SPT PS result = " << PSresult.result << endl;
   cout << "press any key to continue" << endl;
   cin.get();

   return 0;
} // end of main program

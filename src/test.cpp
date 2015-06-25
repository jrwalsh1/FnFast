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
   
   vector<ThreeVector> mom = {k1, k2, k3, k4, k5};
   double F5 = kernels.Fn(mom);
   double F5sym = kernels.Fn_sym(mom);
   cout << "F5, F5sym = " << F5 << ", " << F5sym << endl;

   LinearPowerSpectrumAnalytic PL(1);

   PowerSpectrum PS(kOneLoop, &PL);
   Bispectrum BS(kOneLoop, &PL);
   Bispectrum TS(kOneLoop, &PL);

   return 0;
} // end of main program

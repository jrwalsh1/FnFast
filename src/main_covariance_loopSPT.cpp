//------------------------------------------------------------------------------
// FnFast
// main program to compute the 1-loop SPT corrections to the covariance matrix
// given the magnitudes of the two momenta
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Trispectrum.hpp"
#include "LinearPowerSpectrumAnalytic.hpp"
#include "LinearPowerSpectrumCAMB.hpp"

using namespace std;

// main routine
int main(int argc, const char* argv[])
{
   // get arguments
   if (argc != 4) {
      cout << "Usage: bin/main_covariance_loopSPT k(magnitude) kp(magnitude) output_file" << endl;
      exit(1);
   }
   double k = atof(argv[1]);
   double kp = atof(argv[2]);
   string outfilename = argv[3];

   // the linear power spectrum
   // LinearPowerSpectrumAnalytic PL(1);
   LinearPowerSpectrumCAMB PLcamb("data/LIdata.txt");

   // the counterterms (default is all coefficients = 0)
   EFTcoefficients coeffs;
   // set up the trispectrum calculation
   Trispectrum TS(kOneLoop, &PLcamb, &coeffs);

   // do the calculation
   double TSresult_loopSPT = TS.cov_loopSPT(k, kp);
   
   // save the result
   ofstream outfile(outfilename.c_str());
   outfile.precision(30);
   outfile << fixed << TSresult_loopSPT << endl;
   outfile.close();

   return 0;
} // end of main program

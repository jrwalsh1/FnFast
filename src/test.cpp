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
   Trispectrum TS(kOneLoop, &PL);

   Diagram* P11 = PS[PowerSpectrum::P11];
   Diagram* P31 = PS[PowerSpectrum::P31];
   Diagram* P22 = PS[PowerSpectrum::P22];
   cout << "symmetry factors: " << endl;
   cout << "P11: " << P11->symmetry_factor() << endl;
   cout << "P31: " << P31->symmetry_factor() << endl;
   cout << "P22: " << P22->symmetry_factor() << endl;

   Diagram* B211 = BS[Bispectrum::B211];
   Diagram* B411 = BS[Bispectrum::B411];
   Diagram* B321a = BS[Bispectrum::B321a];
   Diagram* B321b = BS[Bispectrum::B321b];
   Diagram* B222 = BS[Bispectrum::B222];
   cout << "symmetry factors: " << endl;
   cout << "B211: " << B211->symmetry_factor() << endl;
   cout << "B411: " << B411->symmetry_factor() << endl;
   cout << "B321a: " << B321a->symmetry_factor() << endl;
   cout << "B321b: " << B321b->symmetry_factor() << endl;
   cout << "B222: " << B222->symmetry_factor() << endl;

   Diagram* T3111 = TS[Trispectrum::T3111];
   Diagram* T2211 = TS[Trispectrum::T2211];
   Diagram* T5111 = TS[Trispectrum::T5111];
   Diagram* T4211a = TS[Trispectrum::T4211a];
   Diagram* T4211b = TS[Trispectrum::T4211b];
   Diagram* T3311a = TS[Trispectrum::T3311a];
   Diagram* T3311b = TS[Trispectrum::T3311b];
   Diagram* T3221a = TS[Trispectrum::T3221a];
   Diagram* T3221b = TS[Trispectrum::T3221b];
   Diagram* T3221c = TS[Trispectrum::T3221c];
   Diagram* T2222 = TS[Trispectrum::T2222];
   cout << "symmetry factors: " << endl;
   cout << "T3111: " << T3111->symmetry_factor() << endl;
   cout << "T2211: " << T2211->symmetry_factor() << endl;
   cout << "T5111: " << T5111->symmetry_factor() << endl;
   cout << "T4211a: " << T4211a->symmetry_factor() << endl;
   cout << "T4211b: " << T4211b->symmetry_factor() << endl;
   cout << "T3311a: " << T3311a->symmetry_factor() << endl;
   cout << "T3311b: " << T3311b->symmetry_factor() << endl;
   cout << "T3221a: " << T3221a->symmetry_factor() << endl;
   cout << "T3221b: " << T3221b->symmetry_factor() << endl;
   cout << "T3221c: " << T3221c->symmetry_factor() << endl;
   cout << "T2222: " << T2222->symmetry_factor() << endl;

   cout << "-----------" << endl;

   return 0;
} // end of main program

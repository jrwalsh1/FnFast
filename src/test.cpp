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
#include "EFTkernels.hpp"
#include "ThreeVector.hpp"
#include "PowerSpectrum.hpp"
#include "Bispectrum.hpp"
#include "Covariance.hpp"
#include "LinearPowerSpectrumAnalytic.hpp"
#include "LinearPowerSpectrumCAMB.hpp"
#include "LabelMap.hpp"

using namespace fnfast;

// main routine
int main()
{
   LinearPowerSpectrumAnalytic PL(1);
   LinearPowerSpectrumCAMB PLcamb("data/LIdata.txt");

   PowerSpectrum PS(Order::kOneLoop);
   PS.set_qmax(12.);
   PS.set_seed(37);
   Bispectrum BS(Order::kOneLoop);
   BS.set_qmax(12.);
   BS.set_seed(37);
   Covariance CV(Order::kOneLoop);
   CV.set_qmax(12.);
   CV.set_seed(37);

   double qmag = 0.8;
   double qcostheta = 0.9;
   double qphi = 0;
   ThreeVector q(qmag * sqrt(1 - qcostheta*qcostheta) * cos(qphi), qmag * sqrt(1 - qcostheta*qcostheta) * sin(qphi), qmag * qcostheta);

   ThreeVector k1(0.5, 0, 0.1);
   double k1mag = k1.magnitude();
   ThreeVector k2(0.2, 0.4, -0.3);
   double k2mag = k2.magnitude();
   double theta12 = acos(k1*k2 / (k1.magnitude() * k2.magnitude()));

/*
   std::cout << "---------- computing 1-loop power spectrum ----------" << std::endl;
   SPTkernels* kernelsSPT = new SPTkernels();
   LabelMap<Vertex, KernelBase*> kernels_SPT {{Vertex::v1, kernelsSPT}, {Vertex::v2, kernelsSPT}};
   IntegralResult PSresult = PS.oneLoop(k1mag, kernels_SPT, &PL);
   std::cout << "1 loop SPT PS result = " << PSresult.result << std::endl;
   std::cout << "press any key to continue" << std::endl;
   std::cin.get();
*/

/*
    std::cout << "------------ computing 1-loop EFT power spectrum ------------" << std::endl;
    SPTkernels* kernelsSPT = new SPTkernels();
    LabelMap<Vertex, KernelBase*> kernels_SPT {{Vertex::v1, kernelsSPT}, {Vertex::v2, kernelsSPT}};

    EFTcoefficients coeffs; //All coefficients are initialized to zero
    std::cout<<coeffs.description()<<std::endl; //Description of EFT operators/coefficients
    coeffs[EFTcoefficients::cs]= 10; //Set a value for speed of sound
    EFTkernels* kernelsEFT = new EFTkernels(coeffs);
    LabelMap<Vertex, KernelBase*> kernels_EFT {{Vertex::v1, kernelsEFT}, {Vertex::v2, kernelsSPT}};
    double PSresultEFT = PS.treeEFT(k1mag, kernels_EFT, &PL);
    std::cout << "EFT PS result = " << PSresultEFT << std::endl;
    std::cout << "press any key to continue" << std::endl;
    std::cin.get();
*/

/*
   std::cout << "------------ computing 1-loop bispectrum ------------" << std::endl;
   SPTkernels* kernelsSPT = new SPTkernels();
   LabelMap<Vertex, KernelBase*> kernels_SPT_BS {{Vertex::v1, kernelsSPT}, {Vertex::v2, kernelsSPT}, {Vertex::v3, kernelsSPT}};
   double BStree = BS.tree(k1mag, k2mag, theta12, kernels_SPT_BS, &PL);
   std::cout << "tree level SPT BS result = " << BStree << std::endl;
   IntegralResult BSresult = BS.oneLoop(k1mag, k2mag, theta12, kernels_SPT_BS, &PL);
   std::cout << "1 loop SPT BS result = " << BSresult.result << std::endl;
   std::cout << "press any key to continue" << std::endl;
   std::cin.get();
*/


   std::cout << "------------ computing 1-loop EFT bispectrum ------------" << std::endl;
   SPTkernels* kernelsSPT = new SPTkernels();
   LabelMap<Vertex, KernelBase*> kernels_SPT_BS {{Vertex::v1, kernelsSPT}, {Vertex::v2, kernelsSPT}, {Vertex::v3, kernelsSPT}};

   EFTcoefficients coeffs; //All coefficients are initialized to zero by default
   std::cout<<coeffs.description()<<std::endl; //Description of EFT operators/coefficients
   coeffs[EFTcoefficients::cs]= 10;   //Set a value for speed of sound
   coeffs[EFTcoefficients::c1]= -2;   //Set a value for NLO coefficients
   coeffs[EFTcoefficients::c2]= 1.2;  //Set a value for NLO coefficients
   coeffs[EFTcoefficients::c3]= -0.5; //Set a value for NLO coefficients
   coeffs[EFTcoefficients::t2]= 5;    //Set a value for NLO coefficients
   coeffs[EFTcoefficients::t3]= 12;   //Set a value for NLO coefficients
   EFTkernels* kernelsEFT = new EFTkernels(coeffs);
   LabelMap<Vertex, KernelBase*> kernels_EFT_BS {{Vertex::v1, kernelsEFT}, {Vertex::v2, kernelsSPT}, {Vertex::v3, kernelsSPT}};
   double BSresultEFT = BS.treeEFT(k1mag, k2mag, theta12, kernels_EFT_BS, &PL);
   std::cout << "EFT BS result = " << BSresultEFT << std::endl;
   std::cout << "press any key to continue" << std::endl;
   std::cin.get();


/*
   std::cout << "------------ computing 1-loop covariance ------------" << std::endl;
   SPTkernels* kernelsSPT = new SPTkernels();
   LabelMap<Vertex, KernelBase*> kernels_SPT_CV {{Vertex::v1, kernelsSPT}, {Vertex::v2, kernelsSPT}, {Vertex::v3, kernelsSPT}, {Vertex::v4, kernelsSPT}};
   //IntegralResult CVtree = CV.tree(k1mag, k2mag, kernels_SPT_CV, &PL);
   //std::cout << "tree SPT CV result = " << CVtree.result << std::endl;
   IntegralResult CVresult = CV.oneLoop(k1mag, k2mag, kernels_SPT_CV, &PL);
   std::cout << "1 loop SPT CV result = " << CVresult.result << std::endl;
   std::cout << "press any key to continue" << std::endl;
   std::cin.get();
*/

/*
   std::cout << "------------ computing 1-loop EFT covariance ------------" << std::endl;
   SPTkernels* kernelsSPT = new SPTkernels();
   LabelMap<Vertex, KernelBase*> kernels_SPT_CV {{Vertex::v1, kernelsSPT}, {Vertex::v2, kernelsSPT}, {Vertex::v3, kernelsSPT}, {Vertex::v4, kernelsSPT}};
   EFTcoefficients coeffs; //All coefficients are initialized to zero by default
   std::cout<<coeffs.description()<<std::endl; //Description of EFT operators/coefficients
   coeffs[EFTcoefficients::cs]= 10;   //Set a value for speed of sound
   coeffs[EFTcoefficients::c1]= -2;   //Set a value for NLO coefficients
   coeffs[EFTcoefficients::c2]= 1.2;  //Set a value for NLO coefficients
   coeffs[EFTcoefficients::c3]= -0.5; //Set a value for NLO coefficients
   coeffs[EFTcoefficients::t2]= 5;    //Set a value for NLO coefficients
   coeffs[EFTcoefficients::t3]= 12;   //Set a value for NLO coefficients
   coeffs[EFTcoefficients::d1]= 1;    //Set a value for NNLO coefficients
   coeffs[EFTcoefficients::d2]= 2;    //Set a value for NNLO coefficients
   coeffs[EFTcoefficients::d3]= 3;    //Set a value for NNLO coefficients
   coeffs[EFTcoefficients::d4]= 4;    //Set a value for NNLO coefficients
   coeffs[EFTcoefficients::d5]= 5;    //Set a value for NNLO coefficients
   coeffs[EFTcoefficients::d6]= 6;    //Set a value for NNLO coefficients
   coeffs[EFTcoefficients::d7]= 7;    //Set a value for NNLO coefficients
   EFTkernels* kernelsEFT = new EFTkernels(coeffs);
   LabelMap<Vertex, KernelBase*> kernels_EFT_CV {{Vertex::v1, kernelsEFT}, {Vertex::v2, kernelsSPT}, {Vertex::v3, kernelsSPT},{Vertex::v4, kernelsSPT}};
   IntegralResult CVresultEFT = CV.treeEFT(k1mag, k2mag, kernels_EFT_CV, &PL);
   std::cout << "EFT CV result = " << CVresultEFT.result << std::endl;
   std::cout << "press any key to continue" << std::endl;
   std::cin.get();
*/


   return 0;
} // end of main program

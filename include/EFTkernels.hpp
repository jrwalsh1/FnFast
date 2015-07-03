//------------------------------------------------------------------------------
/// \file EFTkernels.hpp
//
// Author(s):
//    Daniele Bertolini
//
// Copyright:
//    Copyright (C) 2015  LBL
//
//    This file is part of the EFTofLSS library. EFTofLSS is distributed under the
//    terms of the GNU General Public License version 3 (GPLv3), see the COPYING
//    file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Definition of EFTkernels
//------------------------------------------------------------------------------

#ifndef EFT_KERNELS_HPP
#define EFT_KERNELS_HPP

#include <numeric>
#include <vector>
#include <map>
#include <string>
#include <sstream>


#include "KernelBase.hpp"
#include "SPTkernels.hpp"
#include "ThreeVector.hpp"


using namespace std;

//------------------------------------------------------------------------------
/**
 * \namespace EFTkernels
 *
 * \brief Defines the EFT kernels.
 *
 */
//------------------------------------------------------------------------------
//EFT counterterms coefficients

class EFTcoefficients {
    
    public:
    
        //default constructor: initialize all coefficients to zero
        EFTcoefficients(){ for(unsigned int i=0; i<_labels.size(); i++) _coeff_value[_labels[i]]=0;}
        virtual ~EFTcoefficients(){}
    
        //Labels for counterterm coefficients
        enum Labels {cs,ch,c1,c2,c3,c4,c5,c6,ch1,ch2,ch3,d1,d2,d3};
    
        //Input/read values
        double& operator[](EFTcoefficients::Labels label) {return _coeff_value[label];}
    
        //Alternative way of inputing values
        void set_coefficients(map<EFTcoefficients::Labels, double> coeff_value) {_coeff_value=coeff_value;}
    
        //General description
        string description();
    
        //Coefficients description
        string description(EFTcoefficients::Labels label);
    
        //Print coefficient values
        void print_all_coefficients(){ for(unsigned int i=0; i<_labels.size(); i++) cout<<_labels_names[_labels[i]]<<" = "<<_coeff_value[_labels[i]]<<endl;}
    
    private:
    
        map<EFTcoefficients::Labels, double> _coeff_value;
    
        vector<EFTcoefficients::Labels> _labels = {cs,ch,c1,c2,c3,c4,c5,c6,ch1,ch2,ch3,d1,d2,d3};
        vector<string> _labels_names = {"cs","ch","c1","c2","c3","c4","c5","c6","ch1","ch2","ch3","d1","d2","d3"};
    
};


//------------------------------------------------------------------------------
//EFT kernels


class EFTkernels : public KernelBase {
   
    public:
    
        EFTkernels(){}
        EFTkernels(EFTcoefficients& coefficients): _coefficients(&coefficients) {}
    
        void set_coefficients(EFTcoefficients& coefficients) {_coefficients=&coefficients;}
    
        //Coefficients
        double cF_1(int n);
        double cF_2(int n);
        double cG_1(int n);
        double cG_2(int n);


        double alpha(ThreeVector p1, ThreeVector p2){return _sptkernels.alpha(p1,p2);}      ///< kernel function alpha
        double beta(ThreeVector p1, ThreeVector p2){return _sptkernels.beta(p1,p2);}        ///< kernel function alpha

        double Fn(vector<ThreeVector> p);                  ///< EFT kernel Fn (q1, ..., qn)
        double Gn(vector<ThreeVector> p);                  ///< EFT kernel Gn (q1, ..., qn)
    
    private:
    
        EFTcoefficients* _coefficients;
        SPTkernels _sptkernels;
    
        //Counterterm shapes
        vector<double> _lo_shapes(ThreeVector p);
        vector<double> _nlo_shapes(ThreeVector p1, ThreeVector p2);
        vector<double> _nnlo_shapes(ThreeVector p1, ThreeVector p2, ThreeVector p3);
    
        //Helper function
        //Returns the dot product between two vectors
        double _dot_product(vector<double> a, vector<double> b);
    
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

#endif // EFT_KERNELS_HPP

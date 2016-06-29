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

namespace fnfast {

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
        enum Labels {cs,c1,c2,c3,t2,t3,d1,d2,d3,d4,d5,d6};

        //Input/read values
        double& operator[](EFTcoefficients::Labels label) {return _coeff_value[label];}

        //Alternative way of inputing values
        void set_coefficients(std::map<EFTcoefficients::Labels, double> coeff_value) {_coeff_value=coeff_value;}

        //General description
        std::string description();

        //Print coefficient values
        void print_all_coefficients(){ for(unsigned int i=0; i<_labels.size(); i++) std::cout<<_labels_names[_labels[i]]<<" = "<<_coeff_value[_labels[i]]<<std::endl;}

    private:
        std::map<EFTcoefficients::Labels, double> _coeff_value;

        std::vector<EFTcoefficients::Labels> _labels = {cs,c1,c2,c3,t2,t3,d1,d2,d3,d4,d5,d6};
        std::vector<std::string> _labels_names = {"cs","c1","c2","c3","t2","t3","d1","d2","d3","d4","d5","d6"};
};

//------------------------------------------------------------------------------
//EFT kernels


class EFTkernels : public KernelBase {
    public:
        EFTkernels(){}
        EFTkernels(EFTcoefficients& coefficients): _coefficients(&coefficients) {}

        void set_coefficients(EFTcoefficients& coefficients) {_coefficients=&coefficients;}

        //Recursion coefficients
        double cF_E(int n);  ///< F-kernels, with theta-Euler-equation operators
        double cG_E(int n);  ///< G-kernels, with theta-Euler-equation operators
        double cF_C(int n);  ///< F-kernels, with Continuity-equation operators
        double cG_C(int n);  ///< G-kernels, with Continuity-equation operators

        double alpha(const ThreeVector& p1, const ThreeVector& p2) {return _sptkernels.alpha(p1,p2);}      ///< kernel function alpha
        double beta(const ThreeVector& p1, const ThreeVector& p2) {return _sptkernels.beta(p1,p2);}        ///< kernel function beta

        ThreeVector alphaOmega(const ThreeVector& p1, const ThreeVector& p2);   ///< vorticity kernel function alpha
        ThreeVector betaOmega(const ThreeVector& p1, const ThreeVector& p2);    ///< vorticity kernel function beta

        double Fn(const std::vector<ThreeVector>& p);                  ///< EFT kernel Fn (q1, ..., qn)
        double Gn(const std::vector<ThreeVector>& p);                  ///< EFT kernel Gn (q1, ..., qn)

        double Fn_sym(const std::vector<ThreeVector>& p);              ///< symmetrized EFT kernel Fn (q1, ..., qn)
        double Gn_sym(const std::vector<ThreeVector>& p);              ///< symmetrized EFT kernel Gn (q1, ..., qn)

    private:
        EFTcoefficients* _coefficients;
        SPTkernels _sptkernels;

        //EFT operators/shapes and stress-tensor
        std::vector<ThreeVector> _LO_shapes(ThreeVector p);
        std::vector<ThreeVector> _NLO_shapes(ThreeVector p1, ThreeVector p2);
        std::vector<ThreeVector> _NNLO_shapes(ThreeVector p1, ThreeVector p2, ThreeVector p3);
        ThreeVector _tau(const std::vector<ThreeVector>& p);

        //Vorticity
        ThreeVector _omega(const std::vector<ThreeVector>& p);

        //Helper function
        //Returns the dot product between a vector of coefficients and a vector of shapes
        ThreeVector _dot_product(std::vector<double> a, std::vector<ThreeVector> b);
};

////////////////////////////////////////////////////////////////////////////////
// Inline Declarations
////////////////////////////////////////////////////////////////////////////////

} // namespace fnfast

#endif // EFT_KERNELS_HPP

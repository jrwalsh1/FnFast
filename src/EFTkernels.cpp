//------------------------------------------------------------------------------
/// \file EFTkernels.cpp
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
//    Implementation of functions in EFTkernels
//------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>

#include "ThreeVector.hpp"
#include "SPTkernels.hpp"
#include "EFTkernels.hpp"

namespace fnfast {

//------------------------------------------------------------------------------
//EFT coefficients
//------------------------------------------------------------------------------

std::string EFTcoefficients::description()
{
    std::stringstream stream;

   stream << "*******************************************"<<std::endl
   << "FnFast includes leading effective operators (i.e. scaling as k^2/k_NL^2) and up to three-field order."<<std::endl
   << "Refer to arXiv:1512.07630 for the basis used here and details about the construction of the EFT operators."<<std::endl
   << "F1^tilde, F2^tilde, and F3^tilde have one, three, and eight independent operators. The operators are chosen as those corresponding to c_s, c^{DeltaDelta}_{1,2,3}, c^{ThetaTheta}_{2,3}, and c^{DeltaDeltaDelta}_{1...6}."<< std::endl
   << "For simplicity, we have renamed here c_s -> cs, c^{DeltaDelta}_{1,2,3} -> c_{1,2,3}, c^{ThetaTheta}_{2,3} -> t_{2,3} and c^{DeltaDeltaDelta}_{1...6} -> d_{1...6}."<<std::endl
   << "*******************************************";
   return stream.str();
}

//------------------------------------------------------------------------------
//EFT kernels
//------------------------------------------------------------------------------
   
//Recursion coefficients
double EFTkernels::cF_E(int n)
{
    return (-2.) / (2 * n * n + 9 * n + 7);
}

double EFTkernels::cG_E(int n)
{
    return (-2*n - 4.) / (2 * n * n + 9 * n + 7);
}
   
double EFTkernels::cF_C(int n)
{
    return (2*n + 5.) / (2 * n * n + 9 * n + 7);
}

double EFTkernels::cG_C(int n)
{
    return 3. / (2 * n * n + 9 * n + 7);
}

//Vorticity EOM kernels
ThreeVector EFTkernels::alphaOmega(const ThreeVector& p1, const ThreeVector& p2)
{
   double eps = 1e-12;
   ThreeVector kernel(0.,0.,0.);
   
   if(p1*p1 > eps) kernel = crossProduct(p2,p1) / (p1*p1);
   
   return kernel;
}
   
ThreeVector EFTkernels::betaOmega(const ThreeVector& p1, const ThreeVector& p2)
{
   double eps = 1e-12;
   ThreeVector kernel(0.,0.,0.);
   
   if(p1*p1 > eps && p2*p2 > eps) kernel = (p2*p2 + 2*p1*p2) / ((p1*p1)*(p2*p2)) * crossProduct(p1,p2);
      
   return kernel;
}

//------------------------------------------------------------------------------
//Build Ftilde, Gtilde kernels in three steps. At each order:
//1. Write shapes = k_i * tau_ij.
//2. Write tau = 1/(1+delta) * shapes.
//3. Write Ftilde, Gtilde kernels.
//------------------------------------------------------------------------------

//shapes = k_i * tau_ij
//LO shapes
std::vector<ThreeVector> EFTkernels::_LO_shapes(ThreeVector p)
{
    return(std::vector<ThreeVector> {p});
}
   
//------------------------------------------------------------------------------
//NLO shapes
std::vector<ThreeVector> EFTkernels::_NLO_shapes(ThreeVector p1, ThreeVector p2)
{
    double eps = 1e-12;
    ThreeVector shape1(0.,0.,0.), shape2(0.,0.,0.), shape3(0.,0.,0.);
    ThreeVector p = p1+p2;

    //Shape1
    shape1=p;

    //Shape2
    if(p1*p1 > eps) shape2 = p*p1 / (p1*p1)* p1;

    //Shape3
    if(p1*p1 > eps && p2*p2 > eps) shape3 = (p*p1)*(p1*p2) / (2*(p1*p1)*(p2*p2)) * p2 + (p*p2)*(p1*p2) / (2*(p1*p1)*(p2*p2)) * p1;

    return(std::vector<ThreeVector> {shape1,shape2,shape3});
}
   
//------------------------------------------------------------------------------
//NNLO shapes
std::vector<ThreeVector> EFTkernels::_NNLO_shapes(ThreeVector p1, ThreeVector p2, ThreeVector p3)
{
   double eps = 1e-12;
   ThreeVector shape1(0.,0.,0.), shape2(0.,0.,0.), shape3(0.,0.,0.), shape4(0.,0.,0.), shape5(0.,0.,0.), shape6(0.,0.,0.), shape7(0.,0.,0.);
   ThreeVector p = p1+p2+p3;

   //Shape1, Shape2, Shape3
   shape1 = p;
   if(p1*p1 > eps) shape2 = p*p1 / (p1*p1)* p1;
   if(p1*p1 > eps && p2*p2 > eps) shape3 = (p*p1)*(p1*p2) / (2*(p1*p1)*(p2*p2)) * p2 + (p*p2)*(p1*p2) / (2*(p1*p1)*(p2*p2)) * p1;

   
   //Shape4
   if(p1*p1 > eps && p2*p2 > eps) shape4 = (p1*p2)*(p1*p2) / ((p1*p1)*(p2*p2)) * p;
   
   //Shape5
   if(p1*p1 > eps && p2*p2 > eps && p3*p3 > eps) shape5 = (p2*p3)*(p2*p3)*(p*p1) / ((p1*p1)*(p2*p2)*(p3*p3)) * p1;
   
   //Shape6
   if(p1*p1 > eps && p2*p2 > eps && p3*p3 > eps) shape6 = (p1*p3)*(p2*p3)*(p*p1) / (2*(p1*p1)*(p2*p2)*(p3*p3)) * p2 + (p1*p3)*(p2*p3)*(p*p2) / (2*(p1*p1)*(p2*p2)*(p3*p3)) * p1;
   
   return(std::vector<ThreeVector> {shape1,shape2,shape3,shape4,shape5,shape6});
}
   
//------------------------------------------------------------------------------
//tau = 1/(1+delta) * shapes
ThreeVector EFTkernels::_tau(const std::vector<ThreeVector>& p)
{
   int n = p.size();
   ThreeVector tau(0.,0.,0.);
   
   double cs =  (*_coefficients)[EFTcoefficients::cs];
   std::vector<double> c = {(*_coefficients)[EFTcoefficients::c1],(*_coefficients)[EFTcoefficients::c2],(*_coefficients)[EFTcoefficients::c3]};
   std::vector<double> t = {0.,(*_coefficients)[EFTcoefficients::t2],(*_coefficients)[EFTcoefficients::t3]};
   std::vector<double> d = {(*_coefficients)[EFTcoefficients::d1],(*_coefficients)[EFTcoefficients::d2],(*_coefficients)[EFTcoefficients::d3],(*_coefficients)[EFTcoefficients::d4],(*_coefficients)[EFTcoefficients::d5],(*_coefficients)[EFTcoefficients::d6]};
   
   
   //Handle trivial cases
   if(n==0 || n>3) { return tau; }
   
   if(n==1) { ThreeVector p1 = p[0]; tau = cs*_LO_shapes(p1)[0]; }
   
   if(n==2) {
      ThreeVector p1 = p[0];
      ThreeVector p2 = p[1];
      
      std::vector<ThreeVector> P2={p2};

      tau = cs*_sptkernels.Fn_sym({p1,p2})*_LO_shapes(p1+p2)[0]+_dot_product(c,_NLO_shapes(p1,p2))+_dot_product(t,_NLO_shapes(p1,p2))-_tau(P2);
   }
   
   if(n==3) {
      ThreeVector p1 = p[0];
      ThreeVector p2 = p[1];
      ThreeVector p3 = p[2];
      
      std::vector<ThreeVector> P23={p2,p3};
      std::vector<ThreeVector> P3={p3};
      
      tau = cs*_sptkernels.Fn_sym({p1,p2,p3})*_LO_shapes(p1+p2+p3)[0]+_dot_product(c,_NLO_shapes(p1,p2+p3))*_sptkernels.Fn_sym({p2,p3})+_dot_product(c,_NLO_shapes(p1+p2,p3))*_sptkernels.Fn_sym({p1,p2})+_dot_product(t,_NLO_shapes(p1,p2+p3))*_sptkernels.Gn_sym({p2,p3})+_dot_product(t,_NLO_shapes(p1+p2,p3))*_sptkernels.Gn_sym({p1,p2})+_dot_product(d,_NNLO_shapes(p1,p2,p3))-_tau(P23)-_sptkernels.Fn_sym({p1,p2})*_tau(P3);
   }
   
   return tau;
}
   
//------------------------------------------------------------------------------
//vorticity
ThreeVector EFTkernels::_omega(const std::vector<ThreeVector>& p)
{
   int n = p.size();
   ThreeVector omega(0.,0.,0.);
   
   //Handle trivial cases
   if (n == 0 || n == 1) { return omega; }
   if (n > 2) { std::cout<<"Vorticity is not available at the order specified"<<std::endl; return omega; }
   
   if (n == 2) { omega = 2. / 9 * crossProduct(p[0]+p[1],_tau(p)); }
   
   return omega;
}

//------------------------------------------------------------------------------
//kernels
double EFTkernels::Fn(const std::vector<ThreeVector>& p)
{
    int n = p.size();
    double Fnval=0;

    //Handle trivial cases
    if (n == 0) { return Fnval; }
    if (n > 3) { std::cout<<"There is no EFT kernel available at the order specified"<<std::endl; return Fnval; }

    if (n==1) { Fnval = cF_E(1) * (p[0]*_tau(p)) ;}

    if (n==2) {
       
       std::vector<ThreeVector> P1={p[0]};
       std::vector<ThreeVector> P2={p[1]};
       
       Fnval = cF_C(2) * alpha(p[0],p[1]) * (Gn_sym(P1) + Fn_sym(P2)) - cF_E(2) * beta(p[0],p[1]) * (Gn_sym(P1) + Gn_sym(P2)) + cF_E(2) * (p[0]+p[1])*_tau(p);
    }

    if (n==3) {
       
       std::vector<ThreeVector> P1={p[0]};
       std::vector<ThreeVector> P23={p[1],p[2]};
       std::vector<ThreeVector> P12={p[0],p[1]};
       std::vector<ThreeVector> P3={p[2]};

       Fnval = cF_C(3) * alpha(p[0],p[1]+p[2]) * (Gn_sym(P1)*_sptkernels.Fn_sym(P23) + Fn_sym(P23)) + cF_C(3) * alpha(p[0]+p[1],p[2]) * (Fn_sym(P3)*_sptkernels.Gn_sym(P12) + Gn_sym(P12)) - cF_E(3) * beta(p[0],p[1]+p[2]) * (Gn_sym(P1)*_sptkernels.Gn_sym(P23) + Gn_sym(P23)) - cF_E(3) * beta(p[0]+p[1],p[2]) * (Gn_sym(P3)*_sptkernels.Gn_sym(P12) + Gn_sym(P12)) + cF_E(3) * (p[0]+p[1]+p[2])*_tau(p) + cF_C(3) * alphaOmega(p[0]+p[1],p[2]) * _omega(P12) + cF_E(3) * betaOmega(p[0]+p[1],p[2]) * _omega(P12);
    }
       
    return Fnval;
}
   
//------------------------------------------------------------------------------
double EFTkernels::Gn(const std::vector<ThreeVector>& p)
{
   int n = p.size();
   double Gnval=0;
      
   //Handle trivial cases
   if (n == 0) { return Gnval; }
   if (n > 3) { std::cout<<"There is no EFT kernel available at the order specified"<<std::endl; return Gnval; }
   
   if (n==1) { Gnval = cG_E(1) * (p[0]*_tau(p)) ;}
      
   if (n==2) {
         
      std::vector<ThreeVector> P1={p[0]};
      std::vector<ThreeVector> P2={p[1]};
         
      Gnval = cG_C(2) * alpha(p[0],p[1]) * (Gn_sym(P1) + Fn_sym(P2)) - cG_E(2) * beta(p[0],p[1]) * (Gn_sym(P1) + Gn_sym(P2)) + cG_E(2) * (p[0]+p[1])*_tau(p);
   }
      
   if (n==3) {
         
      std::vector<ThreeVector> P1={p[0]};
      std::vector<ThreeVector> P23={p[1],p[2]};
      std::vector<ThreeVector> P12={p[0],p[1]};
      std::vector<ThreeVector> P3={p[2]};
         
      Gnval = cG_C(3) * alpha(p[0],p[1]+p[2]) * (Gn_sym(P1)*_sptkernels.Fn_sym(P23) + Fn_sym(P23)) + cG_C(3) * alpha(p[0]+p[1],p[2]) * (Fn_sym(P3)*_sptkernels.Gn_sym(P12) + Gn_sym(P12)) - cG_E(3) * beta(p[0],p[1]+p[2]) * (Gn_sym(P1)*_sptkernels.Gn_sym(P23) + Gn_sym(P23)) - cG_E(3) * beta(p[0]+p[1],p[2]) * (Gn_sym(P3)*_sptkernels.Gn_sym(P12) + Gn_sym(P12)) + cG_E(3) * (p[0]+p[1]+p[2])*_tau(p) + cG_C(3) * alphaOmega(p[0]+p[1],p[2]) * _omega(P12) + cG_E(3) * betaOmega(p[0]+p[1],p[2]) * _omega(P12);
   }
      
   return Gnval;
}
   
//------------------------------------------------------------------------------
double EFTkernels::Fn_sym(const std::vector<ThreeVector>& p)
{
   double value = 0;
   int nperm = 0; // count the permutations
   // use the next_permutation algorithm together with the comparison operator
   // in ThreeVector to generate permutations
   std::vector<ThreeVector> pperm = p;
   std::sort(pperm.begin(),pperm.end());
   do {
      nperm++;
      value += Fn(pperm);
   } while (std::next_permutation(pperm.begin(), pperm.end()));

   return value / nperm;
}

//------------------------------------------------------------------------------
double EFTkernels::Gn_sym(const std::vector<ThreeVector>& p)
{
   double value = 0;
   int nperm = 0; // count the permutations
   // use the next_permutation algorithm together with the comparison operator
   // in ThreeVector to generate permutations
   std::vector<ThreeVector> pperm = p;
   std::sort(pperm.begin(),pperm.end());
   do {
      nperm++;
      value += Gn(pperm);
   } while (std::next_permutation(pperm.begin(), pperm.end()));

   return value / nperm;
}
   
//------------------------------------------------------------------------------
//Helper function
   
ThreeVector EFTkernels::_dot_product(std::vector<double> a, std::vector<ThreeVector> b)
{
   ThreeVector dot(0.,0.,0.);
   if(a.size() == b.size()){for(unsigned int i=0; i<a.size(); i++){dot+=a[i]*b[i];}}
   else { throw std::invalid_argument( "Received invalid argument in dot product. Size of vectors does not match" );}
      
   return dot;
}

} // namespace fnfast
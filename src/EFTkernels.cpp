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


//------------------------------------------------------------------------------
//EFT coefficients
//------------------------------------------------------------------------------


string EFTcoefficients::description()
{    
    stringstream stream;
    
    stream << "There are 14 coefficients in total. Call description(coefficient name) for details."<<endl
    <<"LO:   2 coefficients  -> cs,ch"<<endl
    <<"NLO:  9 coefficients  -> c1,c2,c3,c4,c5,c6,ch1,ch2,ch3"<<endl
    <<"NNLO: 3 coefficients  -> d1,d2,d3";
    
    return stream.str();
}

//------------------------------------------------------------------------------

string EFTcoefficients::description(EFTcoefficients::Labels label)
{
    string str = "Coefficient of ";
    
    switch(label){
            
        case EFTcoefficients::cs:
            str+= "k^2 x delta in Euler equation (speed of sound)";
            break;
        case EFTcoefficients::ch:
            str+= "k^2 x delta in continuity equation (heat capacity)";
            break;
        case EFTcoefficients::c1:
            str+= "k^2 x delta^2 in Euler equation (k=k1+k2)";
            break;
        case EFTcoefficients::c2:
            str+= "k^2(-1/3+(k1.k2)^2/(k1^2k2^2)) x delta^2 in Euler equation (k=k1+k2)";
            break;
        case EFTcoefficients::c3:
            str+= "-k^2/6+(k1.k2)/2(k.k1/k1^2) x delta^2 + perms in Euler equation (k=k1+k2)";
            break;
        case EFTcoefficients::c4:
            str+= "k^2 x delta x theta in Euler equation (k=k1+k2)";
            break;
        case EFTcoefficients::c5:
            str+= "k^2(-1/3+(k1.k2)^2/(k1^2k2^2)) x delta x theta in Euler equation (k=k1+k2)";
            break;
        case EFTcoefficients::c6:
            str+= "-k^2/6+(k1.k2)/2(k.k1/k1^2) x delta x theta + perms in Euler equation (k=k1+k2)";
            break;
        case EFTcoefficients::ch1:
            str+= "k^2 x delta^2 in continuity equation (k=k1+k2)";
            break;
        case EFTcoefficients::ch2:
            str+= "k^2(-1/3+(k1.k2)^2/(k1^2k2^2)) x delta^2 in continuity equation (k=k1+k2)";
            break;
        case EFTcoefficients::ch3:
            str+= "-k^2/6+(k1.k2)/2(k.k1/k1^2) x delta^2 + perms in continuity equation (k=k1+k2)";
            break;
        case EFTcoefficients::d1:
            str+= "(k.k1)^2(k2.k3)^2/(k1^2k2^2k3^2) x delta^3 + perms in Euler equation (k=k1+k2+k3)";
            break;
        case EFTcoefficients::d2:
            str+= "(k.k1)(k.k2)(k1.k3)(k2.k3)/(k1^2k2^2k3^2) x delta^3 + perms in Euler equation (k=k1+k2+k3)";
            break;
        case EFTcoefficients::d3:
            str+= "(k.k1)^2/k1^2 x delta^3 + perms in Euler equation (k=k1+k2+k3)";
            break;
    }
    return str;
}

//------------------------------------------------------------------------------
//EFT kernels
//------------------------------------------------------------------------------


double EFTkernels::cF_1(int n)
{
    return (2*n + 5.) / (2 * n * n + 9 * n + 7);
}

//------------------------------------------------------------------------------
double EFTkernels::cF_2(int n)
{
    return (-2.) / (2 * n * n + 9 * n + 7);

}

//------------------------------------------------------------------------------
double EFTkernels::cG_1(int n)
{
    return (3.) / (2 * n * n + 9 * n + 7);
}

//------------------------------------------------------------------------------
double EFTkernels::cG_2(int n)
{
    return (-2*n - 4.) / (2 * n * n + 9 * n + 7);
}

//------------------------------------------------------------------------------

//LO shapes
vector<double> EFTkernels::_lo_shapes(ThreeVector p)
{
    vector<double> shapes;
    shapes.push_back(p*p);
    
    return shapes;
}

//NLO shapes
vector<double> EFTkernels::_nlo_shapes(ThreeVector p1, ThreeVector p2)
{
    double eps = 1e-12;
    double shape1,shape2,shape3;
    ThreeVector p = p1+p2;
    
    //Shape1 doesn't need IR regulation
    shape1=p*p;
    
    //Shape2
    if(p1*p1 < eps || p2*p2 < eps) shape2 = 0;
    else shape2=(p*p)*(-1./3+(p1*p2)*(p1*p2)/((p1*p1)*(p2*p2)));
    
    //Shape3
    if(p1*p1 < eps || p2*p2 < eps) shape3 = 0;
    else shape3=-(p*p)/6+(p1*p2)/2*(p*p1/(p1*p1)+p*p2/(p2*p2));
    
    vector<double> shapes={shape1,shape2,shape3};
    
    return shapes;
}

//NNLO shapes
vector<double> EFTkernels::_nnlo_shapes(ThreeVector p1, ThreeVector p2, ThreeVector p3)
{
    double eps = 1e-12;
    double shape1,shape2,shape3;
    ThreeVector p = p1+p2+p3;
    
    //Shape1
    if(p1*p1 < eps || p2*p2 < eps || p3*p3 <eps) shape1=0;
    else shape1=((p*p1)*(p*p1)*(p2*p3)*(p2*p3)+(p*p2)*(p*p2)*(p1*p3)*(p1*p3)+(p*p3)*(p*p3)*(p1*p2)*(p1*p2))/((p1*p1)*(p2*p2)*(p3*p3));
    
    //Shape2
    if(p1*p1 < eps || p2*p2 < eps || p3*p3 <eps) shape2=0;
    else shape2=((p*p1)*(p*p2)*(p1*p3)*(p2*p3)+(p*p1)*(p*p3)*(p1*p2)*(p2*p3)+(p*p2)*(p*p3)*(p1*p2)*(p1*p3))/((p1*p1)*(p2*p2)*(p3*p3));
    
    //Shape3
    if(p1*p1 < eps || p2*p2 < eps || p3*p3 <eps) shape3=0;
    else shape3=(p*p1)*(p*p1)/(p1*p1)+(p*p2)*(p*p2)/(p2*p2)+(p*p3)*(p*p3)/(p3*p3);
    
    vector<double> shapes = {shape1,shape2,shape3};
    
    return shapes;
}

//------------------------------------------------------------------------------

double EFTkernels::_dot_product(vector<double> a, vector<double> b)
{
    double dot=0;
    if(a.size() == b.size()){for(unsigned int i=0; i<a.size(); i++){dot+=a[i]*b[i];}}
    else { throw std::invalid_argument( "Received invalid argument in dot product. Size of vectors does not match" );}
    
    return dot;
}

//------------------------------------------------------------------------------
double EFTkernels::Fn(vector<ThreeVector> p)
{
    int n = p.size();
   
    //Handle trivial cases
    if (n == 0) { return 0; }
    if (n > 3) { cout<<"There is no EFT kernel available at the order specified"<<endl; return 0; }
    
    double Fnval=0;
    
    if (n==1) { Fnval= cF_2(1)*(*_coefficients)[EFTcoefficients::cs]*_lo_shapes(p[0])[0]+cF_1(1)*(*_coefficients)[EFTcoefficients::ch]*_lo_shapes(p[0])[0];}
   
    if (n==2) {
        
        vector<double> c_vals_2={(*_coefficients)[EFTcoefficients::c1]-(*_coefficients)[EFTcoefficients::c4],(*_coefficients)[EFTcoefficients::c2]-(*_coefficients)[EFTcoefficients::c5],(*_coefficients)[EFTcoefficients::c3]-(*_coefficients)[EFTcoefficients::c6]};
        vector<double> ch_vals_2={(*_coefficients)[EFTcoefficients::ch1],(*_coefficients)[EFTcoefficients::ch2],(*_coefficients)[EFTcoefficients::ch3]};
        
        vector<ThreeVector> p0={p[0]};
        vector<ThreeVector> p1={p[1]};
        vector<ThreeVector> p01={p[0],p[1]};

        
        Fnval= cF_1(2)*alpha(p[0],p[1])*(Gn(p0)+Fn(p1))-cF_2(2)*beta(p[0],p[1])*(Gn(p0)+Gn(p1))
        +cF_2(2)*(*_coefficients)[EFTcoefficients::cs]*_lo_shapes(p[0]+p[1])[0]*_sptkernels.Fn(p01)
        +cF_1(2)*(*_coefficients)[EFTcoefficients::ch]*_lo_shapes(p[0]+p[1])[0]*_sptkernels.Fn(p01)
        +cF_2(2)*_dot_product(c_vals_2,_nlo_shapes(p[0],p[1]))
        +cF_1(2)*_dot_product(ch_vals_2,_nlo_shapes(p[0],p[1]));
    }
    
    
    if (n==3) {
        
        vector<double> cdelta_vals_3={(*_coefficients)[EFTcoefficients::c1],(*_coefficients)[EFTcoefficients::c2],(*_coefficients)[EFTcoefficients::c3]};
        vector<double> ctheta_vals_3={(*_coefficients)[EFTcoefficients::c4],(*_coefficients)[EFTcoefficients::c5],(*_coefficients)[EFTcoefficients::c6]};
        vector<double> ch_vals_3={(*_coefficients)[EFTcoefficients::ch1],(*_coefficients)[EFTcoefficients::ch2],(*_coefficients)[EFTcoefficients::ch3]};
        vector<double> d_vals_3={(*_coefficients)[EFTcoefficients::d1],(*_coefficients)[EFTcoefficients::d2],(*_coefficients)[EFTcoefficients::d3]};
        
        vector<ThreeVector> p0={p[0]};
        vector<ThreeVector> p1={p[1]};
        vector<ThreeVector> p2={p[2]};
        vector<ThreeVector> p01={p[0],p[1]};
        vector<ThreeVector> p12={p[1],p[2]};
        vector<ThreeVector> p012={p[0],p[1],p[2]};
        
        
        Fnval= cF_1(3)*alpha(p[0],p[1]+p[2])*(Gn(p0)*_sptkernels.Fn(p12)+Fn(p12))
        +cF_1(3)*alpha(p[0]+p[1],p[2])*(Gn(p01)+_sptkernels.Gn(p01)*Fn(p2))
        -cF_2(3)*beta(p[0],p[1]+p[2])*(Gn(p0)*_sptkernels.Gn(p12)+Gn(p12))
        -cF_2(3)*beta(p[0]+p[1],p[2])*(Gn(p01)+_sptkernels.Gn(p01)*Gn(p2))
        +cF_2(3)*(*_coefficients)[EFTcoefficients::cs]*_lo_shapes(p[0]+p[1]+p[2])[0]*_sptkernels.Fn(p012)
        +cF_1(3)*(*_coefficients)[EFTcoefficients::ch]*_lo_shapes(p[0]+p[1]+p[2])[0]*_sptkernels.Fn(p012)
        +cF_2(3)*_dot_product(cdelta_vals_3,_nlo_shapes(p[0],p[1]+p[2]))*_sptkernels.Fn(p12)
        -cF_2(3)*_dot_product(ctheta_vals_3,_nlo_shapes(p[0],p[1]+p[2]))*_sptkernels.Gn(p12)
        +cF_2(3)*_dot_product(cdelta_vals_3,_nlo_shapes(p[0]+p[1],p[2]))*_sptkernels.Fn(p01)
        -cF_2(3)*_dot_product(ctheta_vals_3,_nlo_shapes(p[0]+p[1],p[2]))*_sptkernels.Fn(p01)
        +cF_1(3)*_dot_product(ch_vals_3,_nlo_shapes(p[0],p[1]+p[2]))*_sptkernels.Fn(p12)
        +cF_1(3)*_dot_product(ch_vals_3,_nlo_shapes(p[0]+p[1],p[2]))*_sptkernels.Fn(p01)
        +cF_2(3)*_dot_product(d_vals_3,_nnlo_shapes(p[0],p[1],p[2]));
    }

    return Fnval;
}

//------------------------------------------------------------------------------
double EFTkernels::Gn(vector<ThreeVector> p)
{
    int n = p.size();
    
    //Handle trivial cases
    if (n == 0) { return 0; }
    if (n > 2) { cout<<"There is no EFT kernel available at the order specified"<<endl; return 0; }
    
    double Gnval=0;
    
    if (n==1) { Gnval= cG_2(1)*(*_coefficients)[EFTcoefficients::cs]*_lo_shapes(p[0])[0]+cG_1(1)*(*_coefficients)[EFTcoefficients::ch]*_lo_shapes(p[0])[0];}
    
    if (n==2) {
        
        vector<double> c_vals_2={(*_coefficients)[EFTcoefficients::c1]-(*_coefficients)[EFTcoefficients::c4],(*_coefficients)[EFTcoefficients::c2]-(*_coefficients)[EFTcoefficients::c5],(*_coefficients)[EFTcoefficients::c3]-(*_coefficients)[EFTcoefficients::c6]};
        vector<double> ch_vals_2={(*_coefficients)[EFTcoefficients::ch1],(*_coefficients)[EFTcoefficients::ch2],(*_coefficients)[EFTcoefficients::ch3]};
        
        vector<ThreeVector> p0={p[0]};
        vector<ThreeVector> p1={p[1]};
        vector<ThreeVector> p01={p[0],p[1]};
        
        Gnval= cG_1(2)*alpha(p[0],p[1])*(Gn(p0)+Fn(p1))-cG_2(2)*beta(p[0],p[1])*(Gn(p0)+Gn(p1))
        +cG_2(2)*(*_coefficients)[EFTcoefficients::cs]*_lo_shapes(p[0]+p[1])[0]*_sptkernels.Fn(p01)
        +cG_1(2)*(*_coefficients)[EFTcoefficients::ch]*_lo_shapes(p[0]+p[1])[0]*_sptkernels.Fn(p01)
        +cG_2(2)*_dot_product(c_vals_2,_nlo_shapes(p[0],p[1]))
        +cG_1(2)*_dot_product(ch_vals_2,_nlo_shapes(p[0],p[1]));
    }
    
   return Gnval;
}

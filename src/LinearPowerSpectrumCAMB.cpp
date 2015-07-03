//------------------------------------------------------------------------------
/// \file LinearPowerSpectrumCAMB.cpp
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
//    Implementation of objects in LinearPowerSpectrumCAMB
//------------------------------------------------------------------------------

#include "LinearPowerSpectrumCAMB.hpp"

using namespace std;

//------------------------------------------------------------------------------
LinearPowerSpectrumCAMB::LinearPowerSpectrumCAMB(const string& input_file): _input_file(input_file), _accel_ptr(NULL), _spline_ptr(NULL)
{
    ifstream file;
    file.open(_input_file);
    
    if(!file.good()) { cout << "LinearPowerSpectrumCAMB : I can't open the requested file, " << input_file << endl; }
    
    else {
        
        // Read file into vectors
        int npts = 0;
        while(!file.eof()) {
            double k,P;
            file >> k >> P;
            _kvec.push_back(k);
            _Pvec.push_back(P);
            npts++;
        }
        
        // Define patches at low and high k
        // Use first and last 10 points to get tails behavior
        double k_low[10], P_low[10], k_high[10], P_high[10];
       
        for(size_t i=0; i<10; i++)
        {
            k_low[i] = log(_kvec[i]);
            P_low[i] = log(_Pvec[i]);
            int j = _kvec.size() - 10 - i;
            k_high[i] = log(_kvec[j]);
            P_high[i] = log(_Pvec[j]);
        }
        
        // Fitting function returns linear fit coefficients and covariance and chisq parameters
        double cov1_low,cov2_low,cov3_low,chisq_low,cov1_high,cov2_high,cov3_high,chisq_high;
        
        gsl_fit_linear (k_low, 1, P_low, 1, 10, &_c0_low, &_c1_low, &cov1_low, &cov2_low, &cov3_low, &chisq_low);
        gsl_fit_linear (k_high, 1, P_high, 1, 10, &_c0_high, &_c1_high, &cov1_high, &cov2_high, &cov3_high, &chisq_high);
        
        // Generate a set of points (equally spaced in log)
        int npts_patch = 10;
        vector<double> k_low_patch = _log_gen(_kvec.front() / 10, _kvec.front(), npts_patch);
        vector<double> k_high_patch = _log_gen(_kvec.back(), _kvec.back() * 10, npts_patch);
        
        vector<double> P_low_patch,P_high_patch;
        
        // Define patches
        for(size_t i = 0; i < k_low_patch.size(); i++)
        {
            P_low_patch.push_back( exp(_c0_low) * pow(k_low_patch[i],_c1_low) );
            P_high_patch.push_back( exp(_c0_high) * pow(k_high_patch[i],_c1_high) );
        }
        
        // Merge patches
        _kvec_patches.reserve(k_low_patch.size() + _kvec.size() - 2 + k_high_patch.size());
        _Pvec_patches.reserve(P_low_patch.size() + _Pvec.size() - 2 + P_high_patch.size());
        
        _kvec_patches = k_low_patch;
        _kvec_patches.insert(_kvec_patches.end(), _kvec.begin() + 1, _kvec.end() - 1);
        _kvec_patches.insert(_kvec_patches.end(), k_high_patch.begin(), k_high_patch.end());
        
        _Pvec_patches = P_low_patch;
        _Pvec_patches.insert(_Pvec_patches.end(), _Pvec.begin() + 1, _Pvec.end() - 1);
        _Pvec_patches.insert(_Pvec_patches.end(), P_high_patch.begin(), P_high_patch.end());
        
        // Arrays for interpolation
        double *k_vals;
        double *P_vals;
        
        int npts_tot = _kvec_patches.size();
        
        k_vals = new double[npts_tot];
        P_vals = new double[npts_tot];
        
        // Fill in arrays
        copy(_kvec_patches.begin(),_kvec_patches.end(),k_vals);
        copy(_Pvec_patches.begin(),_Pvec_patches.end(),P_vals);
        
        // Allocate interpolation pointers and initialize interpolation
        _accel_ptr = gsl_interp_accel_alloc();
        _spline_ptr = gsl_spline_alloc (gsl_interp_cspline, npts_tot);
        gsl_spline_init(_spline_ptr, k_vals, P_vals, npts_tot);
    }
    
    file.close();
}

//------------------------------------------------------------------------------
double LinearPowerSpectrumCAMB::operator()(double x)
{
   double res = 0;

   if( x >= 0 && x < _kvec_patches.front()) res = exp(_c0_low) * pow(x,_c1_low);  // Patch at low k
   if( x >= _kvec_patches.front() && x < _kvec_patches.back()) res = gsl_spline_eval(_spline_ptr, x, _accel_ptr); // Interpolated function
   if( x >= _kvec_patches.back()) res = exp(_c0_high) * pow(x,_c1_high); // Patch at high k

   return res;
}

//------------------------------------------------------------------------------
vector<double> LinearPowerSpectrumCAMB::_log_gen(double xmin, double xmax, int n)
{
    double logxmin = log(xmin);
    double logxmax = log(xmax);
    
    vector<double> res;
    
    // Use equally spaced exponents
    for(size_t i = 0; i <= n+1; i++)
    {
        double exponent = logxmin + (logxmax - logxmin) / (n + 1) * i;
        res.push_back( exp(exponent) );
    }
    
    return res;
}

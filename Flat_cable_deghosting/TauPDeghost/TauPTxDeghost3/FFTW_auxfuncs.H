/*
 * =====================================================================================
 *
 *       Filename:  FFTW_auxfuncs.H
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/16/2015 02:14:04 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thinh P. Tran (), tpthinh@fairfield.com.vn, trphthinh@gmail.com
 *   Organization:  Fairfield Vietnam Ltd.
 *
 * =====================================================================================
 */


#ifndef FFTWauxfuncs_H
#define	FFTWauxfuncs_H

#include <iostream>
#include <cstdlib>
#include <complex>
#include <fftw3.h>

using namespace std;

void uniform_fftw1d(int length, fftw_complex * in, fftw_complex * out, int sign);
void uniform_forward_fftw1d_drc(int length, double* in, complex<double>* out);
void uniform_backward_fftw1d_dcr(int length, complex<double> * in, double * out);


#endif	/* FFTWauxfuncs_H */

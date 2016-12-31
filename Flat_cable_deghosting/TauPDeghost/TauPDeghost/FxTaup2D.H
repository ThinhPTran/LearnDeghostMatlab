/* 
 * File:   FxTaup2D.H
 * Author: dragon
 *
 * Created on December 12, 2013, 11:10 AM
 */

#ifndef FXTAUP2D_H
#define	FXTAUP2D_H


#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include "Physics_constant.H"
#include "fftw3.h"
#include "Taup2D_info.H"



using namespace std;


class FxTaup2D {
public:
  FxTaup2D();
  ~FxTaup2D();
  

// Public method
public:
  void initialize(TauP2D_info_t info);
  void execute();
  void finalize();
  
  
// Private method
private:
  void forward_fx_taup_transform();
  void inverse_fx_taup_transform();
  void initialize_x();
  void initialize_p();
  void uniform_fftw1d (
                        int length,
                        fftw_complex * in,
                        fftw_complex * out,
                        int sign
                       );
  void transf_from_t_to_f(complex<double> * output);
  void transf_from_tau_to_f(complex<double> * output);
  void transf_from_f_to_tau(complex<double> * input);
  void transf_from_f_to_t(complex<double> * input);
  void recalculating_A(int freq_idx, complex<double> * A);
  void multiply_A_by_freq_slice(complex<double>* in, complex<double>* out, complex<double>* A);
  void linear_sys_solver ( complex<double>* in, complex<double>* out, complex<double>* A);
  void forward_operator(complex<double> * input, complex<double> * output, complex<double>* A);
  void adjoint_operator(complex<double> * input, complex<double> * output, complex<double>* A);
  double norm_vec(complex<double> * vec, int vec_length);
  
  
// Private variable
private:
  int nt;
  int ntau;
  int nx;
  int np;
  int nf;
  
  double dt;
  double dtau;
  double dx;
  double dp;
  double df;
  
  double tmin;
  double xmin;
  double pmin;  
  
  double T;
  double Tau;
  
  taup_mode_t mode;
  int iter_cg;
  
  double * x;
  double * p;
  
  double * in;
  double * out;
  

};

#endif	/* FXTAUP2D_H */


/* 
 * File:   TxTaup2D.H
 * Author: dragon
 *
 * Created on December 12, 2013, 11:10 AM
 */

#ifndef TXTAUP2D_H
#define	TXTAUP2D_H


#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "Physics_constant.H"
#include "fftw3.h"
#include "Taup2D_info.H"



using namespace std;


class TxTaup2D {
public:
  TxTaup2D();
  ~TxTaup2D();
  

// Public method
public:
  void initialize(TauP2D_info_t info);
  void execute();
  void finalize();
  
  
// Private method
private:
  void forward_tx_taup_transform();
  void inverse_tx_taup_transform();
  void inverse_slant_stack();
  void calculate_dependent_parms();
  void initialize_x();
  void initialize_p();
  void initialize_omega();
  void initialize_filter();
  void uniform_fftw1d (
                        int length,
                        fftw_complex * in,
                        fftw_complex * out,
                        int sign
                      );
  void uniform_forward_fftw1d_drc 
                       (
                        int length,
                        double * in,
                        complex<double> * out
                       );
  void uniform_backward_fftw1d_dcr 
                       (
                        int length,
                        complex<double> * in,
                        double * out
                       );
  void transf_from_t_to_f(double * tracein, complex<double> * traceout, int nx, int nt);
  void transf_from_f_to_t(complex<double> * tracein, double * traceout, int nx, int nt);
  void apply_rhofilter();
  void filter_one_trace(complex<double> * trace);
  double norm_vec(complex<double> * vec, int vec_length);
  void set_max_value();
  double absmax(double * input, int lgth);
  void rescale_output();
  
  
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
  double domega;
  
  double tmin;
  double xmin;
  double pmin; 
  double fmax; 
  int iwmax; 
  
  double * maxval;
  
  double T;
  double Tau;
  
  taup_mode_t mode;
  int iter_cg;
  
  double * omega;
  double * x;
  double * p;
  double * filter;
  
  double * in;
  double * out;
  

};

#endif	/* TXTAUP2D_H */


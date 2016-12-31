/* 
 * File:   TauPDeghost.H
 * Author: dragon
 *
 * Created on May 7, 2014, 8:37 AM
 */

#ifndef TAUP2DDEGHOST_H
#define	TAUP2DDEGHOST_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "Physics_constant.H"
#include "Taup2D_info.H"
#include "fftw3.h"
#include "TxTauP2D.H"
#include "TauP2DDeghost_info.H"
#include "FFTW_auxfuncs.H"


using namespace std;


class TauP2DDeghost {
  
public:
  TauP2DDeghost();
  ~TauP2DDeghost();
  
  
// Public method
public:
  void initialize(TauP2DDeghost_info_t info);
  void execute();
  void finalize();
  
  
// Private method
private:
  double maxabs(complex<double> * vec, int length);
  double maxabs(double * vec, int length);
  void get_input_info(TauP2DDeghost_info_t info);
  void check_input_info();
  void calculate_dependent_parameters();
  void display_info_for_checking();
  void allocate_main_mem();
  void precompute_variables();
  void compute_p();
  void compute_omega();
  void compute_exkx();
  void compute_exomega();
  void compute_phys_filter();
  void get_exdata();
  void get_output_data();
  void biwrite(string file_name, double * data, int M, int nt);
  void biwrite_complex(string file_name, complex<double> * data, int M, int nt);
  void phys_filtering();
  void uniform_fftwf2d(int nt, int nx, complex<double> * in, complex<double> * out);
  void uniform_ifftwf2d(int nt, int nx, complex<double> * in, complex<double> * out);
  void fwd_taup_transform(int nt, int nx, int np, double dt, double dx, double dp, 
                          double tmin, double xmin, double pmin, double * in, double * out);
  void bwd_taup_transform(int nt, int nx, int np, double dt, double dx, double dp, 
                          double tmin, double xmin, double pmin, double * in, double * out);
  void deghost();
  void cleanup(double * taup_img_in, double * taup_img_out);
  void fanalyse_spec_fwd(double* input, complex<double> * fspecout);
  void extanalyse_spec_fwd(double * input, double * tspecout);
  void calculate_freq_hz(double * freq_hz);
  void fourier_img(double * input, complex<double> * finput, int nt);
  void inv_fourier(complex<double> * finput, double * output, int nt);
  void calculate_wavelet(double * freq_hz, double * wavelet);
  void excalculate_wavelet(double * freq_hz, double * wavelet);
  void taup_deghost(double * taup_img_in, double * taup_img_out);
  void taup_deghost(complex<double> * ftaup_img, double depth);
  void taup_deghost_tau_to_f(double * taup_img, complex<double> * ftaup_img);
  void taup_deghost_f_to_tau(complex<double> * ftaup_img, double * taup_img);
  void calculate_exmissing_part(double * exinput, double * exinv);
  
  void scale_back_for_fix_display();
  void add_back_missing_part();
  void reshape_spectrum();
  
  
  
  
private:
  // Input info
  double vwater;
  double recover_eps;
  double src_depth;
  double rec_depth;
  
  int nt;
  int nx;
  
  int npad;
  
  double dt;
  double df;
  double domega;
  double dx;
  double dp;
  
  double fx;
  
  double * in;
  double * out;
  
  
  // Calculated info
  int np;
  int nkx;
  
  int exnt;
  int exnx;
  int exnkx;
  
  
  double dkx;
  
  double exdt;
  double exdf;
  double exdomega;
  double exdx;
  double exdkx;
  
  
  double ft;
  double fp;
  double fomega;
  
  double exfx;
  double exft;
  double exfomega;
  
  
  double pmax;
  
  int nfb;
  double fmin;
  double fmax;
  double res;
  
  double mfapp;
  
  int ifscale;
  int ifaddback;
  int ifreshape;
  
  double p_control;
  double taup_maxval;
  
  
  // Precomputing variables
  double * p;
  double * f;
  double * omega;
  double * exkx;
  double * exf;
  double * exomega;
  double * phys_filter;
  
  
  // exdata
  complex<double> * exdata;
  complex<double> * exmissing_part;

};



#endif	/* TAUP2DDEGHOST_H */


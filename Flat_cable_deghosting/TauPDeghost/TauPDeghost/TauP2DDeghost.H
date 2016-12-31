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
#include "FxTaup2D.H"
#include "TauP2DDeghost_info.H"
#include "fftw3.h"


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
  void get_input_info(TauP2DDeghost_info_t info);
  void calculate_dependent_parameters();
  void display_info_for_checking();
  void allocate_main_mem();
  void precompute_variables();
  void compute_p();
  void compute_kx();
  void compute_exomega();
  void compute_phys_filter();
  void get_exdata();
  void output_data();
  void biwrite(string file_name, double * data, int M, int nt);
  void biwrite_complex(string file_name, complex<double> * data, int M, int nt);
  void phys_filtering();
  void uniform_fftwf1d(int Nt, fftw_complex * in, fftw_complex * out, int sign);
  void uniform_fftwf2d(int nt, int nx, complex<double> * in, complex<double> * out);
  void uniform_ifftwf2d(int nt, int nx, complex<double> * in, complex<double> * out);
  void fwd_taup_transform(int nt, int nx, int np, double dt, double dx, double dp, 
                          double tmin, double xmin, double pmin, double * in, double * out);
  void bwd_taup_transform(int nt, int nx, int np, double dt, double dx, double dp, 
                          double tmin, double xmin, double pmin, double * in, double * out);
  void deghost();
  void taup_deghost(double * taup_img);
  void taup_deghost(complex<double> * ftaup_img, double depth);
  void taup_deghost_tau_to_f(double * taup_img, complex<double> * ftaup_img);
  void taup_deghost_f_to_tau(complex<double> * ftaup_img, double * taup_img);
  
  
  
private:
  // Input info
  double vwater;
  double recover_eps;
  double src_depth;
  double rec_depth;
  
  int nt;
  int nx;
  
  double dt;
  double dx;
  double dp;
  
  double fx;
  
  double * in;
  double * out;
  
  
  // Calculated info
  int np;
  int nkx;
  int exnt;
  
  double dkx;
  double exdt;
  double exdf;
  double exdomega;
  
  double ft;
  double fp;
  double fomega;
  
  double pmax;
  
  
  // Precomputing variables
  double * p;
  double * kx;
  double * exomega;
  double * phys_filter;
  
  
  // exdata
  complex<double> * exdata;
  

};



#endif	/* TAUP2DDEGHOST_H */


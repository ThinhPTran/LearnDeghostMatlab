/* 
 * File:   General_deghost_engine.H
 * Author: thinh
 *
 * Created on September 13, 2013, 11:52 AM
 */

#ifndef GENERAL_DEGHOST_ENGINE_H
#define	GENERAL_DEGHOST_ENGINE_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <fftw3.h>
#include "ghost_data.H"
#include <omp.h>


using namespace std;


class General_deghost_engine {
  
public:
  General_deghost_engine();
  ~General_deghost_engine();
  void process(ghost_data * ghost_data);
  
  
private:
  void initialize();
  void execute(ghost_data * ghostdata);
  void finalize();
  void print_vec_complex(complex<float> * vec, int length);
  template<class T> 
  void print_vec(T * vec, int length);
  void print_vec_float(float * vec, int length);
  void calculate_basic_info();
  void display_info();
  float fmaxf(int length, float * vec);
  float fmaxabsf(int length, float * vec);
  float fmaxabscomplex(int length, complex<float> * vec);
  int ind_max(int length, float * vec);
  template<class T>
  T mean(int length, T * vec);
  template<class T>
  T meanabs(int length, T * vec);
  void uniform_fftw1d(int length, fftwf_complex * in, fftwf_complex * out, int sign);
  void initialize_trace_in_trace_out();
  void calculate_F2_and_denominator();
  void initialize_Fourier_image();
  void initialize_check_array();
  void calculate_Fourier_image();
  void calculate_check_array();
  void calculate_tau();
  void calculate_raw_output();
  void calculate_output(ghost_data * ghostdata);
  
  
private:
  static const complex<float> I;
  static const float Pi;
  static const float nscriteriastandard;
  
  
  // Variables that must be initialized
  int first_run;
  int recompute; // Recompute F2 and denominator
  int nt;
  float dt;
  float zmax;
  float zmin;
  float vwater;
  float epsilon;
  float addback;
  
  
  // Variables that can be calculated 
  float df;
  float tmax;
  float tmin;
  float fmax;
  float fmin;
  
  int ntcheck;
  float tmincheck;
  float tmaxcheck;
  float avgcheck;
  float input_average;
  float tau;
  float nscriteria;
  
  // Variables may be used for many executions
  complex<float> * F2;
  complex<float> * denominator;
  
  
  // Variables are used for only one execution locally
  complex<float> * Fourier_image;
  float * check;
  
  
  // Variables are used for only one execution given by user.
  float * trace_in;
  float * trace_out;  

};


#endif	/* GENERAL_DEGHOST_ENGINE_H */


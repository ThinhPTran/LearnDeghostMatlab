/* 
 * File:   General_deghost_engine.C
 * Author: thinh
 * 
 * Created on September 13, 2013, 11:52 AM
 */

#include "General_deghost_engine.H"
#include "ghost_data.H"


const complex<float> General_deghost_engine::I = complex<float>(0,1);
const float General_deghost_engine::Pi =  3.14159265358979323846;
const float General_deghost_engine::nscriteriastandard = 1.3;


General_deghost_engine::General_deghost_engine() {
  
  // 
  first_run = 1;
  recompute = 0;  int nt;
  dt = 0.0;
  zmax = 0.0;
  zmin = 0.0;
  vwater = 0.0;
  epsilon = 0.00001;
  
  
  // Variables that can be calculated 
  df = 0.0;
  tmax = 0.0;
  tmin = 0.0;
  fmax = 0.0;
  fmin = 0.0;
  
  ntcheck = 0;
  tmincheck = 0.0;
  tmaxcheck = 0.0;
  avgcheck = 0.0;
  input_average = 0;
  tau = 0;
  nscriteria = 0;
  
  
  // Used for many segment with the same length
  F2 = NULL;
  denominator = NULL;
  
  
  // Delete after using for one segment
  Fourier_image = NULL;
  check = NULL;
  
  
  // Given by the user
  trace_in = NULL;
  trace_out = NULL;  

  
}


General_deghost_engine::~General_deghost_engine() {
  
  if (F2 != NULL) {
    delete [] F2;
    F2 = NULL;
  }
  
  if (denominator != NULL) {
    delete [] denominator;
    denominator = NULL;
  }
  
    
  if ( trace_in != NULL ){
    delete [] trace_in;
    trace_in = NULL;
  }
  
  
  if ( trace_out != NULL ) {
    delete [] trace_out;
    trace_out = NULL;
  }
  
}


void General_deghost_engine::print_vec_complex
(
  complex<float> * vec, 
  int length
) 
{  
  
  cout << "Vector: ";
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (!(i_iter%4)) { 
      cout << endl;
    }
    cout << setw(10) << vec[i_iter] << "    ";
  }
  cout << endl;
  
}


void General_deghost_engine::print_vec_float
(
  float* vec, 
  int length
) 
{
  
  cout << "Vector: ";
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (!(i_iter%4)) { 
      cout << endl;
    }
    cout << setw(10) << vec[i_iter] << "    ";
  }
  cout << endl;
  
}


template<class T> 
void General_deghost_engine::print_vec
(
  T* vec, 
  int length
) 
{
  
  cout << "Vector: ";
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (!(i_iter%4)) { 
      cout << endl;
    }
    cout << setw(10) << vec[i_iter] << "    ";
  }
  cout << endl;
  
}


void General_deghost_engine::calculate_basic_info() {
  
  df = 1.0/((nt)*dt);
  tmax = (nt-1)*dt;
  tmin = 0;
  fmax = (nt-1)*df;
  fmin = 0;
  
  tmaxcheck = 1.5*zmax*2.0/vwater;
  tmincheck = zmin*2.0/vwater;
  ntcheck = floor(tmaxcheck/dt) + 1;
  
}


void General_deghost_engine::display_info() {
  
  cout << "********************************************************************" << endl;
  cout << "          General_deghost_engine::display_info()                    " << endl;
  cout << "********************************************************************" << endl;
  cout << "Given variables " << endl;
  cout << "first_run: " << first_run << endl;
  cout << "recompute: " << recompute << endl; // Recompute F2 and denominator
  cout << "nt: " << nt << endl;
  cout << "dt: " << dt << endl;
  cout << "zmin: " << zmin << endl;
  cout << "zmax: " << zmax << endl;
  cout << "vwater: " << vwater << endl;
  cout << "epsilon: " << epsilon << endl;
  cout << "addback: " << addback << endl;
   
  cout << endl;
  cout << "Calculated variables " << endl;
  cout << "df: " << df << endl;
  cout << "tmin: " << tmin << endl;
  cout << "tmax: " << tmax << endl;
  cout << "fmin: " << fmin << endl;
  cout << "fmax: " << fmax << endl;
  cout << "ntcheck: " << ntcheck << endl;
  cout << "tmincheck: " << tmincheck << endl;
  cout << "tmaxcheck: " << tmaxcheck << endl;
  cout << "average: " << input_average << endl;
  cout << "********************************************************************" << endl;
  
}


template<class T>
T General_deghost_engine::mean(int length, T* vec) {
  
  T ret = 0;
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    ret += vec[i_iter]; 
  }
  
  ret /= length;
  
  return ret;
  
}


template<class T>
T General_deghost_engine::meanabs(int length, T* vec) {
  
  T ret = 0; 
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    ret += abs(vec[i_iter]);
  }
  
  ret /= length;
  
  return ret;
  
}


float General_deghost_engine::fmaxf
(
  int length, 
  float* vec
) 
{
  
  float ret = vec[0];
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (vec[i_iter] > ret) {
      ret = vec[i_iter];
    }
  }
  
  return ret;
  
}


float General_deghost_engine::fmaxabsf
(
  int length, 
  float* vec
) 
{
  
  float ret = abs(vec[0]);
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (abs(vec[i_iter]) > ret) {
      ret = abs(vec[i_iter]);
    }
  }
  
  return ret;
  
}


float General_deghost_engine::fmaxabscomplex
(
  int length,
  complex<float>* vec
) 
{
  
  float ret = abs(vec[0]);
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (abs(vec[i_iter]) > ret) {
      ret = abs(vec[i_iter]);
    }
  }
  
  return ret;
  
}


int General_deghost_engine::ind_max
(
  int length, 
  float* vec
) 
{
  
  int ind = 0;
  float max = vec[0];
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    if (vec[i_iter] > max) {
      max = vec[i_iter];
      ind = i_iter;
    }
  }
  
  return ind;
  
}


//FFT function
void General_deghost_engine::uniform_fftw1d
(
int length,
fftwf_complex * in,
fftwf_complex * out,
int sign
)
{

  fftwf_plan p;

  // Preprocessing check
  if ((sign != FFTW_BACKWARD) && (sign != FFTW_FORWARD)) {
    cout << "Wrong option of sign!!!" << endl;
    exit(-1);
  }

#pragma omp critical
  p = fftwf_plan_dft_1d(length, in , out, sign, FFTW_ESTIMATE);
  fftwf_execute(p);

  // Post-processing check
  if (sign == FFTW_BACKWARD) {
      for (int i_iter = 0; i_iter < length; i_iter++) {
        out[i_iter][0] = out[i_iter][0]/(1.0*length);
        out[i_iter][1] = out[i_iter][1]/(1.0*length);
      }
  }

#pragma omp critical
  fftwf_destroy_plan(p);

}


void General_deghost_engine::initialize_Fourier_image() {
  
  if (Fourier_image!=NULL) {
    delete [] Fourier_image;
    Fourier_image = NULL;
  }
  
  Fourier_image = new complex<float>[nt];
  memset(Fourier_image,0x00,nt*sizeof(complex<float>));
  
}


void General_deghost_engine::initialize_check_array() {
  
  if (check!=NULL) {
    delete [] check;
    check = NULL;
  }
  
  check = new float[ntcheck];
  memset(check,0x00,ntcheck*sizeof(float));
  
}


void General_deghost_engine::calculate_Fourier_image() {
  
  int i_iter;
  fftwf_complex * tmp_in = NULL;
  fftwf_complex * tmp_out = NULL;
  
  tmp_in = (fftwf_complex *) fftwf_malloc(nt*sizeof(fftwf_complex));
  tmp_out = (fftwf_complex *) fftwf_malloc(nt*sizeof(fftwf_complex));
  
  
  // Initialize 
  for (i_iter = 0; i_iter < nt; i_iter++) {
    tmp_in[i_iter][0] = trace_in[i_iter];
    tmp_in[i_iter][1] = 0.0;
  }
  
  
  uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_FORWARD);
  
  
  // Get output
  for (i_iter = 0; i_iter < nt; i_iter++) {
    Fourier_image[i_iter] = complex<float>(tmp_out[i_iter][0],tmp_out[i_iter][1]);
  }
  
  
  fftwf_free(tmp_in);
  fftwf_free(tmp_out);
  
}


void General_deghost_engine::calculate_check_array() {
  
  int i_iter,j_iter;
  int findxcheck = floor(tmincheck/dt);
  float inversedmean;
  float maxvalue;
  complex<float> * fP = new complex<float>[nt];
  float * P = new float[nt];
  float * G = new float[nt];
  float * test = new float[nt];
  fftwf_complex * tmp_in = (fftwf_complex *) fftwf_malloc(nt*sizeof(fftwf_complex));
  fftwf_complex * tmp_out = (fftwf_complex *) fftwf_malloc(nt*sizeof(fftwf_complex));
  
  
  for (i_iter = 1; i_iter < ntcheck; i_iter++) {
    
    // Calculate fP
    for (j_iter = 0; j_iter < nt; j_iter++) {
      fP[j_iter] = Fourier_image[j_iter]*F2[i_iter*nt+j_iter]/(denominator[i_iter*nt+j_iter] + complex<float>(2.0,0.0));
    }
    
    // Keep it symmetry 
    for (j_iter = nt-1; j_iter > nt/2; j_iter--) {
      fP[j_iter] = conj(fP[nt - j_iter]);
    }
    
    for (j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = fP[j_iter].real();
      tmp_in[j_iter][1] = fP[j_iter].imag();
    }
    uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_BACKWARD);
    for (j_iter = 0; j_iter < nt; j_iter++) {
      P[j_iter] = tmp_out[j_iter][0];
    }
    
    // Calculate G
    for (j_iter = 0; j_iter < nt; j_iter++) {
      G[j_iter] = trace_in[j_iter] - P[j_iter];
    }
    
    for (j_iter = 0; j_iter < nt; j_iter++) {
      test[j_iter] = (G[j_iter] + P[(j_iter+i_iter)%nt])/input_average; 
    }
    
    
    inversedmean = 1.0/meanabs(nt,test);
    maxvalue = fmaxf(nt, P);
    maxvalue /= input_average;
    check[i_iter] = inversedmean + maxvalue;
    
  }
  
  
  fftwf_free(tmp_in);
  fftwf_free(tmp_out);
  delete [] test;
  delete [] G;
  delete [] P;
  delete [] fP;
  
}


void General_deghost_engine::calculate_raw_output() {
  
  int j_iter;
  complex<float> * fP = new complex<float>[nt];
  float * P = new float[nt];
  float * G = new float[nt];
  float * normweight = new float[nt];
  float maxweight;
  float max1;
  float max2;
  float scale;
  fftwf_complex * tmp_in = (fftwf_complex *) fftwf_malloc(nt*sizeof(fftwf_complex));
  fftwf_complex * tmp_out = (fftwf_complex *) fftwf_malloc(nt*sizeof(fftwf_complex));
  int ind = ind_max(ntcheck, check);
  
    
  // Calculate fP
  for (j_iter = 0; j_iter < nt; j_iter++) {
    fP[j_iter] = Fourier_image[j_iter]*F2[ind*nt+j_iter]/(denominator[ind*nt+j_iter] + complex<float>(epsilon,0.0));
  }

  // Keep it symmetry 
  for (j_iter = nt-1; j_iter > nt/2; j_iter--) {
    fP[j_iter] = conj(fP[nt - j_iter]);
  }

  for (j_iter = 0; j_iter < nt; j_iter++) {
    tmp_in[j_iter][0] = fP[j_iter].real();
    tmp_in[j_iter][1] = fP[j_iter].imag();
  }
  uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_BACKWARD);
  for (j_iter = 0; j_iter < nt; j_iter++) {
    P[j_iter] = tmp_out[j_iter][0];
  }
  
  
  
  for (j_iter = 0; j_iter < nt; j_iter++) {
    trace_out[j_iter] = (P[j_iter] + addback*trace_in[j_iter])/(1.0+addback);
  }
  
  
  fftwf_free(tmp_in);
  fftwf_free(tmp_out);
  delete [] P;
  delete [] G;
  delete [] normweight;
  delete [] fP;
  
}


void General_deghost_engine::calculate_output() {

  
  if ( (tau > tmaxcheck)||(tau < tmincheck) ) {
    for (int i_iter = 0; i_iter < nt; i_iter++) {
      trace_out[i_iter] = trace_in[i_iter];
    }
  } else {
    calculate_raw_output();
  } 
  
  
}


void General_deghost_engine::calculate_tau() {
  
  // Estimate Tau and show
  int ind = ind_max(ntcheck, check);
  float maxabscheck = fmaxabsf(ntcheck, check);
  float meanabscheck = meanabs(ntcheck, check);
  
  tau = dt*ind;
  nscriteria = maxabscheck/meanabscheck;
  
  
  cout << endl;
  cout << "Tau: " << tau << endl;
//  cout << "tmincheck: " << tmincheck << endl;
//  cout << "tmaxcheck: " << tmaxcheck << endl;
//  cout << "nscriteria: " << nscriteria << endl;
//  cout << endl;
  
}


void General_deghost_engine::calculate_F2_and_denominator() {
 
  
  // if F2 already exists, free the mem
  if ( F2 != NULL ) {
    delete [] F2;
    F2 = NULL;
  }
  
  // if denominator exists, free the mem
  if ( denominator != NULL ) {
    delete [] denominator;
    denominator = NULL;
  }
  
  
  F2 = new complex<float>[nt*ntcheck];
  denominator = new complex<float>[nt*ntcheck];
  
  
  for (int i_iter = 0; i_iter < ntcheck; i_iter++) {
    float tmptau = i_iter*dt;
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      float freq = df*j_iter;
      float omega = 2.0*Pi*freq;
      F2[i_iter*nt + j_iter] = 1.0f - exp(I*omega*tmptau);
      denominator[i_iter*nt + j_iter] = 2.0f - 2.0f*cos(omega*tmptau);
    }
  }
  
  
}


void General_deghost_engine::process(ghost_data * ghostdata) {
  
  int nt_old = nt;
  float dt_old = dt;
  
  // Get info
  nt = ghostdata->get_nt();
  dt = ghostdata->get_dt();
  zmin = ghostdata->get_zmin();
  zmax = ghostdata->get_zmax();
  vwater = ghostdata->get_vwater();
  epsilon = ghostdata->get_epsilon();
  addback = ghostdata->get_addback();
  
  
  // If nt or dt change, we have to recompute F2 and denominator
  if ((nt_old != nt)||(abs(dt-dt_old)>0.001)) {
    recompute = 1;
  }
  
  
  initialize();
  
  ghostdata->get_input_trace(trace_in, nt);
  
  input_average = meanabs(nt, trace_in);
//  display_info();
  
  execute();
  
  ghostdata->set_output_trace(trace_out,nt);
  
  finalize();
  
  
}


void General_deghost_engine::initialize_trace_in_trace_out() {
  
//  cout << "General_deghost_engine::initialize_trace_in_trace_out()" << endl;
  
  
  if (trace_in!=NULL) {
    delete [] trace_in;
    trace_in = NULL;
  }
  trace_in = new float[nt];
  
  if (trace_out!=NULL) {
    delete [] trace_out;
    trace_out = NULL;
  }
  trace_out = new float[nt];
  
  
}


void General_deghost_engine::initialize() {
    
  calculate_basic_info();

    
  // Prepare F2 and denominator
  if ((first_run)||(recompute)) {
    initialize_trace_in_trace_out();
    calculate_F2_and_denominator();
  }  
  
}


void General_deghost_engine::execute() {
 
  
  // Initialize all important arrays
  initialize_Fourier_image();
  initialize_check_array();
  
  
  // Calculate all important arrays
  calculate_Fourier_image();
  calculate_check_array();
  
  
  calculate_tau();
  calculate_output();
  
  
}


void General_deghost_engine::finalize() {
 
  
  first_run = 0;
  recompute = 0;
  
  
  if ( Fourier_image != NULL ) {
    delete [] Fourier_image;
    Fourier_image = NULL;
  } 
  
    
  if ( check != NULL ) {
    delete [] check;
    check = NULL;
  }
    
  
}




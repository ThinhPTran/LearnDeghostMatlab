/* 
 * File:   FxTaup2D.C
 * Author: dragon
 * 
 * Created on December 12, 2013, 11:10 AM
 */

#include "FxTaup2D.H"

FxTaup2D::FxTaup2D() {
  
  nt = 0;
  ntau = 0;
  nx = 0;
  np = 0;
  nf = 0;
  
  dt = 0;
  dtau = 0;
  dx = 0;
  dp = 0;
  df = 0;
  
  tmin = 0;
  xmin = 0;
  pmin = 0;  
  
  T = 0;
  Tau = 0;
  
  mode = UNDEFINED;
  iter_cg = 1000;
  
  x = NULL;
  p = NULL;
  
  in = NULL;
  out = NULL;
  
}


FxTaup2D::~FxTaup2D() {
  
  if (x!=NULL) {
    delete [] x;
    x = NULL;
  }
    
  
  if (p!=NULL) {
    delete [] p;
    p = NULL;
  }
  
}


void FxTaup2D::initialize(TauP2D_info_t info) {
  
  nt = info.nt;
  ntau = info.nt;
  nx = info.nx;
  np = info.np;
  nf = info.nt;
  
  dt = info.dt;
  dtau = info.dt;
  dx = info.dx;
  dp = info.dp;
  df = 1.0/(nt*dt);
  
  tmin = info.tmin;
  xmin = info.xmin;
  pmin = info.pmin;
  
 
  T = nt*dt;
  Tau = ntau*dtau;
  
  in = info.in;
  out = info.out;
  
  
  mode = info.mode;
  
  
  cout << "*********************************************************************" << endl;
  cout << "                        FxTaup2D::initialize                         " << endl;
  cout << "*********************************************************************" << endl;
  cout << "nt: " << nt << endl;
  cout << "ntau: " << ntau << endl;
  cout << "nx: " << nx << endl;
  cout << "np: " << np << endl;
  cout << "nf: " << nf << endl;
  cout << "dt: " << dt << endl;
  cout << "dtau: " << dtau << endl;
  cout << "dx: " << dx << endl;
  cout << "dp: " << dp << endl;
  cout << "df: " << df << endl;
  cout << "tmin: " << tmin << endl;
  cout << "xmin: " << xmin << endl;
  cout << "pmin: " << pmin << endl;
  cout << "T: " << T << endl;
  cout << "Tau: " << Tau << endl;
  switch(mode) {
    case FWD_TAUP:
      cout << "mode: " << "Forward Taup Transform." << endl;
      break;
    case INV_TAUP: 
      cout << "mode: " << "Inverse Taup Transform." << endl;
      break;
    default:
      cout << "Unidentified mode: " << mode << endl;
      exit(-1);
      break;
  }
//  cout << "x: " << x << endl;
//  cout << "p: " << p << endl;
  cout << "*********************************************************************" << endl;
  
  
  initialize_x();
  initialize_p();
  
  
}


void FxTaup2D::execute() {

  switch(mode) {
    case FWD_TAUP:
      forward_fx_taup_transform();
      break;
    case INV_TAUP: 
      inverse_fx_taup_transform();
      break;
    default:
      cout << "Unidentified mode: " << mode << endl;
      exit(-1);
      break;
  }
  
}


void FxTaup2D::forward_fx_taup_transform() {
  
  int hf = floor(nf/2);
  complex<double> * frequency_slice_in = new complex<double>[nx*nf];
  complex<double> * frequency_slice_out = new complex<double>[np*nf];
  
  
  transf_from_t_to_f(frequency_slice_in);
  
  
  // Loop over frequency slices 
#pragma omp parallel for
  for (int i_iter = 0; i_iter < hf; i_iter++) {
    
    complex<double> * tmp_slice_in = new complex<double>[nx];
    complex<double> * tmp_slice_out = new complex<double>[np];
    complex<double> * A = new complex<double>[np*nx];
    
    
    // Take one frequency slice
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      tmp_slice_in[j_iter] = frequency_slice_in[j_iter*nf + i_iter];  
    }
    
    
    for (int j_iter = 0; j_iter < np; j_iter++) {
      tmp_slice_out[j_iter] = 0.0;
    }
    
    
    // Recalculating A
    recalculating_A(i_iter, A);
    
    
    // Multiply A by tmp_slice_in
    multiply_A_by_freq_slice(tmp_slice_in, tmp_slice_out, A);
    
    
    // Write the output to output array
    for (int j_iter = 0; j_iter < np; j_iter++) {
      frequency_slice_out[j_iter*nf + i_iter] = tmp_slice_out[j_iter];
    }
    
    
    delete [] tmp_slice_out;
    delete [] tmp_slice_in;
    delete [] A;
  
  }
  
  
  // Symmetrize the output frequency slices
  for (int i_iter = 1; i_iter < hf; i_iter++) {
    for (int j_iter = 0; j_iter < np ; j_iter++) {
      frequency_slice_out[j_iter*nf + nf - i_iter] 
              = conj(frequency_slice_out[j_iter*nf + i_iter]); 
    }
  }
  
  
  
  transf_from_f_to_tau(frequency_slice_out); 
  
  
  delete [] frequency_slice_out;
  delete [] frequency_slice_in; 

  
}


void FxTaup2D::inverse_fx_taup_transform() {
  
  int hf = floor(nf/2);
  complex<double> * frequency_slice_in = new complex<double>[np*nf];
  complex<double> * frequency_slice_out = new complex<double>[nx*nf];
  
  
  for (int i_iter = 0; i_iter < nx*nf; i_iter++) {
    frequency_slice_out[i_iter] = 0.0;
  }
  
  
  transf_from_tau_to_f(frequency_slice_in);
  
  
  // Loop over frequency slices 
#pragma omp parallel for
  for (int i_iter = 0; i_iter < hf; i_iter++) {
    
    complex<double> * tmp_slice_in = new complex<double>[np];
    complex<double> * tmp_slice_out = new complex<double>[nx];
    complex<double> * A = new complex<double>[np*nx];
    
    
    // Take one frequency slice
    for (int j_iter = 0; j_iter < np; j_iter++) {
      tmp_slice_in[j_iter] = frequency_slice_in[j_iter*nf + i_iter];  
    }
    
    
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      tmp_slice_out[j_iter] = 0.0;
    }
    
    
    // Recalculating A
    recalculating_A(i_iter, A);
    
    
    // Solve a linear system
    linear_sys_solver(tmp_slice_in, tmp_slice_out, A);
    
    
    
    // Write the output to output array
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      frequency_slice_out[j_iter*nf + i_iter] = tmp_slice_out[j_iter];
    }
    
    
    delete [] tmp_slice_out;
    delete [] tmp_slice_in;
    delete [] A;
  
  }
  
  
//  // Symmetrize the output frequency slices
//  for (int i_iter = 1; i_iter < hf; i_iter++) {
//    for (int j_iter = 0; j_iter < nx ; j_iter++) {
//      frequency_slice_out[j_iter*nf + nf - i_iter] 
//              = conj(frequency_slice_out[j_iter*nf + i_iter]); 
//    }
//  }
  
  
  transf_from_f_to_t(frequency_slice_out);
  
  
  delete [] frequency_slice_out;
  delete [] frequency_slice_in; 
  
}


void FxTaup2D::finalize() {
  
  if (x!=NULL) {
    delete [] x;
    x = NULL;
  }
    
  
  if (p!=NULL) {
    delete [] p;
    p = NULL;
  } 
  
}


void FxTaup2D::initialize_x() {
  
  if (x==NULL) {
    x = new double[nx];
  } else {
    delete [] x;
    x = NULL;
    x = new double[nx];
  }
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    x[i_iter] = xmin + i_iter*dx;
  }
  
}


void FxTaup2D::initialize_p() {
  
  if (p==NULL) {
    p = new double[np];
  } else {
    delete [] p;
    p = NULL;
    p = new double[np];
  }
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    p[i_iter] = pmin + i_iter*dp;
  }
  
}


void FxTaup2D::uniform_fftw1d
(
int length,
fftw_complex * in,
fftw_complex * out,
int sign
)
{

  int i_iter;
  fftw_plan p;

  // Preprocessing check
  if ((sign != FFTW_BACKWARD) && (sign != FFTW_FORWARD)) {
    printf("Wrong option of sign!!!\n");
    exit(-1);
  }
  
  
#pragma omp critical
  p = fftw_plan_dft_1d(length, in , out, sign, FFTW_ESTIMATE);
  fftw_execute(p);

  // Post-processing check
  switch (sign) {
    case FFTW_BACKWARD:
      for (i_iter = 0; i_iter < length; i_iter++) {
        out[i_iter][0] = out[i_iter][0]/(1.0*length);
        out[i_iter][1] = out[i_iter][1]/(1.0*length);
      }
      break;

    case FFTW_FORWARD:
      break;
  }

#pragma omp critical
  fftw_destroy_plan(p);

}


void FxTaup2D::transf_from_t_to_f(complex<double>* output) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nf*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = in[i_iter*nt+j_iter];
      tmp_in[j_iter][1] = 0;
    }
    
    
    uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_FORWARD);
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      output[i_iter*nf + j_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
    }
    
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void FxTaup2D::transf_from_tau_to_f(complex<double>* output) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nf*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = in[i_iter*nt+j_iter];
      tmp_in[j_iter][1] = 0;
    }
    
    
    uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_FORWARD);
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      output[i_iter*nf + j_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
    }
    
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void FxTaup2D::transf_from_f_to_tau(complex<double>* input) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nf*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = input[i_iter*nt+j_iter].real();
      tmp_in[j_iter][1] = input[i_iter*nt+j_iter].imag();
    }
    
    
    uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_BACKWARD);
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nf + j_iter] = tmp_out[j_iter][0];
    }
    
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void FxTaup2D::transf_from_f_to_t(complex<double>* input) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nf*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = input[i_iter*nt+j_iter].real();
      tmp_in[j_iter][1] = input[i_iter*nt+j_iter].imag();
    }
    
    
    uniform_fftw1d(nt, tmp_in, tmp_out, FFTW_BACKWARD);
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nf + j_iter] = tmp_out[j_iter][0];
    }
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void FxTaup2D::recalculating_A(int freq_idx, complex<double> * A) {
  
  double fid = freq_idx;
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    for (int j_iter = 0; j_iter < np; j_iter++) {
      A[i_iter*np + j_iter] = exp(I*2.0*MY_PI*fid*df*x[i_iter]*p[j_iter]);
    }
  }
  
}


void FxTaup2D::multiply_A_by_freq_slice
(
  complex<double>* in, 
  complex<double>* out,
  complex<double>* A
) 
{
  
  memset((void*)out, 0x00, np*sizeof(complex<double>));
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      out[i_iter] += A[j_iter*np+i_iter]*in[j_iter]; 
    }
  }
  
  
}


void FxTaup2D::linear_sys_solver
(
complex<double>* input,
complex<double>* output,
complex<double>* A
)
{
  
  int count = 0;
  

  complex<double> * S = new complex<double>[np];
  
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    S[i_iter] = input[i_iter];
  }
 

  double f = 0.0;
  double rms = 0.0;
  double anpha = 0.0;
  double beta = 0.0;
  double rr_old = 0.0;
  double rr_new = 0.0;
  double err_old = 0.0;
  double err_new = 0.0;
  complex<double> * tmp_r = new complex<double>[nx];
  complex<double> * tmp_p = new complex<double>[nx];
  complex<double> * tmp_q = new complex<double>[np];


  // r = b - AU = b - AO = F^h[u]
  adjoint_operator(S, tmp_r, A);
  
  // p = r
  for (int i_iter = 0; i_iter <   nx; i_iter++) {
    tmp_p[i_iter] = tmp_r[i_iter];
  }

  
  //rr_new = norm_vec(r)
  rr_new = norm_vec(tmp_r, nx);
  err_old = norm_vec(S, np);


  for (int i_iter = 0; i_iter < iter_cg; i_iter++) {
    
    count = i_iter;
    
    
    // q = F[p]
    forward_operator(tmp_p, tmp_q, A);

    // anpha = rr_new/||q||^2
    anpha = rr_new/(norm_vec(tmp_q, np));

    
    // U = U + anpha*p
    for (int j_iter = 0; j_iter <   nx; j_iter++) {
      output[j_iter] += anpha*tmp_p[j_iter];
    }
   

    // u = u - anpha*q
    for (int j_iter = 0; j_iter < np; j_iter++) {
      S[j_iter] -= anpha*tmp_q[j_iter];
    }
    

    // r = F^h[u]
    adjoint_operator(S, tmp_r, A);

    // rr_old = rr_new
    rr_old = rr_new;

    // rr_new = ||r||^2
    rr_new = norm_vec(tmp_r,nx);

    // beta = rr_new/rr_old
    beta = rr_new/rr_old;

    // p = r + beta.p
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      tmp_p[j_iter] = tmp_r[j_iter] + beta*tmp_p[j_iter];
    }

    // Error term computation
    err_new = norm_vec(S, np);
    f = 200.0*abs(err_old - err_new)/(err_old + err_new);
    rms = sqrt(err_new);
    err_old = err_new;

    
    // if rms < 0.001; break; end
    if ( rms < 0.001 ) break;
    if ( f < 0.01 ) break;

  }


//  cout << "       " << "iter_cg: " << count << endl;
//  cout << "           " << "rms: " << rms << endl;
//  cout << "             " << "f: " << f << endl;
  
  delete [] tmp_r;
  delete [] tmp_p;
  delete [] tmp_q;
  delete [] S;

  
}


void FxTaup2D::forward_operator
(
  complex<double> * input,
  complex<double> * output,
  complex<double>* A
)
{

  memset((void*)output,0x00,np*sizeof(complex<double>));
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      output[i_iter] += (A[j_iter*np + i_iter]*input[j_iter]);
    }
  }
  
}


void FxTaup2D::adjoint_operator
(
  complex<double> * input,
  complex<double> * output,
  complex<double> * A
)
{

  memset((void*)output,0x00,nx*sizeof(complex<double>));
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    for (int j_iter = 0; j_iter < np; j_iter++) {
      output[i_iter] += (conj(A[i_iter*np + j_iter])*input[j_iter]);
    }
  }
  
}


double FxTaup2D::norm_vec
(
complex<double> * vec,
int vec_length
)
{
  double result = 0.0;
  for (int i_iter = 0; i_iter < vec_length; i_iter++) {
    result += real(vec[i_iter]*conj(vec[i_iter]));
  }

  return result;
}



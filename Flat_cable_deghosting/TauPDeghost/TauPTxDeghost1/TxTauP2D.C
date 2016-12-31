/* 
 * File:   TxTaup2D.C
 * Author: dragon
 * 
 * Created on December 12, 2013, 11:10 AM
 */

#include "TxTauP2D.H"

TxTaup2D::TxTaup2D() {
  
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
  domega = 0;
  
  tmin = 0;
  xmin = 0;
  pmin = 0;  
  fmax = 0;
  iwmax = 0;
  
  T = 0;
  Tau = 0;
  
  mode = UNDEFINED;
  iter_cg = 1000;
  
  omega = NULL;
  x = NULL;
  p = NULL;
  
  filter = NULL;
  
  in = NULL;
  out = NULL;
  
}


TxTaup2D::~TxTaup2D() {
  
  if (omega!=NULL) {
    delete [] omega;
    omega = NULL;
  }
  
  if (x!=NULL) {
    delete [] x;
    x = NULL;
  }   
  
  if (p!=NULL) {
    delete [] p;
    p = NULL;
  }
  
  if (filter!=NULL) {
    delete [] filter;
    filter = NULL;
  }
  
}


void TxTaup2D::calculate_dependent_parms() {
  
  T = nt*dt;
  Tau = ntau*dtau;
  domega = 2.0*M_PI*df;
  iwmax = int(fmax/df);
  
}


void TxTaup2D::initialize(TauP2D_info_t info) {
  
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
  fmax = info.fmax;
  
  in = info.in;
  out = info.out;
  
  mode = info.mode;
  
  
  calculate_dependent_parms();
  
  
  cout << "*********************************************************************" << endl;
  cout << "                        TxTaup2D::initialize                         " << endl;
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
  cout << "domega: " << domega << endl;
  cout << "tmin: " << tmin << endl;
  cout << "xmin: " << xmin << endl;
  cout << "pmin: " << pmin << endl;
  cout << "fmax: " << fmax << endl;
  cout << "iwmax: " << iwmax << endl;
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
  initialize_omega();
  initialize_filter();
  
  
}


void TxTaup2D::execute() {

  switch(mode) {
    case FWD_TAUP:
      forward_tx_taup_transform();
      break;
    case INV_TAUP: 
      inverse_tx_taup_transform();
      break;
    default:
      cout << "Unidentified mode: " << mode << endl;
      exit(-1);
      break;
  }
  
}


void TxTaup2D::forward_tx_taup_transform() {

  /* loop over slopes */
  for (int ip=0; ip<np; ip++) {

    /* initialize output array */
    for (int it=0; it<nt; it++)
            out[ip*nt+it]=0.0;

    /* loop over traces */
    for (int ix=0; ix<nx; ix++) {
      
      int fit,lit;            /* first and last time samples */
      int id;                 /* auxiliary variables */
      float dfrac,delay;      /* more auxiliary variables */

      /* compute two point interpolator */
      delay=p[ip]*x[ix]/dt;
      if (delay>=0) {
              id = (int)delay;
              fit = id+1;
              lit = nt-1;
      } else {
              id = (int)delay-1;
              fit = 1;
              lit = nt+id;
      }
      dfrac = delay-id;

      /* compute the actual slant stack */
      for (int it=fit; it<lit; it++) {
              out[ip*nt+it-id]+=fabs(dx)*
                  (in[ix*nt+it]+
                   dfrac*(in[ix*nt+it+1]-in[ix*nt+it]));
      }
    }
  }
  
}


void TxTaup2D::inverse_slant_stack() {
  
  /* loop over slopes */
  for (int ix=0; ix<nx; ix++) {

    /* initialize output array */
    for (int it=0; it<nt; it++)
            out[ix*nt+it]=0.0;

    /* loop over traces */
    for (int ip=0; ip<np; ip++) {
      
      int fit,lit;            /* first and last time samples */
      int id;                 /* auxiliary variables */
      float dfrac,delay;      /* more auxiliary variables */

      /* compute two point interpolator */
      delay=-p[ip]*x[ix]/dt;
      if (delay>=0) {
              id = (int)delay;
              fit = id+1;
              lit = nt-1;
      } else {
              id = (int)delay-1;
              fit = 1;
              lit = nt+id;
      }
      dfrac = delay-id;

      /* compute the actual slant stack */
      for (int it=fit; it<lit; it++) {
              out[ix*nt+it-id]+=fabs(dp)*
                  (in[ip*nt+it]+
                   dfrac*(in[ip*nt+it+1]-in[ip*nt+it]));
      }
    }
  }
  
}


void TxTaup2D::filter_one_trace(complex<double>* trace) {
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    trace[i_iter] *= filter[i_iter];
  }
  
}


void TxTaup2D::apply_rhofilter() {
  
  complex<double> * tmp_out = new complex<double>[nx*nt];
  
  transf_from_t_to_f(out, tmp_out, nx, nt);
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    filter_one_trace(&tmp_out[i_iter*nt]);
  }
  
  transf_from_f_to_t(tmp_out, out, nx, nt);
  
  delete [] tmp_out;
  
}


void TxTaup2D::inverse_tx_taup_transform() {
  
  inverse_slant_stack(); 
  apply_rhofilter();
  
}


void TxTaup2D::finalize() {
  
  if (omega!=NULL) {
    delete [] omega;
    omega = NULL;
  }
  
  if (x!=NULL) {
    delete [] x;
    x = NULL;
  }   
  
  if (p!=NULL) {
    delete [] p;
    p = NULL;
  }
  
  if (filter!=NULL) {
    delete [] filter;
    filter = NULL;
  } 
  
}


void TxTaup2D::initialize_x() {
  
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


void TxTaup2D::initialize_p() {
  
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


void TxTaup2D::initialize_omega() {
  
  int hnf = floor(nt/2) + 1; 
  
  if (omega==NULL) {
    omega = new double[nf];
  } else {
    delete [] omega;
    omega = NULL;
    omega = new double[nf];
  }
  
  omega[0] = 0;
  
  for (int i_iter = 1; i_iter < hnf; i_iter++) {
    omega[i_iter] = i_iter*domega;
    omega[nt-i_iter] = -omega[i_iter];
  }
  
}


void TxTaup2D::initialize_filter() {
  
  int hnf = floor(nf/2);
  double wmax = iwmax*domega;
  double wnyq = hnf*domega;
  double wscale;
  
  if (wnyq > wmax) {
    wscale = wmax / (wnyq-wmax);
  } else {
    wscale = 1.0;
  }
  
  
  if (filter==NULL) {
    filter = new double[nf];
  } else {
    delete [] filter;
    filter = NULL;
    filter = new double[nf];
  }
  
  
  
  filter[0] = 0;
  
  for (int i_iter = 1; i_iter <= hnf; i_iter++) {
    if (i_iter < iwmax) {
      filter[i_iter] = omega[i_iter];
      filter[nf-i_iter] = filter[i_iter];
    } else {
      filter[i_iter] = wscale*(wnyq-omega[i_iter]);
      filter[nf-i_iter] = filter[i_iter];
    }
  }
  
  
}


void TxTaup2D::uniform_fftw1d
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


void TxTaup2D::uniform_forward_fftw1d_drc
(
  int length, 
  double* in, 
  complex<double>* out
) 
{
  
  fftw_plan p; 
  fftw_complex * tmp_in = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  fftw_complex * tmp_out = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  
  // initialize input trace
  for (int i_iter = 0; i_iter < length; i_iter++) {
    tmp_in[i_iter][0] = in[i_iter];
    tmp_in[i_iter][1] = 0;
  }
  
#pragma omp critical
  p = fftw_plan_dft_1d(length, tmp_in , tmp_out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  
  // output data
  for (int i_iter = 0; i_iter < length; i_iter++) {
    out[i_iter] = complex<double>(tmp_out[i_iter][0], tmp_out[i_iter][1]);
  }
  

#pragma omp critical
  fftw_destroy_plan(p);
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TxTaup2D::uniform_backward_fftw1d_dcr 
(
  int length,
  complex<double> * in,
  double * out
) 
{
  
  fftw_plan p; 
  fftw_complex * tmp_in = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  fftw_complex * tmp_out = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  
  // initialize input trace
  for (int i_iter = 0; i_iter < length; i_iter++) {
    tmp_in[i_iter][0] = in[i_iter].real();
    tmp_in[i_iter][1] = in[i_iter].imag();
  }
  
#pragma omp critical
  p = fftw_plan_dft_1d(length, tmp_in , tmp_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  
  // output data
  for (int i_iter = 0; i_iter < length; i_iter++) {
    out[i_iter] = tmp_out[i_iter][0];
  }
  

#pragma omp critical
  fftw_destroy_plan(p);
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TxTaup2D::transf_from_t_to_f(double * tracein, complex<double>* traceout, int nx, int nt) {
  
//#pragma omp parallel for
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    uniform_forward_fftw1d_drc(nt, &tracein[i_iter*nt],&traceout[i_iter*nt]);
  }
  
}


void TxTaup2D::transf_from_f_to_t(complex<double>* tracein, double* traceout, int nx, int nt) {
  
//#pragma omp parallel for
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    uniform_backward_fftw1d_dcr(nt, &tracein[i_iter*nt], &traceout[i_iter*nt]);
  }
  
}


double TxTaup2D::norm_vec
(
complex<double> * vec,
int vec_length
)
{
  
  double result = 0.0;
  for (int i_iter = 0; i_iter < vec_length; i_iter++) {
    result += norm(vec[i_iter]);
  }

  return result;
}



/*
 * =====================================================================================
 *
 *       Filename:  FFTW_auxfuncs.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/16/2015 10:05:46 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thinh P. Tran (), tpthinh@fairfield.com.vn, trphthinh@gmail.com
 *   Organization:  Fairfield Vietnam Ltd.
 *
 * =====================================================================================
 */

#include "FFTW_auxfuncs.H"



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  uniform_fftw1d
 *  Description:  My first aux. for using fftw function. It has been extensively used 
 *                on my old classes but there is still an unconvenience of malloc mem
 *                before using this function. This function will be considered obsolete
 *                in near future.
 *   Parameters: 
 *       length:  (int) the length of input and output vector.
 *           in:  (fftw_complex *) input vector.
 *          out:  (fftw_complex *) output vector.
 *         sign:  (int) Forward or Inverse Transform. 
 *
 * =====================================================================================
 */

void uniform_fftw1d
(
int length,
fftw_complex * in,
fftw_complex * out,
int sign
)
{

  int i_iter;
  fftw_plan p;

  //Preprocessing check
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

}		/* -----  end of function uniform_fftw1d  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  uniform_forward_fftw1d_drc
 *  Description:  Forward FFTW transform
 *
 *   Parameters: 
 *       length:  (int) the length of input and output vector.
 *           in:  (double *) input vector.
 *          out:  (complex<double> *) output vector.
 *
 * =====================================================================================
 */

void uniform_forward_fftw1d_drc
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

}		/* -----  end of function uniform_forward_fftw1d_drc  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  uniform_backward_fftw1d_dcr
 *  Description:  Backward FFTW transform
 *
 *   Parameters: 
 *       length:  (int) the length of input and output vector.
 *           in:  (complex<double> *) input vector.
 *          out:  (double *) output vector.
 *
 * =====================================================================================
 */


void uniform_backward_fftw1d_dcr
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
 
}		/* -----  end of function uniform_backward_fftw1d_dcr  ----- */



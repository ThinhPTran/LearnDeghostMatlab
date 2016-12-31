/* 
 * File:   TauPDeghost.C
 * Author: dragon
 * 
 * Created on May 7, 2014, 8:37 AM
 */

#include "TauP2DDeghost.H"


TauP2DDeghost::TauP2DDeghost() {
  
  // Input info
  vwater = 0;
  recover_eps = 0;
  src_depth = 0;
  rec_depth = 0;
  
  nt = 0;
  nx = 0;
  
  dt = 0;
  dx = 0;
  dp = 0;
  
  fx = 0;
  
  in = NULL;
  out = NULL;
  
  
  // Calculated info
  np = 0;
  nkx = 0;
  exnt = 0;
  
  dkx = 0;
  exdt = 0;
  exdf = 0;
  
  ft = 0;
  fp = 0;
  fomega = 0;
  
  
  p = NULL;
  kx = NULL;
  exomega = NULL;
  phys_filter = NULL;
  exdata = NULL;
  
  
  
}


TauP2DDeghost::~TauP2DDeghost() {
  
  // Input info
  vwater = 0;
  recover_eps = 0;
  src_depth = 0;
  rec_depth = 0;
  
  nt = 0;
  nx = 0;
  
  dt = 0;
  dx = 0;
  dp = 0;
  
  fx = 0;
  
  in = NULL;
  out = NULL;
    
  
  // Calculated info
  np = 0;
  nkx = 0;
  exnt = 0;
  
  dkx = 0;
  exdt = 0;
  exdf = 0;
  
  ft = 0;
  fp = 0;
  fomega = 0;
  
  
  p = NULL;
  kx = NULL;
  exomega = NULL;
  phys_filter = NULL;
  exdata = NULL;
  
  
}


void TauP2DDeghost::get_input_info(TauP2DDeghost_info_t info) {
  
  vwater = info.vwater;
  recover_eps = info.recover_eps;
  src_depth = info.src_depth;
  rec_depth = info.rec_depth;
  
  nt = info.nt;
  nx = info.nx;
  dt = info.dt;
  dx = info.dx;
  dp = info.dp;
  fx = info.fx;
  in = info.in;
  out = info.out;
  
}


void TauP2DDeghost::calculate_dependent_parameters() {
  
  np = floor(2.0/(vwater*dp)) + 1;
  if (!(np%2)) {
    np = np+1;
  }
  
  nkx = nx;
  exnt = 2*nt;
  
  dkx = 2.0*MY_PI/((nx*dx));
  exdt = dt;
  exdf = 1.0/(exnt*exdt);
  exdomega = 2.0*MY_PI*exdf;
  
  ft = 0.0;
  fp = -floor(np/2)*dp;
  fomega = 0.0;
  
  pmax = floor(np/2)*dp;
  
}


void TauP2DDeghost::display_info_for_checking() {
  
  cout << "*********************************************************************" << endl;
  cout << "                    TauP2DDeghost::Initialize                        " << endl;
  cout << "*********************************************************************" << endl;
  cout << "vwater: " << vwater << endl;
  cout << "recover_eps: " << recover_eps << endl;
  cout << "src depth: " << src_depth << endl;
  cout << "rec depth: " << rec_depth << endl;
  cout << "nt: " << nt << endl;
  cout << "nx: " << nx << endl;
  cout << "nkx: " << nkx << endl;
  cout << "np: " << np << endl;
  cout << "exnt: " << exnt << endl;
  cout << "dt: " << dt << endl;
  cout << "dx: " << dx << endl;
  cout << "dkx: " << dkx << endl;
  cout << "dp: " << dp << endl;
  cout << "exdt: " << exdt << endl;
  cout << "exdf: " << exdf << endl;
  cout << "exdomega: " << exdomega << endl;
  cout << "ft: " << ft << endl;
  cout << "fx: " << fx << endl;
  cout << "fp: " << fp << endl;
  cout << "fomega: " << fomega << endl;
  cout << "pmax: " << pmax << endl;
//  cout << "in: " << in << endl;
//  cout << "out: " << out << endl;
  cout << "*********************************************************************" << endl;
  
}


void TauP2DDeghost::allocate_main_mem() {
  
  p = new double[np];
  kx = new double[nx];
  exomega = new double[exnt];
  phys_filter = new double[exnt*nx];
  exdata = new complex<double>[nx*exnt];
  
}


void TauP2DDeghost::compute_p() {
  
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    p[i_iter] = fp + i_iter*dp;
  }
  
  
}


void TauP2DDeghost::compute_kx() {
  
  
  for (int i_iter = 0; i_iter <= floor(nx/2); i_iter++) {
    kx[i_iter] = i_iter*dkx;
  }

  for (int i_iter = floor(nx/2) + 1; i_iter < nx; i_iter++) {
    kx[i_iter] = -kx[nx-i_iter];
  }
  
  
}


void TauP2DDeghost::compute_exomega() {
  
    
  for (int i_iter = 0; i_iter <= floor(exnt/2); i_iter++) {
    exomega[i_iter] = i_iter*exdomega;
  }

  for (int i_iter = floor(exnt/2) + 1; i_iter < exnt; i_iter++) {
    exomega[i_iter] = -exomega[exnt-i_iter];
  }
 
  
}


void TauP2DDeghost::compute_phys_filter() {
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    for (int j_iter = 0; j_iter < nkx; j_iter++) {
      double tmp = (exomega[i_iter]*exomega[i_iter])/(vwater*vwater) - kx[j_iter]*kx[j_iter];
      if (tmp <= 0) {
        phys_filter[j_iter*exnt+i_iter] = 0.0;
      } else {
        phys_filter[j_iter*exnt+i_iter] = 1.0;
      }
    }  
  }
  
}


void TauP2DDeghost::precompute_variables() {
  
  compute_p();
  compute_kx();
  compute_exomega();
  compute_phys_filter();
  
}


void TauP2DDeghost::get_exdata() {
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      exdata[i_iter*exnt+j_iter] = in[i_iter*nt+j_iter];
    }
  }
  
}


void TauP2DDeghost::output_data() {
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nt+j_iter] = exdata[i_iter*exnt+j_iter].real();
    }
  }
  
}


void TauP2DDeghost::biwrite(string file_name, double * data, int M, int nt) {
  
  ofstream myfile;
  int num_elems = M*nt;
  int data_size = num_elems*sizeof(float);
  float * memblock = new float[num_elems];
  
  
  for (int i_iter = 0; i_iter < num_elems; i_iter++) {
    memblock[i_iter] = data[i_iter];
  }
  
  
  myfile.open(file_name.c_str(), ios::binary | ios::out);
  
  if (myfile.is_open()) {
    myfile.write((char*)memblock, data_size);
    myfile.close();
  } else {
    cout << "Can't open file " << file_name << endl;
    exit(-1);
  }
  
  
  delete [] memblock;
  
}


void TauP2DDeghost::biwrite_complex(string file_name, complex<double> * data, int M, int nt) {
  
  ofstream myfile;
  int num_elems = M*nt;
  int data_size = num_elems*sizeof(float);
  float * tmp_float = new float[num_elems];
  
  
  for (int i_iter = 0; i_iter < num_elems; i_iter++) {
    tmp_float[i_iter] = abs(data[i_iter]);
  }
  
  
  myfile.open(file_name.c_str(), ios::binary | ios::out);
  
  if (myfile.is_open()) {
    myfile.write((char*)tmp_float, data_size);
    myfile.close();
  } else {
    cout << "Can't open file " << file_name << endl;
    exit(-1);
  }
  
  delete [] tmp_float;
  
}


void TauP2DDeghost::uniform_fftwf1d(int Nt, fftw_complex * in, fftw_complex * out, int sign) {
  
  int i_iter;
  fftw_plan p;

  // Preprocessing check
  if ((sign != FFTW_BACKWARD) && (sign != FFTW_FORWARD)) {
    printf("Wrong option of sign!!!\n");
    exit(-1);
  }
  
  
#pragma omp critical
  p = fftw_plan_dft_1d(Nt, in , out, sign, FFTW_ESTIMATE);
  fftw_execute(p);

  // Post-processing check
  switch (sign) {
    case FFTW_BACKWARD:
      for (i_iter = 0; i_iter < Nt; i_iter++) {
        out[i_iter][0] = out[i_iter][0]/(1.0*Nt);
        out[i_iter][1] = out[i_iter][1]/(1.0*Nt);
      }
      break;

    case FFTW_FORWARD:
      break;
  }

#pragma omp critical
  fftw_destroy_plan(p);
  
}


void TauP2DDeghost::uniform_fftwf2d(int nt, int nx, complex<double>* in, complex<double>* out) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = in[i_iter*nt+j_iter].real();
      tmp_in[j_iter][1] = in[i_iter*nt+j_iter].imag();
    }
    
    
    uniform_fftwf1d(nt, tmp_in, tmp_out, FFTW_FORWARD);
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nt + j_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
    }
    
  }
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
  
  tmp_in = (fftw_complex *) fftw_malloc(nx*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nx*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      tmp_in[j_iter][0] = out[j_iter*nt+i_iter].real();
      tmp_in[j_iter][1] = out[j_iter*nt+i_iter].imag();
    }
    
    
    uniform_fftwf1d(nx, tmp_in, tmp_out, FFTW_FORWARD);
    
    
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      out[j_iter*nt + i_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
    }
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TauP2DDeghost::uniform_ifftwf2d(int nt, int nx, complex<double>* in, complex<double>* out) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_in[j_iter][0] = in[i_iter*nt+j_iter].real();
      tmp_in[j_iter][1] = in[i_iter*nt+j_iter].imag();
    }
    
    
    uniform_fftwf1d(nt, tmp_in, tmp_out, FFTW_BACKWARD);
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nt + j_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
    }
    
  }
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
  
  tmp_in = (fftw_complex *) fftw_malloc(nx*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nx*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    
    // Take one trace
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      tmp_in[j_iter][0] = out[j_iter*nt+i_iter].real();
      tmp_in[j_iter][1] = out[j_iter*nt+i_iter].imag();
    }
    
    
    uniform_fftwf1d(nx, tmp_in, tmp_out, FFTW_BACKWARD);
    
    
    for (int j_iter = 0; j_iter < nx; j_iter++) {
      out[j_iter*nt + i_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]).real();
    }
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TauP2DDeghost::initialize(TauP2DDeghost_info_t info) {
  
  // Input info
  get_input_info(info);
  
  
  // Calculate dependent parameters
  calculate_dependent_parameters();
  
  
  // Display info for checking
  display_info_for_checking();
  
  
  // Allocate main memory
  allocate_main_mem();
  
  
  // Precompute variables
  precompute_variables();
  
  
  // Get exinput
  get_exdata();
//  biwrite_complex("exinput", exdata, nx, exnt);

  
}


void TauP2DDeghost::phys_filtering() {
  
  int M = exnt*nx;
  complex<double> * exdata_fk = new complex<double>[exnt*nx];
  
  uniform_fftwf2d(exnt, nx, exdata, exdata_fk);
  
  // Keep symmetry
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata_fk[i_iter] = exdata_fk[i_iter]*phys_filter[i_iter];
  }
  
  
  for (int i_iter = 1; i_iter<nt; i_iter++) { 
    exdata_fk[exnt-i_iter] = conj(exdata_fk[i_iter]);
  }


  for (int i_iter = 1; i_iter < nt; i_iter++) {
    for (int j_iter = 1; j_iter < nkx; j_iter++) {
      exdata_fk[(nkx-j_iter)*exnt+exnt-i_iter] = conj(exdata_fk[j_iter*exnt+i_iter]);
    }
  }
  
  uniform_ifftwf2d(exnt, nx, exdata_fk, exdata);
//  biwrite_complex("exdata_phys_filter",exdata,nx, exnt);
  
  
  delete [] exdata_fk;
  
}


void TauP2DDeghost::fwd_taup_transform
(
  int nt, int nx, int np, double dt, double dx, double dp, 
  double tmin, double xmin, double pmin, double * in, double * out
) 
{
  
  TauP2D_info_t info;
  FxTaup2D taup_engine;
  
  // Forward TauP transform
  info.nt = nt;
  info.nx = nx;
  info.np = np;
  info.dt = dt;
  info.dx = dx;
  info.dp = dp;
  info.tmin = tmin;
  info.xmin = xmin;
  info.pmin = pmin;
  info.mode = FWD_TAUP;
  info.in = in;
  info.out = out;
  
  
  taup_engine.initialize(info);
  taup_engine.execute();
  taup_engine.finalize();
  
}


void TauP2DDeghost::bwd_taup_transform
(
  int nt, int nx, int np, double dt, double dx, double dp, 
  double tmin, double xmin, double pmin, double * in, double * out
) 
{
  
  TauP2D_info_t info;
  FxTaup2D taup_engine;
  
  // Forward TauP transform
  info.nt = nt;
  info.nx = nx;
  info.np = np;
  info.dt = dt;
  info.dx = dx;
  info.dp = dp;
  info.tmin = tmin;
  info.xmin = xmin;
  info.pmin = pmin;
  info.mode = INV_TAUP;
  info.in = in;
  info.out = out;
  
  
  taup_engine.initialize(info);
  taup_engine.execute();
  taup_engine.finalize();
  
}


void TauP2DDeghost::taup_deghost_tau_to_f(double* taup_img, complex<double>* ftaup_img) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    double test = 0.8 - (p[i_iter]*p[i_iter]*vwater*vwater);
    
    if (test > 0) {
    
      // Take one trace
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        tmp_in[j_iter][0] = taup_img[i_iter*exnt+j_iter];
        tmp_in[j_iter][1] = 0.0;
      }


      uniform_fftwf1d(exnt, tmp_in, tmp_out, FFTW_FORWARD);


      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        ftaup_img[i_iter*exnt + j_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
      }
    
    } else {
      
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        taup_img[i_iter*exnt+j_iter] = 0.0;
        ftaup_img[i_iter*exnt + j_iter] = 0.0;
      }
      
    }
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TauP2DDeghost::taup_deghost_f_to_tau(complex<double>* ftaup_img, double* taup_img) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
  
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    double test = 0.8 - (p[i_iter]*p[i_iter]*vwater*vwater);
    
    if (test > 0) {
    
//      for (int j_iter = 1; j_iter < floor(exnt/2); j_iter++) {
//        ftaup_img[i_iter*exnt+exnt-j_iter] = conj(ftaup_img[i_iter*exnt+j_iter]);
//      }

      // Take one trace
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        tmp_in[j_iter][0] = ftaup_img[i_iter*exnt+j_iter].real();
        tmp_in[j_iter][1] = ftaup_img[i_iter*exnt+j_iter].imag();
      }


      uniform_fftwf1d(exnt, tmp_in, tmp_out, FFTW_BACKWARD);


      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        taup_img[i_iter*exnt + j_iter] = tmp_out[j_iter][0];
      }
    
    }
    
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TauP2DDeghost::taup_deghost(complex<double> * ftaup_img, double depth) {
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    double test = 1.0 - (p[i_iter]*p[i_iter]*vwater*vwater);
    
    if (test > 0) {
      
      double tau = 2.0*depth*sqrt(1.0-(p[i_iter]*p[i_iter]*vwater*vwater))/vwater;
      complex<double> * F2 = new complex<double>[exnt];
      complex<double> * denominator = new complex<double>[exnt];
      complex<double> * fdata = new complex<double>[exnt];
      
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        F2[j_iter] = 1.0 -exp(I*tau*exomega[j_iter]);
        denominator[j_iter] = 2.0 - 2.0*cos(tau*exomega[j_iter]) + recover_eps;
        fdata[j_iter] = ftaup_img[i_iter*exnt+j_iter];
      }
      
      
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        fdata[j_iter] = fdata[j_iter]*F2[j_iter]/denominator[j_iter];
        ftaup_img[i_iter*exnt+j_iter] = fdata[j_iter];
      }
      
      
      delete [] F2;
      delete [] denominator;
      delete [] fdata;
      
    }
    
  }
  
}


void TauP2DDeghost::taup_deghost(double* taup_img) {
  
  complex<double> * ftaup_img = new complex<double>[np*exnt]; 
  
  taup_deghost_tau_to_f(taup_img, ftaup_img);
  
  
  taup_deghost(ftaup_img, src_depth);
  taup_deghost(ftaup_img, rec_depth);
  
  
  taup_deghost_f_to_tau(ftaup_img, taup_img);
  
  delete [] ftaup_img;
  
}


void TauP2DDeghost::deghost() {
  
  
  int M = exnt*nx;
  double * exdata_double = new double[exnt*nx];
  double * taup_imag = new double[exnt*np];
  
  
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata_double[i_iter] = exdata[i_iter].real();
  }
  
  
//  biwrite("exdata_double", exdata_double, nx, exnt);
  
  
  // Forward TauP transform
  fwd_taup_transform(exnt, nx, np, dt, dx, dp, ft, fx, fp, exdata_double, taup_imag);
  
  
  biwrite("taup_imag_before_deghost", taup_imag, np, exnt);
  
  
  // TauPdeghost
  taup_deghost(taup_imag);
  
  
  biwrite("taup_imag_after_deghost", taup_imag, np, exnt);
  
  
  for (int i_iter = 0; i_iter < nx*exnt; i_iter++) {
    exdata_double[i_iter] = 0.0;
  }
  
  // Backward TauP transform
  bwd_taup_transform(exnt, nx, np, dt, dx, dp, ft, fx, fp, taup_imag, exdata_double);
  
  
//  biwrite("taup_inv", exdata_double, nx, exnt);
  
  
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata[i_iter] = exdata_double[i_iter];
  }
  
  
  delete [] exdata_double;
  delete [] taup_imag;
  
}


void TauP2DDeghost::execute() {
  
  // Phys filter
  phys_filtering();
  
  // Deghost
  deghost();
  
  // Phys filter
  phys_filtering();
  
  // Output
  output_data();
  
}


void TauP2DDeghost::finalize() {
  
  delete [] p;
  delete [] kx;
  delete [] exomega;
  delete [] phys_filter;
  delete [] exdata;
  
}





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
  f = NULL;
  omega = NULL;
  exkx = NULL;
  exf = NULL;
  exomega = NULL;
  phys_filter = NULL;
  exdata = NULL;
  exmissing_part = NULL;
  
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
  f = NULL;
  omega = NULL;
  exkx = NULL;
  exf = NULL;
  exomega = NULL;
  phys_filter = NULL;
  exdata = NULL;
  exmissing_part = NULL;
  
}


double TauP2DDeghost::maxabs(double* vec, int length) {
  
  double ret;
  
  ret = abs(vec[0]);
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    
    double test = abs(vec[i_iter]);
    
    if (test > ret) {
      ret = test;
    }
  }
  
  return ret;
  
}


double TauP2DDeghost::maxabs(complex<double>* vec, int length) {
  
  double ret;
  
  ret = abs(vec[0]);
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    
    double test = abs(vec[i_iter]);
    
    if (test > ret) {
      ret = test;
    }
  }
  
  return ret;
  
}


void TauP2DDeghost::get_input_info(TauP2DDeghost_info_t info) {
  
  vwater = info.vwater;
  recover_eps = info.recover_eps;
  src_depth = info.src_depth;
  rec_depth = info.rec_depth;
  
  nt = info.nt;
  nx = info.nx;
  np = info.np;
  npad = info.npad;
  dt = info.dt;
  dx = info.dx;
  fx = info.fx;
  
  nfb = 15;
  fmin = info.fmin;
  fmax = info.fmax;
  res = 5.0;
  
  p_control = 1.0;
  
  in = info.in;
  out = info.out;
  
}


void TauP2DDeghost::check_input_info() {
    
  if (recover_eps<=0) {
    cout << "recover_eps must be positive!!!" << endl;
    exit(-1);
  }
  
  if (src_depth<0) {
    cout << "src_depth must not be negative!!!" << endl;
    exit(-1);
  }
  
  if (rec_depth<0) {
    cout << "rec_depth must not be negative!!!" << endl;
    exit(-1);
  }
  
  if (nt<=0) {
    cout << "nt must be positive!!!" << endl;
    exit(-1);
  } 
  
  if (nx<=0) {
    cout << "nx must be positive!!!" << endl;
    exit(-1);
  }
  
  if (np<=0) {
    cout << "np must be positive!!!" << endl;
    exit(-1);
  }
  
  if (npad<0) {
    cout << "npad must not be negative!!!" << endl;
    exit(-1);
  }
  
  if (dt<=0) {
    cout << "dt must be positive!!!" << endl;
    exit(-1);
  }
  
  if (dx<=0) {
    cout << "dx must be positive!!!" << endl;
    exit(-1);
  }

  if (fmin<=0) {
    cout << "fmin must be positive!!!" << endl;
    exit(-1);
  }
  
  if (fmax<=0) {
    cout << "fmax must be positive!!!" << endl;
    exit(-1);
  }
  
  if (fmax<fmin) {
    cout << "fmax must be larger than fmin!!!" << endl;
    exit(-1);
  }
    
}


void TauP2DDeghost::calculate_dependent_parameters() {
  
  // N
  if (!(np%2)) {
    np = np+1;
  }
  nkx = nx;
  
  exnt = 2*nt;
  exnx = nx + 2*npad;
  exnkx = exnx;
  
  
  // D
  dp = 2.0/(vwater*(np-1));
  df = 1.0/(nt*dt);
  domega=2.0*MY_PI*df;
  dkx = 2.0*MY_PI/((nx*dx));
  
  exdt = dt;
  exdf = 1.0/(exnt*exdt);
  exdomega = 2.0*MY_PI*exdf;
  exdx = dx;
  exdkx = 2.0*MY_PI/(exnx*exdx);
  
  
  // F
  ft = 0.0;
  fp = -floor(np/2)*dp;
  fomega = 0.0;
  
  exft = ft;
  exfx = fx-npad*exdx;
  exfomega = fomega;
  
  
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
  cout << "exnx: " << exnx << endl;
  cout << "exnkx: " << exnkx << endl;
  cout << "dt: " << dt << endl;
  cout << "df: " << df << endl;
  cout << "domega: " << domega << endl;
  cout << "dx: " << dx << endl;
  cout << "dkx: " << dkx << endl;
  cout << "dp: " << dp << endl;
  cout << "exdt: " << exdt << endl;
  cout << "exdf: " << exdf << endl;
  cout << "exdomega: " << exdomega << endl;
  cout << "exdx: " << exdx << endl;
  cout << "exdkx: " << exdkx << endl;
  cout << "ft: " << ft << endl;
  cout << "fx: " << fx << endl;
  cout << "fp: " << fp << endl;
  cout << "fomega: " << fomega << endl;
  cout << "exft: " << exft << endl;
  cout << "exfx: " << exfx << endl;
  cout << "exfomega: " << exfomega << endl;
  cout << "pmax: " << pmax << endl;
  cout << "fmin: " << fmin << endl;
  cout << "fmax: " << fmax << endl;
  cout << "pcontrol: " << p_control << endl;
//  cout << "in: " << in << endl;
//  cout << "out: " << out << endl;
  cout << "*********************************************************************" << endl;
  
}


void TauP2DDeghost::allocate_main_mem() {
  
  p = new double[np];
  f = new double[nt];
  omega = new double[nt];
  exkx = new double[exnkx];
  exf = new double[exnt];
  exomega = new double[exnt];
  phys_filter = new double[exnt*exnkx];
  exdata = new complex<double>[exnx*exnt];
  exmissing_part = new complex<double>[exnx*exnt];
  
}


void TauP2DDeghost::compute_p() {
  
  for (int i_iter = 0; i_iter < np; i_iter++) {
    p[i_iter] = fp + i_iter*dp;
  }
  
}


void TauP2DDeghost::compute_omega() {
  
    
  for (int i_iter = 0; i_iter <= floor(nt/2); i_iter++) {
    omega[i_iter] = i_iter*domega;
    f[i_iter] = i_iter*df;
  }

  for (int i_iter = floor(nt/2) + 1; i_iter < nt; i_iter++) {
    omega[i_iter] = -omega[nt-i_iter];
    f[i_iter] = -f[nt-i_iter];
  }
 
  
}


void TauP2DDeghost::compute_exkx() {
  
  
  for (int i_iter = 0; i_iter <= floor(exnx/2); i_iter++) {
    exkx[i_iter] = i_iter*exdkx;
  }

  for (int i_iter = floor(exnx/2) + 1; i_iter < exnx; i_iter++) {
    exkx[i_iter] = -exkx[exnx-i_iter];
  }
  
  
}


void TauP2DDeghost::compute_exomega() {
  
    
  for (int i_iter = 0; i_iter <= floor(exnt/2); i_iter++) {
    exomega[i_iter] = i_iter*exdomega;
    exf[i_iter] = i_iter*exdf;
  }

  for (int i_iter = floor(exnt/2) + 1; i_iter < exnt; i_iter++) {
    exomega[i_iter] = -exomega[exnt-i_iter];
    exf[i_iter] = -exf[exnt-i_iter];
  }
 
  
}


void TauP2DDeghost::compute_phys_filter() {
  
  for (int i_iter = 0; i_iter < exnt; i_iter++) {
    for (int j_iter = 0; j_iter < exnkx; j_iter++) {
      double tmp = (exomega[i_iter]*exomega[i_iter])/(vwater*vwater) - exkx[j_iter]*exkx[j_iter];
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
  compute_omega();
  compute_exkx();
  compute_exomega();
  compute_phys_filter();
  
}


void TauP2DDeghost::get_exdata() {
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      exdata[(npad+i_iter)*exnt+j_iter] = in[i_iter*nt+j_iter];
    }
  }
  
}


void TauP2DDeghost::get_output_data() {
  
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nt+j_iter] = exdata[(npad+i_iter)*exnt+j_iter].real();
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
  
  
  // Check input info (Can be commented)
  check_input_info();
  
  
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
//  biwrite_complex("exinput", exdata, exnx, exnt);

  
}


void TauP2DDeghost::phys_filtering() {
  
  int Nh = floor(exnt/2) + 1;
  int M = exnt*exnx;
  complex<double> * exdata_fk = new complex<double>[exnt*exnx];
  
  uniform_fftwf2d(exnt, exnx, exdata, exdata_fk);
  
  
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata_fk[i_iter] = exdata_fk[i_iter]*phys_filter[i_iter];
  }
  
  
  // Keep symmetry
  for (int i_iter = 1; i_iter< Nh; i_iter++) { 
    exdata_fk[exnt-i_iter] = conj(exdata_fk[i_iter]);
  }


  for (int i_iter = 1; i_iter < Nh; i_iter++) {
    for (int j_iter = 1; j_iter < exnkx; j_iter++) {
      exdata_fk[(exnkx-j_iter)*exnt+exnt-i_iter] = conj(exdata_fk[j_iter*exnt+i_iter]);
    }
  }
  
  uniform_ifftwf2d(exnt, exnx, exdata_fk, exdata);
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
  TxTaup2D taup_engine;
  
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
  info.fmax = fmax; 
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
  TxTaup2D taup_engine;
  
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
  info.fmax = fmax;
  info.mode = INV_TAUP;
  info.in = in;
  info.out = out;
  
  
  taup_engine.initialize(info);
  taup_engine.execute();
  taup_engine.finalize();
  
}


void TauP2DDeghost::taup_deghost_tau_to_f(double* taup_img, complex<double>* ftaup_img) {
  
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    double test = p_control - (p[i_iter]*p[i_iter]*vwater*vwater);
    
    if (test > 0) {
      
      fftw_complex * tmp_in;
      fftw_complex * tmp_out;

      tmp_in = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
      tmp_out = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));

      // Take one trace
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        tmp_in[j_iter][0] = taup_img[i_iter*exnt+j_iter];
        tmp_in[j_iter][1] = 0.0;
      }


      uniform_fftwf1d(exnt, tmp_in, tmp_out, FFTW_FORWARD);


      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        ftaup_img[i_iter*exnt + j_iter] = complex<double>(tmp_out[j_iter][0],tmp_out[j_iter][1]);
      }
      
      fftw_free(tmp_in);
      fftw_free(tmp_out);
    
    } else {
      
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        taup_img[i_iter*exnt+j_iter] = 0.0;
        ftaup_img[i_iter*exnt + j_iter] = 0.0;
      }
      
    }
    
  }

  
}


void TauP2DDeghost::taup_deghost_f_to_tau(complex<double>* ftaup_img, double* taup_img) {
  
  int Nh = floor(exnt/2)+1;
  
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    double test = p_control - (p[i_iter]*p[i_iter]*vwater*vwater);
    
    if (test > 0) {
      
      fftw_complex * tmp_in;
      fftw_complex * tmp_out;

      tmp_in = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
      tmp_out = (fftw_complex *) fftw_malloc(exnt*sizeof(fftw_complex));
    
      for (int j_iter = 1; j_iter < Nh; j_iter++) {
        ftaup_img[i_iter*exnt+exnt-j_iter] = conj(ftaup_img[i_iter*exnt+j_iter]);
      }

      // Take one trace
      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        tmp_in[j_iter][0] = ftaup_img[i_iter*exnt+j_iter].real();
        tmp_in[j_iter][1] = ftaup_img[i_iter*exnt+j_iter].imag();
      }


      uniform_fftwf1d(exnt, tmp_in, tmp_out, FFTW_BACKWARD);


      for (int j_iter = 0; j_iter < exnt; j_iter++) {
        taup_img[i_iter*exnt + j_iter] = tmp_out[j_iter][0];
      }
      
      
      fftw_free(tmp_in);
      fftw_free(tmp_out);
      
    
    }
    
  }
  
  
}


void TauP2DDeghost::taup_deghost(complex<double> * ftaup_img, double depth) {
  
#pragma omp parallel for
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


void TauP2DDeghost::taup_deghost(double* taup_img_in, double* taup_img_out) {
  
  complex<double> * ftaup_img = new complex<double>[np*exnt]; 
  
  taup_deghost_tau_to_f(taup_img_in, ftaup_img);
  
  
  taup_deghost(ftaup_img, src_depth);
  taup_deghost(ftaup_img, rec_depth);
  
  
  taup_deghost_f_to_tau(ftaup_img, taup_img_out);
  
  delete [] ftaup_img;
  
}


void TauP2DDeghost::calculate_freq_hz(double* freq_hz) {
  
  double fratio = fmin/fmax;
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    double ww = double(i_iter)/double(nfb-1);
    freq_hz[i_iter] = fmax*pow(fratio,ww);
  }
  
}


void TauP2DDeghost::fourier_img(double* input, complex<double>* finput, int nt) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
   
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    tmp_in[i_iter][0] = input[i_iter];
    tmp_in[i_iter][1] = 0.0;
  }


  uniform_fftwf1d(nt, tmp_in, tmp_out, FFTW_FORWARD);


  for (int i_iter = 0; i_iter < nt; i_iter++) {
    finput[i_iter] = complex<double>(tmp_out[i_iter][0],tmp_out[i_iter][1]);
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TauP2DDeghost::inv_fourier(complex<double>* finput, double* output, int nt) {
  
  fftw_complex * tmp_in;
  fftw_complex * tmp_out;
  
  tmp_in = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
  tmp_out = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));
   
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    tmp_in[i_iter][0] = finput[i_iter].real();
    tmp_in[i_iter][1] = finput[i_iter].imag();
  }


  uniform_fftwf1d(nt, tmp_in, tmp_out, FFTW_BACKWARD);


  for (int i_iter = 0; i_iter < nt; i_iter++) {
    output[i_iter] = tmp_out[i_iter][0];
  }
  
  
  fftw_free(tmp_in);
  fftw_free(tmp_out);
  
}


void TauP2DDeghost::calculate_wavelet(double* freq_hz, double* wavelet) {
  
  int M = nt*nfb;
  int Nh = floor(nt/2)+1;
  double * sum_wavelet = new double[nt];
//  double * test_sum = new double[nt];
 
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    sum_wavelet[i_iter] = 0.0;
  }
  
  for (int i_iter = 0; i_iter < M; i_iter++) {
    wavelet[i_iter] = 0.0;
  }
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    for (int j_iter = 0; j_iter < Nh; j_iter++) {
      double rarg = (1.0-f[j_iter]/freq_hz[i_iter])*res;
      if (rarg > 8) {
        rarg = 8; 
      } 
      wavelet[i_iter*nt+j_iter] = exp(-0.5*rarg*rarg);
      sum_wavelet[j_iter] += wavelet[i_iter*nt+j_iter];
    }
  }
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    for (int j_iter = 0; j_iter < Nh; j_iter++) {
      wavelet[i_iter*nt+j_iter] = wavelet[i_iter*nt+j_iter]/sum_wavelet[j_iter];
    }
  }
  
    
//  for (int i_iter = 0; i_iter < nfb; i_iter++) {
//    for (int j_iter = 0; j_iter < Nh; j_iter++) {
//      test_sum[j_iter] += wavelet[i_iter*nt+j_iter];
//    }
//  }
  
  
//  delete [] test_sum;
  delete [] sum_wavelet;
  
}


void TauP2DDeghost::excalculate_wavelet(double * freq_hz, double* wavelet) {
  
  int M = exnt*nfb;
  int Nh = floor(exnt/2)+1;
  double * sum_wavelet = new double[exnt];
//  double * test_sum = new double[exnt];
 
  
  for (int i_iter = 0; i_iter < exnt; i_iter++) {
    sum_wavelet[i_iter] = 0.0;
  }
  
  for (int i_iter = 0; i_iter < M; i_iter++) {
    wavelet[i_iter] = 0.0;
  }
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    for (int j_iter = 0; j_iter < Nh; j_iter++) {
      double rarg = (1.0-exf[j_iter]/freq_hz[i_iter])*res;
      if (rarg > 8) {
        rarg = 8; 
      } 
      wavelet[i_iter*exnt+j_iter] = exp(-0.5*rarg*rarg);
      sum_wavelet[j_iter] += wavelet[i_iter*exnt+j_iter];
    }
  }
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    for (int j_iter = 0; j_iter < Nh; j_iter++) {
      wavelet[i_iter*exnt+j_iter] = wavelet[i_iter*exnt+j_iter]/sum_wavelet[j_iter];
    }
  }
  
    
//  for (int i_iter = 0; i_iter < nfb; i_iter++) {
//    for (int j_iter = 0; j_iter < Nh; j_iter++) {
//      test_sum[j_iter] += wavelet[i_iter*exnt+j_iter];
//    }
//  }
  
  
//  for (int i_iter = 0; i_iter < Nh; i_iter++) {
//    cout << "sum_wavelet[" << i_iter << "]: " << sum_wavelet[i_iter] << endl;
//  }
  
  
//  delete [] test_sum;
  delete [] sum_wavelet;
  
}


void TauP2DDeghost::fanalyse_spec_fwd(double* input, complex<double> * fspecout) {
  
  int M = nfb*nt;
  int Nh = floor(nt/2)+1;
  double * freq_hz = new double[nfb];
  complex<double> * finput = new complex<double>[nt];
  double * wavelet = new double[nfb*nt];
  
  
  // Set zero output
  for (int i_iter = 0; i_iter < M; i_iter++) {
    fspecout[i_iter] = 0;
  }
  
  
  fourier_img(input, finput, nt);
  calculate_freq_hz(freq_hz);
  calculate_wavelet(freq_hz, wavelet);
  
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    
    complex<double> * out_wavelet = new complex<double>[nt];
    
    for (int j_iter = 0; j_iter < Nh; j_iter++) {
      out_wavelet[j_iter] = finput[j_iter]*wavelet[i_iter*nt+j_iter];
    }
    
    for (int j_iter = Nh; j_iter < nt; j_iter++) {
      out_wavelet[j_iter] = conj(out_wavelet[nt-j_iter]);
    }
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      fspecout[i_iter*nt+j_iter]=out_wavelet[j_iter];
    }
    
    delete [] out_wavelet;
    
  } 
  
  
  delete [] wavelet;
  delete [] finput;
  delete [] freq_hz;
  
}


void TauP2DDeghost::extanalyse_spec_fwd(double* input, double * tspecout) {
  
  int M = nfb*exnt;
  int Nh = floor(exnt/2)+1;
  double * freq_hz = new double[nfb];
  complex<double> * finput = new complex<double>[exnt];
  double * exwavelet = new double[nfb*exnt];

  
  // Set zero output
  for (int i_iter = 0; i_iter < M; i_iter++) {
    tspecout[i_iter] = 0;
  }
  
  fourier_img(input, finput, exnt);
  calculate_freq_hz(freq_hz);
  excalculate_wavelet(freq_hz, exwavelet);
  
  
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    
    complex<double> * out_wavelet = new complex<double>[exnt];
    
    for (int j_iter = 0; j_iter < Nh; j_iter++) {
      out_wavelet[j_iter] = finput[j_iter]*exwavelet[i_iter*exnt+j_iter];
    }
    
    for (int j_iter = Nh; j_iter < exnt; j_iter++) {
      out_wavelet[j_iter] = conj(out_wavelet[exnt-j_iter]);
    }
    
    inv_fourier(out_wavelet, &tspecout[i_iter*exnt], exnt);
    
    delete [] out_wavelet;
    
  } 
  
  delete [] exwavelet;
  delete [] finput;
  delete [] freq_hz;
  
}


void TauP2DDeghost::cleanup(double* taup_img_in, double* taup_img_out) {
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < np; i_iter++) {
    
    double * tmp_in = new double[exnt];
    double * tmp_out = new double[exnt];
    double * tspecin = new double[exnt*nfb];
    double * tspecout = new double[exnt*nfb];
    
    
    for (int j_iter = 0; j_iter < exnt; j_iter++) {
      tmp_in[j_iter] = taup_img_in[i_iter*exnt+j_iter];
      tmp_out[j_iter] = taup_img_out[i_iter*exnt+j_iter];
    }
    
    
    extanalyse_spec_fwd(tmp_in, tspecin);
    extanalyse_spec_fwd(tmp_out, tspecout);
    
    
    for (int j_iter = 0; j_iter < exnt; j_iter++) {
      
      double sum1 = 0;
      double sum2 = 0;
      
      for (int k_iter = 0; k_iter < nfb; k_iter++) {
        sum1 += tspecin[k_iter*exnt+j_iter];
        sum2 += tspecout[k_iter*exnt+j_iter];
      }
      
      sum1 = abs(sum1);
      sum2 = abs(sum2);
      
      if (sum2>1.0*sum1) {
        for (int k_iter = 0; k_iter < nfb; k_iter++) {
          tspecout[k_iter*exnt+j_iter] = 0;
        }
      }
      
    }
    
    
    for (int j_iter = 0; j_iter < exnt; j_iter++) {
      
      taup_img_out[i_iter*exnt+j_iter] = 0.0;
      
      for (int k_iter = 0; k_iter < nfb; k_iter++) {
        taup_img_out[i_iter*exnt+j_iter] += tspecout[k_iter*exnt+j_iter];
      }
      
    }
    
    
    delete [] tspecin;
    delete [] tspecout;
    delete [] tmp_in;
    delete [] tmp_out;
    
  }
  
}


void TauP2DDeghost::calculate_exmissing_part(double* exinput, double* exinv) {

  int M = exnt*exnx;
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exmissing_part[i_iter] = exinput[i_iter] - exinv[i_iter];
  }
  
}


void TauP2DDeghost::deghost() {
  
  int M = exnt*exnx;
  double * exdata_double = new double[exnt*exnx];
  double * exdata_inv_double = new double[exnt*exnx];
  double * taup_img_in = new double[exnt*np];
  double * taup_img_out = new double[exnt*np];
  
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata_double[i_iter] = exdata[i_iter].real();
  }
  
  
//  biwrite("exdata_double", exdata_double, exnx, exnt);
  
//  cout << "Forward TauP transform starts." << endl;
  
  // Forward TauP transform
  fwd_taup_transform(exnt, exnx, np, exdt, exdx, dp, exft, exfx, fp, exdata_double, taup_img_in);
  
//  cout << "Forward TauP transform ends." << endl;
  
//  biwrite("taup_imag_before_deghost", taup_img_in, np, exnt);
  
  
  // Backward TauP transform
  bwd_taup_transform(exnt, exnx, np, exdt, exdx, dp, exft, exfx, fp, taup_img_in, exdata_inv_double);
  
  
  // Calculate exmissing part
  calculate_exmissing_part(exdata_double, exdata_inv_double);
  
  
//  cout << "TauP Deghost starts." << endl;
  
  // TauPdeghost
  taup_deghost(taup_img_in, taup_img_out);
  
//  cout << "TauP Deghost ends." << endl;
  
//  biwrite("taup_imag_after_deghost", taup_img_out, np, exnt);
  
//  cout << "Cleanup starts." << endl;
  
  // Clean taup output
//  cleanup(taup_img_in, taup_img_out);
  
//  cout << "Cleanup ends." << endl;
  
//  biwrite("taup_imag_after_deghost_clean", taup_img_out, np, exnt);
  
//  cout << "Backward TauP starts." << endl;
  
  // Backward TauP transform
  bwd_taup_transform(exnt, exnx, np, exdt, exdx, dp, exft, exfx, fp, taup_img_out, exdata_double);
  
//  cout << "Backward TauP ends." << endl;
  
//  biwrite("taup_inv", exdata_double, exnx, exnt);
  
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata[i_iter] = exdata_double[i_iter];
  }
  
  
  delete [] exdata_double;
  delete [] exdata_inv_double;
  delete [] taup_img_in;
  delete [] taup_img_out;
  
}


void TauP2DDeghost::reshape_spectrum() {
  
  double * sum_in = new double[nt];
  double * sum_out = new double[nt];
  complex<double> * fspecin = new complex<double>[nt*nfb];
  complex<double> * fspecout = new complex<double>[nt*nfb];
  double * ratio = new double[nfb];
  double * freq_hz = new double[nfb];
  
  
  // Calculate freq_hz
  calculate_freq_hz(freq_hz);
  
  
  // Set zero sum_in and sum_out
#pragma omp parallel for
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    sum_in[i_iter] = 0;
    sum_out[i_iter] = 0;
  }
  
  
  // sum input traces and output traces
#pragma omp parallel for
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    for (int j_iter  = 0; j_iter < nx; j_iter++) {
      sum_in[i_iter] += in[j_iter*nt+i_iter];
      sum_out[i_iter] += out[j_iter*nt+i_iter];
    }
  }
  
  
  // Analyze spec
  fanalyse_spec_fwd(sum_in, fspecin);
  fanalyse_spec_fwd(sum_out, fspecout);
  
  
  // Calculate ratio
#pragma omp parallel for
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    double maxin = maxabs(&fspecin[i_iter*nt],nt);
    double maxout = maxabs(&fspecout[i_iter*nt],nt);
    if (maxout>0) {
      ratio[i_iter] = maxin/maxout;
    } else {
      ratio[i_iter] = 0;
    }
  }
  
  
  // Keep away from low freq part
  for (int i_iter = 0; i_iter < nfb; i_iter++) {
    if ((freq_hz[i_iter]<80)&&(ratio[i_iter]<1.0)) {
      ratio[i_iter]=1.0;
    }
  }
  
  
// Reshape the spec by multiplying by the ratio
#pragma omp parallel for
  for (int i_iter = 0; i_iter < nx; i_iter++) {
    
    double * tmp_out = new double[nt];
    complex<double> * inner_fspecout = new complex<double>[nt*nfb];
    
    
    // Take one trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      tmp_out[j_iter] = out[i_iter*nt+j_iter];
    }
    
    
    // Analyze spec
    fanalyse_spec_fwd(tmp_out,inner_fspecout);
    
    
    // Reshape
    for (int j_iter = 0; j_iter < nfb; j_iter++) {
      for (int k_iter = 0; k_iter < nt; k_iter++) {
        inner_fspecout[j_iter*nt+k_iter]*=ratio[j_iter];
      }
    }
    
    
    // Transform back to time domain
    for (int j_iter = 0; j_iter < nfb; j_iter++) {
      inv_fourier(&inner_fspecout[j_iter*nt], tmp_out, nt);
      for (int k_iter = 0; k_iter < nt; k_iter++) {
        inner_fspecout[j_iter*nt+k_iter] = tmp_out[k_iter];
      }
    }
    
    
    // Set zero this trace
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      out[i_iter*nt+j_iter]=0;  
    }
    
    
    for (int j_iter = 0; j_iter < nt; j_iter++) {
      for (int k_iter = 0; k_iter < nfb; k_iter++) {
        out[i_iter*nt+j_iter]+=real(inner_fspecout[k_iter*nt+j_iter]);
      }
    }
    
    
    delete [] inner_fspecout;
    delete [] tmp_out;
    
  }
  
  
  delete [] ratio;
  delete [] fspecin;
  delete [] fspecout;
  delete [] sum_in;
  delete [] sum_out;
  delete [] freq_hz;
  
}


void TauP2DDeghost::scale_back_for_fix_display() {
  
  int M = nx*nt;
  int exM = exnx*exnt;
  double maxin=maxabs(in,M);
  double maxout=maxabs(exdata,exM);
  double ratio;
  
  
  if (maxin>0) {
    ratio = maxout/maxin;
  } else {
    ratio = 0;
  }
  
  
  if (ratio == 0) {
    cout << "max value of input equals zero" << endl;
    exit(-1);
  }
  
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < exM; i_iter++) {
    exdata[i_iter]/=ratio;
  }
  
}

void TauP2DDeghost::add_back_missing_part() {

  int M = exnx*exnt;
  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < M; i_iter++) {
    exdata[i_iter]+=exmissing_part[i_iter];
  }
  
}


void TauP2DDeghost::execute() {
  
  // Deghost
  deghost();
  
  // Phys filter
  phys_filtering();
  
  // Scale back for fix display
  scale_back_for_fix_display();
  
  // Add back missing part
//  add_back_missing_part();
  
  // get Output
  get_output_data();
  
  // Reshape the spectrum
//  reshape_spectrum();
  
}


void TauP2DDeghost::finalize() {
  
  delete [] p;
  delete [] f;
  delete [] omega;
  delete [] exkx;
  delete [] exf;
  delete [] exomega;
  delete [] phys_filter;
  delete [] exdata;
  delete [] exmissing_part;
  
}





#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <omp.h>
#include "TauP2DDeghost_info.H"
#include "TauP2DDeghost.H"


using namespace std;



void readbinfile
(
  string filename, 
  int M_input, 
  int Nt, 
  double * data
) 
{
  
  ifstream myfile; 
  int num_elems = M_input*Nt;
  int data_size = num_elems*sizeof(float);
  float * memblock = new float[num_elems];
  
  
  myfile.open( filename.c_str(), ios::binary | ios::in );
  
  
  if (myfile.is_open()) {
    
    myfile.read((char*) memblock, data_size);
    myfile.close();
    
    cout << "*********************************************************************" << endl;
    cout << "readbinfile()" << endl;
    cout << "*********************************************************************" << endl;
    cout << "The complete file content is in memory." << std::endl;
    cout << "*********************************************************************" << endl;
    cout << endl;
    cout << endl;
    
  } else {
    std::cout << "Unable to open file" << filename << "!!!" << endl;
    exit(-1);
  }
  
  
//#pragma omp parallel for
  for (int i_iter = 0; i_iter < num_elems; i_iter++) {
    data[i_iter] = memblock[i_iter];
  }
  
      
  delete [] memblock;
  
}


void writebinfile
(
  string filename, 
  int M_input, 
  int Nt, 
  double * data
) 
{
  
  ofstream myfile; 
  int num_elems = M_input*Nt;
  int data_size = num_elems*sizeof(float);
  float * memblock = new float[num_elems];
  
  
//#pragma omp parallel for
  for (int i_iter = 0; i_iter < num_elems; i_iter++) {
    memblock[i_iter] = data[i_iter];
  }
  
  
  myfile.open( filename.c_str(), ios::binary | ios::out );
  
  
  if (myfile.is_open()) {
    
    myfile.write((char*) memblock, data_size);
    myfile.close();
    
    cout << "*********************************************************************" << endl;
    cout << "writebinfile()" << endl;
    cout << "*********************************************************************" << endl;
    cout << "The complete mem is written to file." << std::endl;
    cout << "*********************************************************************" << endl;
    cout << endl;
    cout << endl;
    
  } else {
    std::cout << "Unable to open file" << filename << "!!!" << endl;
    exit(-1);
  }
  
  
  delete [] memblock;
  
}


void get_basic_info
(
  string parm_file_name, 
  string & input_file_name, 
  string & deghost_file_name,
  TauP2DDeghost_info_t & info
) 
{

  ifstream myfile;
  string tmp_string;

  myfile.open(parm_file_name.c_str(), ios::in);


  if (myfile.is_open()) {

    myfile >> tmp_string >> input_file_name;
    myfile >> tmp_string >> deghost_file_name;
    myfile >> tmp_string >> info.vwater;
    myfile >> tmp_string >> info.recover_eps;
    myfile >> tmp_string >> info.src_depth;
    myfile >> tmp_string >> info.rec_depth;
    myfile >> tmp_string >> info.nt;
    myfile >> tmp_string >> info.nx;
    myfile >> tmp_string >> info.np;
    myfile >> tmp_string >> info.npad;
    myfile >> tmp_string >> info.dt;
    myfile >> tmp_string >> info.dx;
    myfile >> tmp_string >> info.fx;
    myfile >> tmp_string >> info.fmin;
    myfile >> tmp_string >> info.fmax;

    myfile.close();
    
  } else {
    cout << "Unable to open file: " << parm_file_name << endl;
    exit(-1);
  }

    
  // Display to check
  cout << "*********************************************************************"<<endl;
  cout << "****************** Basic info ***************************************"<<endl;
  cout << "input file: " << input_file_name << endl;
  cout << "deghost file: " << deghost_file_name << endl;
  cout << "vwater: " << info.vwater << endl;
  cout << "recover_eps: " << info.recover_eps << endl;
  cout << "src_depth: " << info.src_depth << endl;
  cout << "rec_depth: " << info.rec_depth << endl;
  cout << "nt: " << info.nt << endl;
  cout << "nx: " << info.nx << endl;
  cout << "np: " << info.np << endl;
  cout << "npad: " << info.npad << endl;
  cout << "dt: " << info.dt << endl;
  cout << "dx: " << info.dx << endl;
  cout << "fx: " << info.fx << endl;
  cout << "fmin: " << info.fmin << endl;
  cout << "fmax: " << info.fmax << endl;
  cout << "*********************************************************************"<<endl;
  cout << endl;
  
}


int
main() {
  
  // Variable for measuring time consuming
  double start, end, timediff;
  
  string input_parm_file = "input.inf";
  string input_file_name;
  string deghost_file_name;
  
  
  TauP2DDeghost_info_t deghost_info;
  TauP2DDeghost solver;

  
  double * input_trace=  NULL;
  double * deghost_trace = NULL;
  
  
  get_basic_info(input_parm_file, input_file_name, deghost_file_name, deghost_info);
  
  
  input_trace = new double [deghost_info.nt*deghost_info.nx];
  deghost_trace = new double [deghost_info.nt*deghost_info.nx];
  
  
  deghost_info.in = input_trace;
  deghost_info.out = deghost_trace;
  
  
  readbinfile(input_file_name, deghost_info.nx, deghost_info.nt, input_trace);
  
    
  start = omp_get_wtime();
  

  solver.initialize(deghost_info);
  solver.execute();
  solver.finalize();
  
  
  end = omp_get_wtime();
  
  
  timediff = end - start;
  printf("TauP deghost : %lf s\n\n",timediff);  
  
  
  writebinfile(deghost_file_name, deghost_info.nx, deghost_info.nt, deghost_trace);
  
  
  delete [] input_trace;
  delete [] deghost_trace;
  
  return 0;

}




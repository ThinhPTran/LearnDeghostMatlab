#include <iostream>
#include <fstream>
#include <string>
#include "deghost_basic_info.H"
#include "General_deghost_engine.H"
#include "ghost_data_block.H"
#include "window_deghost_process_basic_info.H"
#include "window_process_deghosting.H"
#include <omp.h>


using namespace std;



void readbinfile
(
  string filename, 
  int M_input, 
  int Nt, 
  float * data
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
  
  
#pragma omp parallel for
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
  float * data
) 
{
  
  ofstream myfile; 
  int num_elems = M_input*Nt;
  int data_size = num_elems*sizeof(float);
  
  
  myfile.open( filename.c_str(), ios::binary | ios::out );
  
  
  if (myfile.is_open()) {
    
    myfile.write((char*) data, data_size);
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
  
}


template<class T> 
void print_vec
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


int main () 
{
  
  string infilename("flat_ghost.bin");
  string outfilename("flat_ghost_out.bin");
  int nt = 3001;
  int n2 = 100;
  float dt = 0.000666;
  int ntwindow = 75;
  int wshift = 1;
  float zmin = 30;
  float zmax = 30;
  float vwater = 5000;
  float epsilon = 0.00000001;
  float addback = 0.1;
  window_deghost_process_basic_info window_process_basic_info;
  float * full_data = new float[n2*nt];
  float * traces_out = new float[n2*nt];
  
  
  // Read all data into mem
  readbinfile(infilename,n2,nt,full_data);
  
  
  window_process_basic_info.set_nt(nt);
  window_process_basic_info.set_dt(dt);
  window_process_basic_info.set_zmax(zmax);
  window_process_basic_info.set_zmin(zmin);
  window_process_basic_info.set_vwater(vwater);
  window_process_basic_info.set_epsilon(epsilon);
  window_process_basic_info.set_addback(addback);
  window_process_basic_info.set_ntwindow(ntwindow);
  window_process_basic_info.set_wshift(wshift);

  
#pragma omp parallel for
  for (int i_iter = 0; i_iter < n2; i_iter++) {
    window_process_deghosting solver;
    solver.set_basic_info(window_process_basic_info);
    solver.set_input_trace(&full_data[i_iter*nt],nt);
    solver.execute();
    solver.get_output_trace(&traces_out[i_iter*nt],nt);
    solver.finalize();
  }
  
  
  writebinfile(outfilename,n2,nt,traces_out);
  
  
  delete [] traces_out;
  delete [] full_data;
  
}

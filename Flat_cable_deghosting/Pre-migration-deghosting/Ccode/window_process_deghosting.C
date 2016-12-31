/* 
 * File:   window_process_deghosting.C
 * Author: thinh
 * 
 * Created on September 20, 2013, 10:20 AM
 */

#include "window_process_deghosting.H"


window_process_deghosting::window_process_deghosting() {
  
  first_run = 1;
  realloc_flag = 0;
  array_of_blocks = NULL;
  trace_in = NULL;
  trace_out = NULL;
  weight = NULL;
  
}


window_process_deghosting::~window_process_deghosting() {
  
  first_run = 0;
  realloc_flag = 0;
  
  
  if (array_of_blocks!=NULL) {
    delete [] array_of_blocks;
    array_of_blocks = NULL;
  }
  
  if (trace_in!=NULL) {
    delete [] trace_in;
    trace_in = NULL;
  }
  
  if (trace_out!=NULL) {
    delete [] trace_out;
    trace_out = NULL;
  }
  
  if (weight!=NULL) {
    delete [] weight;
    weight = NULL;
  }
  
  
}


template<class T> 
void window_process_deghosting::print_vec
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


void window_process_deghosting::display_info() {
  
  cout << "********************************************************************" << endl;
  cout << "           window_process_deghosting::display_info()                " << endl;
  cout << "********************************************************************" << endl;
  cout << "Given variables " << endl;
  cout << "first_run: " << first_run << endl;
  cout << "realloc_flag: " << realloc_flag << endl;
  cout << "nt: " << basic_info.get_nt() << endl;
  cout << "dt: " << basic_info.get_dt() << endl;
  cout << "zmin: " << basic_info.get_zmin() << endl;
  cout << "zmax: " << basic_info.get_zmax() << endl;
  cout << "vwater: " << basic_info.get_vwater() << endl;
  cout << "ntwindow: " << basic_info.get_ntwindow() << endl;
  cout << "wshift: " << basic_info.get_wshift() << endl;
  cout << "nw: " << nw << endl;
  cout << "epsilon: " << basic_info.get_epsilon() << endl;
  cout << "addback: " << basic_info.get_addback() << endl;
  cout << "********************************************************************" << endl;
  
}


void window_process_deghosting::prepare_mem() {
  
  int nt = basic_info.get_nt();
  
  if (trace_in!=NULL) {
    delete [] trace_in;
    trace_in = NULL;
  }
  trace_in = new float[nt];
  memset(trace_in, 0x00, nt*sizeof(float));
  
  
  if (trace_out!=NULL) {
    delete [] trace_out;
    trace_out = NULL;
  }
  trace_out = new float[nt];
  memset(trace_out, 0x00, nt*sizeof(float));
  
  
  if (weight!=NULL) {
    delete [] weight;
    weight=NULL;
  }
  weight = new float[nt];
  memset(weight, 0x00, nt*sizeof(float));
  
  
  if (array_of_blocks!=NULL) {
    delete [] array_of_blocks;
    array_of_blocks = NULL;
  }
  array_of_blocks = new ghost_data_block[nw];
  
}


void window_process_deghosting::calculate_basic_info() {
  
  int nt = basic_info.get_nt();
  int ntwindow = basic_info.get_ntwindow();
  int wshift = basic_info.get_wshift();
  
  nw = floor((nt - ntwindow)/wshift) + 1;
  
}


void window_process_deghosting::set_basic_info(window_deghost_process_basic_info & info) {
  
  int old_nt = basic_info.get_nt();
  int old_ntwindow = basic_info.get_ntwindow();
  int old_wshift = basic_info.get_wshift();
  int new_nt = 0;
  int new_ntwindow = 0;
  int new_wshift = 0;
  
  
  basic_info = info;
  new_nt = basic_info.get_nt();
  new_ntwindow = basic_info.get_ntwindow();
  new_wshift = basic_info.get_wshift();
  
  
  calculate_basic_info();
  
  
  if (
      (old_nt!=new_nt)
    ||(old_ntwindow!=new_ntwindow)
    ||(old_wshift!=new_wshift)
     ) 
  {
    realloc_flag = 1;
  }
  
  
  if (first_run||realloc_flag) {
    prepare_mem();
  }
  
  
  display_info(); 
  
}


void window_process_deghosting::blocking() {
  
//  cout << "call blocking!!!" << endl;
  
  
  int nt = basic_info.get_nt();
  int wshift = basic_info.get_wshift();
  int ntwindow = basic_info.get_ntwindow();
  
  for (int i_iter = 0; i_iter < nw-1; i_iter++) {
    array_of_blocks[i_iter].set_basic_info(basic_info);
    array_of_blocks[i_iter].set_fidx(i_iter*wshift);
    array_of_blocks[i_iter].set_lidx(i_iter*wshift+ntwindow-1);
  }
  array_of_blocks[nw-1].set_basic_info(basic_info);
  array_of_blocks[nw-1].set_fidx((nw-1)*wshift);
  array_of_blocks[nw-1].set_lidx(nt-1);
  
  
  for (int i_iter = 0; i_iter < nw; i_iter++) {
    int fidx = array_of_blocks[i_iter].get_fidx();
    int lidx = array_of_blocks[i_iter].get_lidx();
    int length = lidx - fidx + 1;
//    cout << "Block: " << i_iter << endl;
//    cout << "   fidx: " << fidx << endl;
//    cout << "   lidx: " << lidx << endl;
    array_of_blocks[i_iter].set_input_trace(&trace_in[fidx],length);
  }
  
  
}


void window_process_deghosting::inversed_blocking() {
  
  int nt = basic_info.get_nt();
  float * tmpfloat = new float[nt];
  memset(trace_out,0x00,nt*sizeof(float));
  
  
  for (int i_iter = 0; i_iter < nw; i_iter++) {
    int fidx = array_of_blocks[i_iter].get_fidx();
    int lidx = array_of_blocks[i_iter].get_lidx();
    int length = lidx - fidx + 1;
    array_of_blocks[i_iter].get_output_trace(tmpfloat,length);
    for (int j_iter = 0; j_iter < length; j_iter++) {
      weight[fidx + j_iter] += 1.0;
      trace_out[fidx + j_iter] += tmpfloat[j_iter];
    }
  }
  
  
  for (int i_iter = 0; i_iter < nt; i_iter++) {
    trace_out[i_iter] /= weight[i_iter];
  }
  
  
  memset(weight, 0x00, nt*sizeof(float));
  
    
  delete [] tmpfloat;
    
}


void window_process_deghosting::execute() {
  
  int nt = basic_info.get_nt();
  
  blocking();
  
//  cout << "Begin to call engine.process!!!" << endl;
  
  
  for (int i_iter = 0; i_iter < nw; i_iter++) {
//    cout << "&array_of_blocks[i_iter]: " << &array_of_blocks[i_iter] << endl;
    engine.process(&array_of_blocks[i_iter]);
  }
  
   
  inversed_blocking();
  
//  print_vec(trace_out,nt);
  
}


void window_process_deghosting::finalize() {
  
  first_run = 0;
  realloc_flag = 0;
  
}


void window_process_deghosting::set_input_trace(float* input_trace, int length) {
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    trace_in[i_iter] = input_trace[i_iter];
  }
  
}


void window_process_deghosting::get_output_trace(float* output_trace, int length) {
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    output_trace[i_iter] = trace_out[i_iter];
  }  
  
}






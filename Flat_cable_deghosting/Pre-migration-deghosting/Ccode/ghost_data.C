/* 
 * File:   ghost_data.C
 * Author: thinh
 * 
 * Created on September 19, 2013, 1:15 PM
 */

#include "ghost_data.H"


ghost_data::ghost_data() {
  
  nt  =  0;
  dt = 0;
  zmax = 0;
  zmin = 0;
  vwater = 0;
  
  trace_in = NULL;
  trace_out = NULL;  
  
}


ghost_data::~ghost_data() {
  
  if (trace_in!=NULL) {
    delete [] trace_in;
    trace_in = NULL;
  }
  
  if (trace_out!=NULL) {
    delete [] trace_out;
    trace_out = NULL;
  }
  
}


void ghost_data::set_basic_info(deghost_basic_info & basic_info) {
  
  nt = basic_info.get_nt();
  dt = basic_info.get_dt();
  zmax = basic_info.get_zmax();
  zmin = basic_info.get_zmin();
  vwater = basic_info.get_vwater();
  epsilon = basic_info.get_epsilon();
  addback = basic_info.get_addback();
  
  if (trace_in!=NULL) {
    delete [] trace_in;
    trace_in = NULL;
  }
  trace_in = new float[nt];
  memset(trace_in,0x00,nt*sizeof(float));
  
  
  if (trace_out!=NULL) {
    delete [] trace_out;
    trace_out = NULL;
  }
  trace_out = new float[nt];
  memset(trace_out,0x00,nt*sizeof(float));
  
  
}


int ghost_data::get_nt() {
  return nt;
}


float ghost_data::get_dt() {
  return dt;
}


float ghost_data::get_zmax() {
  return zmax;
}


float ghost_data::get_zmin() {
  return zmin;
}


float ghost_data::get_vwater() {
  return vwater;
}


float ghost_data::get_epsilon() {
  return epsilon;
}


float ghost_data::get_addback() {
  return addback;
}


void ghost_data::set_input_trace(float* input_trace, int length) {
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    trace_in[i_iter] = input_trace[i_iter];
  }
          
}


void ghost_data::get_output_trace(float* output_trace, int length) {
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    output_trace[i_iter] = trace_out[i_iter];
  }
  
}


void ghost_data::get_input_trace(float* input_trace, int length) {
  
  for (int i_iter = 0; i_iter < length; i_iter++) {
    input_trace[i_iter] = trace_in[i_iter];
  }
  
}


void ghost_data::set_output_trace(float* output_trace, int length) {
 
  for (int i_iter = 0; i_iter < length; i_iter++) {
    trace_out[i_iter] = output_trace[i_iter];
  }
  
}




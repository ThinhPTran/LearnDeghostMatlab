/* 
 * File:   ghost_data_block.C
 * Author: thinh
 * 
 * Created on September 19, 2013, 3:51 PM
 */

#include "ghost_data_block.H"
#include "window_deghost_process_basic_info.H"


ghost_data_block::ghost_data_block() {
  
  fidx = 0;
  lidx = 0;
  
}


ghost_data_block::~ghost_data_block() {
  
  fidx = 0;
  lidx = 0;
  
}


void ghost_data_block::set_fidx(int fidx) {
  this->fidx = fidx;
}


void ghost_data_block::set_lidx(int lidx) {
  this->lidx = lidx;
}


int ghost_data_block::get_fidx() {
  return fidx;
}


int ghost_data_block::get_lidx() {
  return lidx;
}


void ghost_data_block::set_basic_info(window_deghost_process_basic_info & basic_info) {
  
  nt = basic_info.get_ntwindow();
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


/* 
 * File:   window_process_deghosting.H
 * Author: thinh
 *
 * Created on September 20, 2013, 10:20 AM
 */

#ifndef WINDOW_PROCESS_DEGHOSTING_H
#define	WINDOW_PROCESS_DEGHOSTING_H

#include "iostream"
#include "cstring"
#include "General_deghost_engine.H"
#include "window_deghost_process_basic_info.H"
#include "ghost_data_block.H"
#include <omp.h>


using namespace std;


class window_process_deghosting {
  
public:
  window_process_deghosting();
  ~window_process_deghosting();
  
  void display_info();
  void set_basic_info(window_deghost_process_basic_info & info);
  void execute();
  void finalize();
  void set_input_trace(float * input_trace, int length);
  void get_output_trace(float * output_trace, int length);
  
private:
  template<class T> 
  void print_vec ( T* vec, int length );
  void blocking();
  void inversed_blocking();
  void calculate_basic_info();
  void prepare_mem();
  
  
private:
  int first_run;
  int realloc_flag;
  General_deghost_engine engine;
  window_deghost_process_basic_info basic_info;
  int nw;
  ghost_data_block * array_of_blocks;
  
  float * weight;
  float * trace_in;
  float * trace_out;
  
};

#endif	/* WINDOW_PROCESS_DEGHOSTING_H */


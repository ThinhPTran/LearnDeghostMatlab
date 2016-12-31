/* 
 * File:   ghost_data.H
 * Author: thinh
 *
 * Created on September 19, 2013, 1:15 PM
 */

#ifndef GHOST_DATA_H
#define	GHOST_DATA_H


#include <iostream>
#include <cstring>
#include "deghost_basic_info.H"


using namespace std;
  

class ghost_data {
  
public:
  ghost_data();
  ~ghost_data();
  
  void set_basic_info(deghost_basic_info & basic_info);
  int get_nt();
  float get_dt();
  float get_zmax();
  float get_zmin();
  float get_vwater();
  float get_epsilon();
  float get_addback();
  void set_input_trace(float * input_trace, int length);
  void get_output_trace(float * output_trace, int length);
  void get_input_trace(float * input_trace, int length);
  void set_output_trace(float * output_trace, int length);
  
  
public:
  int qualified;;
  
  
protected:
  int nt;
  float dt;
  float zmax;
  float zmin;
  float vwater;
  float epsilon;
  float addback;
  
  float * trace_in;
  float * trace_out;  

};


#endif	/* GHOST_DATA_H */


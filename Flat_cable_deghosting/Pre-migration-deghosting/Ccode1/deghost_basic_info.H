/* 
 * File:   deghost_basic_info.H
 * Author: thinh
 *
 * Created on September 19, 2013, 1:57 PM
 */

#ifndef DEGHOST_BASIC_INFO_H
#define	DEGHOST_BASIC_INFO_H


#include <iostream>



class deghost_basic_info {

public:    
  deghost_basic_info();
  ~deghost_basic_info();
  
  void set_nt(int nt);
  void set_dt(float dt);
  void set_zmax(float zmax);
  void set_zmin(float zmin);
  void set_vwater(float vwater);
  void set_epsilon(float epsilon);
  void set_addback(float addback);
  int get_nt();
  float get_dt();
  float get_zmax();
  float get_zmin();
  float get_vwater();
  float get_epsilon();
  float get_addback();
  
    
private:
  int nt;
  float dt;
  float zmax;
  float zmin;
  float vwater;
  float epsilon;
  float addback;
  
};



#endif	/* DEGHOST_BASIC_INFO_H */


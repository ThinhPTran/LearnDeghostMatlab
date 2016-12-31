/* 
 * File:   window_deghost_process_basic_info.C
 * Author: thinh
 * 
 * Created on September 20, 2013, 9:43 AM
 */

#include "window_deghost_process_basic_info.H"

window_deghost_process_basic_info::window_deghost_process_basic_info() {
  
  ntwindow = 0;
  wshift = 0;
  
}



window_deghost_process_basic_info::~window_deghost_process_basic_info() {

  ntwindow = 0;
  wshift = 0;
  
}


void window_deghost_process_basic_info::set_ntwindow(int ntwindow) {
  this->ntwindow = ntwindow;
}


void window_deghost_process_basic_info::set_wshift(int wshift) {
  this->wshift = wshift;
}


int window_deghost_process_basic_info::get_ntwindow() {
  return ntwindow;
}


int window_deghost_process_basic_info::get_wshift() {
  return wshift;
}






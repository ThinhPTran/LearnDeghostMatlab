/* 
 * File:   window_deghost_process_basic_info.H
 * Author: thinh
 *
 * Created on September 20, 2013, 9:43 AM
 */

#ifndef WINDOW_DEGHOST_PROCESS_BASIC_INFO_H
#define	WINDOW_DEGHOST_PROCESS_BASIC_INFO_H

#include "deghost_basic_info.H"


class window_deghost_process_basic_info : public deghost_basic_info {
  
public:
  window_deghost_process_basic_info();
  ~window_deghost_process_basic_info();
  
  void set_ntwindow(int ntwindow);
  void set_wshift(int wshift);
  int get_ntwindow();
  int get_wshift();
  
  
private:
  int ntwindow; // the length of window to process
  int wshift;  // The shift between the new window and the previous one
  
};

#endif	/* WINDOW_DEGHOST_PROCESS_BASIC_INFO_H */


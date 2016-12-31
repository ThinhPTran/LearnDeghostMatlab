/* 
 * File:   ghost_data_block.H
 * Author: thinh
 *
 * Created on September 19, 2013, 3:51 PM
 */

#ifndef GHOST_DATA_BLOCK_H
#define	GHOST_DATA_BLOCK_H


#include <iostream>
#include "window_deghost_process_basic_info.H"
#include "deghost_basic_info.H"
#include "ghost_data.H"


using namespace std;


class ghost_data_block : public ghost_data {
  
public:
  ghost_data_block();
  ~ghost_data_block();
  void set_basic_info(window_deghost_process_basic_info & basic_info);

  void set_fidx(int fidx);
  void set_lidx(int lidx);
  int get_fidx();
  int get_lidx();
  
  
private:
  int fidx;
  int lidx;

};

#endif	/* GHOST_DATA_BLOCK_H */


/* 
 * File:   Taup_info.H
 * Author: dragon
 *
 * Created on December 12, 2013, 11:07 AM
 */

#ifndef TAUP_INFO_H
#define	TAUP_INFO_H


typedef enum {
  FWD_TAUP=1, 
  UNDEFINED = 0,
  INV_TAUP=-1 
} taup_mode_t;


typedef struct TauP2D_info {
  
  int nt;
  int nx;
  int np;
  
  double dt;
  double dx;
  double dp;
  
  double tmin;
  double xmin;
  double pmin;
  double fmax; 
  
  taup_mode_t mode;
  
  double * in;
  double * out;
  
} TauP2D_info_t;


#endif	/* TAUP_INFO_H */


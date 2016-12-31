/* 
 * File:   Taup_info.H
 * Author: dragon
 *
 * Created on December 12, 2013, 11:07 AM
 */

#ifndef TAUP2DDEGHOST_INFO_H
#define	TAUP2DDEGHOST_INFO_H


typedef struct TauP2DDeghost_info {
  
  double vwater;
  double recover_eps;
  double src_depth;
  double rec_depth;
  
  int nt;
  int nx;
  int np;
  int npad;
  
  double dt;
  double dx;
  
  double fx;
  
  int nfb;
  double fmin;
  double fmax;
  double mfapp;
  
  int ifscale;
  int ifaddback;
  int ifreshape;
  
  double * in;
  double * out;
  
} TauP2DDeghost_info_t;


#endif	/* TAUP2DDEGHOST_INFO_H */


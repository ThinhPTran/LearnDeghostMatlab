#include "deghost_basic_info.H"


deghost_basic_info::deghost_basic_info() {
    
  nt = 0;
  dt = 0;
  zmax = 0;
  zmin = 0;
  vwater = 0;
    
}


deghost_basic_info::~deghost_basic_info() {
    
    
    
}


void deghost_basic_info::set_nt(int nt) {
  this->nt = nt;
}


void deghost_basic_info::set_dt(float dt) {
  this->dt = dt;
}


void deghost_basic_info::set_zmax(float zmax) {
  this->zmax = zmax;
}


void deghost_basic_info::set_zmin(float zmin) {
  this->zmin = zmin;  
}
  
  
void deghost_basic_info::set_vwater(float vwater) {
  this->vwater = vwater;
}


void deghost_basic_info::set_epsilon(float epsilon) {
  this->epsilon = epsilon;
}


void deghost_basic_info::set_addback(float addback) {
  this->addback = addback;
}


int deghost_basic_info::get_nt() {
  return nt;
}


float deghost_basic_info::get_dt() {
  return dt;
}


float deghost_basic_info::get_zmax() {
  return zmax;
}


float deghost_basic_info::get_zmin() {
  return zmin;
}


float deghost_basic_info::get_vwater() {
  return vwater;
}


float deghost_basic_info::get_epsilon() {
  return epsilon;
}


float deghost_basic_info::get_addback() {
  return addback;
}






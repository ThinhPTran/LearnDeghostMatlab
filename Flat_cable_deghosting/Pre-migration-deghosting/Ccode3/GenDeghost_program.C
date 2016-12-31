#include <iostream>
#include <fstream>
#include <string>
#include "deghost_basic_info.H"
#include "General_deghost_engine.H"
#include "ghost_data_block.H"
#include "window_deghost_process_basic_info.H"
#include "window_process_deghosting.H"
#include "GenDeghost_program.H"
#include <omp.h>


using namespace std;


int GenDeghost_program 
(
float * trace_in, 
float * trace_out, 
window_deghost_process_basic_info & basic_info
) 
{
   
  int nt = basic_info.get_nt();

  window_process_deghosting solver;
  solver.set_basic_info(basic_info);
  solver.set_input_trace(trace_in,nt);
  solver.execute();
  solver.get_output_trace(trace_out,nt);
  solver.finalize();

  
}

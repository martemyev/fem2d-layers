#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "config.h"
#include "boost/filesystem.hpp"
#include <vector>
#include <string>

//namespace fs = boost::filesystem;


enum TIME_SCHEMES
{
  EXPLICIT,      // explicit scheme
  CRANK_NICOLSON // sort of implicit scheme
};


class Parameters
{
public:
  static const int DIM = 2; // dimension of the task

  int TIME_SCHEME; // explicit, Crank-Nicolson

  double X_BEG, X_END, Z_BEG, Z_END; // 2D computational domain
  double CL; // characteristic length of the triangles

  double TIME_BEG, TIME_END; // time of the beginning of the simulation (usually it's 0), and of the ending
  double TIME_STEP; // time step
  int N_TIME_STEPS; // number of the time steps

  int FE_ORDER; // the order of the FE basis functions
  //int QUAD_ORDER; // the order of quadrature formula

  typedef std::string path; // or boost::filesystem::path

  path MESH_DIR; // path to the directory where the meshes (.msh) are
  path GEO_DIR; // path to the directory where the geometry file (.geo) are
  path RES_DIR; // path to the directory where all the results of the simulation will be kept
  path VTU_DIR; // the name of the directory where the .vtu files with results on some time steps will be kept
  path SOL_DIR; // the name of the directory where the .sol files with solutions on some time steps will be kept
  path TIME_FILE; // the name of the file where the time of calculations will be kept
  path INFO_FILE; // the name of the file where the information about simulation will be kept
  path MESH_FILE; // the name of the file with a mesh

  bool PRINT_VTU; // whether we need to print .vtu files or not
  bool VTU_STEP; // if we need to print .vtu files we can do that not for each time step, but for every (VTU_STEP)-th step
  bool PRINT_INFO; // whether we need to print some info on the screen during the calculations


  Parameters(int argc = 0, char **argv = 0); // constructor
  ~Parameters();

  std::string print() const;


private:
  Parameters(const Parameters&);
  const Parameters& operator =(const Parameters&);

  void default_parameters(); // initialize default parameters
  void read_from_command_line(int argc, char **argv); // read parameters from command line
  void generate_paths(); // generate the paths to some files and directories
  void check_clean_dirs() const; // check the existance and clearance of some directories
};


#endif // PARAMETERS_H

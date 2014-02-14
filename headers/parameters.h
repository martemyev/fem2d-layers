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

  static const int INCL_DOMAIN = 11; // the number of domain that characterizes the inclusions

  int TIME_SCHEME; // explicit, Crank-Nicolson

  double X_BEG, X_END, Y_BEG, Y_END; // 2D computational domain
  double CL; // characteristic length of the triangles

  double TIME_BEG, TIME_END; // time of the beginning of the simulation (usually it's 0), and of the ending
  double TIME_STEP; // time step
  int N_TIME_STEPS; // number of the time steps

  int FE_ORDER; // the order of the FE basis functions
  //int QUAD_ORDER; // the order of quadrature formula

  double SOURCE_FREQUENCY; // the frequency of the Ricket wavelet, that is used as a source of seismic waves
  double SOURCE_SUPPORT; // the support of the source function
  double SOURCE_CENTER_X; // x-coordinate of the source center
  double SOURCE_CENTER_Y; // y-coordinate of the source center

  // equation is:
  // alpha * u'' - div(beta * grad u) = f
  double COEF_ALPHA_1_VALUE; // the value of the homogeneous ALPHA_1 coefficient
  double COEF_ALPHA_2_VALUE; // the value of the homogeneous ALPHA_2 coefficient
  std::string COEF_BETA_1_FILE; // the name of the file where coefficient beta in the main domain (not in inclusions) distribution is represented.
                                // it can be empty string which means that coef_beta_1 is homogeneous
                                // in the main domain. in this case COEF_BETA_1_VALUE is used.
  std::string COEF_BETA_2_FILE; // the same as COEF_BETA_1_FILE but for inclusions
  double COEF_BETA_1_VALUE; // the value of the homogeneous BETA_1 coefficient
  double COEF_BETA_2_VALUE; // the value of the homogeneous BETA_2 coefficient

  double INCL_CENTER_X; // x-coordinate of the center of the inclusion
  double INCL_CENTER_Y; // y-coordinate of the center of the inclusion
  double INCL_RADIUS;   // radius of the circular inclusion

  std::string MESH_DIR; // path to the directory where the meshes (.msh) are
  std::string GEO_DIR; // path to the directory where the geometry file (.geo) are
  std::string RES_DIR; // path to the directory where all the results of the simulation will be kept
  std::string VTU_DIR; // the name of the directory where the .vtu files with results on some time steps will be kept
  std::string SOL_DIR; // the name of the directory where the .sol files with solutions on some time steps will be kept
  std::string TIME_FILE; // the name of the file where the time of calculations will be kept
  std::string INFO_FILE; // the name of the file where the information about simulation will be kept
  std::string MESH_FILE; // the name of the file with a mesh

  bool PRINT_VTU; // whether we need to print .vtu files or not
  bool PRINT_INFO; // whether we need to print some info on the screen during the calculations
  bool SAVE_SOL; // wheather wee need to save the solution of every time step (for further comparison, for example)
  unsigned int VTU_STEP; // if we need to print .vtu files we can do that not for each time step, but for every (VTU_STEP)-th step
  unsigned int SOL_STEP; // if we need to save .dat files we can do that not for each time step, but for every (SOL_STEP)-th step


  Parameters(int argc = 0, char **argv = 0); // constructor
  ~Parameters();

  std::string print() const;

            /**
             * To make some useful things for calculations.
             * For example, create the structure of output directories, files, etc.
             * This procedure was separated from Parameter constructor to
             * avoid time-consuming operations with filesystem in default parameters objects
             * (that are used, for example, in testing procedures).
             * This function has to be called before real work!
             */
  void establish_environment();


private:
  Parameters(const Parameters&);
  const Parameters& operator =(const Parameters&);

  void default_parameters(); // initialize default parameters
  void read_from_command_line(int argc, char **argv); // read parameters from command line
  void generate_paths(); // generate the paths to some files and directories
  void check_clean_dirs() const; // check the existance and clearance of some directories
};


#endif // PARAMETERS_H

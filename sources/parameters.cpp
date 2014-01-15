#include "parameters.h"
#include "config.h"
#include "auxiliary_functions.h"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include <iostream>

namespace po = boost::program_options;


Parameters::Parameters(int argc, char **argv)
{
  default_parameters(); // initialize all parameters
  if (argc > 1)
    read_from_command_line(argc, argv); // change some (or all) parameters from default to custom ones
  generate_paths(); // generate all necessary path to files and directories
  check_clean_dirs(); // check the existance and clearance of some directories
}



void Parameters::default_parameters()
{
  MESH_DIR = "/u/artemyev/projects/tat_gmsfem/brandnew/meshes/";
  GEO_DIR = "/u/artemyev/projects/tat_gmsfem/brandnew/geometries/";

  RES_DIR = ""; // should be changed based on other parameters
  VTU_DIR = "vtu/"; // should be added to RES_DIR after generating of the latter
  SOL_DIR = "sol/"; // should be added to RES_DIR after generating of the latter
  TIME_FILE = "time.txt"; // should be added to RES_DIR after generating of the latter
  INFO_FILE = "info.txt"; // should be added to RES_DIR after generating of the latter

  MESH_FILE = "mesh.msh";  // should be added to MESH_DIR after establishing of the latter (means that MESH_DIR can be changed from parameter file of command line)

  TIME_SCHEME = EXPLICIT;
  X_BEG = Z_BEG = 0.;
  X_END = Z_END = 1.;
  CL = 0.1;
  TIME_BEG = 0.;
  TIME_END = 1.;
  TIME_STEP = 1.;
  N_TIME_STEPS = 1;
  FE_ORDER = 1;
  //QUAD_ORDER = 3;

  PRINT_VTU = 1; // print .vtu files by default
  VTU_STEP = 1; // print the .vtu file on each time step
  PRINT_INFO = 1; // print an information to console on each time step
}



void Parameters::read_from_command_line(int argc, char **argv)
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("mesh", po::value<std::string>(), "name of mesh file")
    ("scheme", po::value<std::string>(), "time scheme")
    ("tend", po::value<double>(), "time ending")
    ("ts", po::value<double>(), "time step")
    ("nt", po::value<unsigned int>(), "number of time steps")
    ("fe", po::value<unsigned int>(), "order of fe basis functions")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    exit(0);
  }

  if (vm.count("mesh"))
    MESH_FILE = vm["mesh"].as<std::string>();

  if (vm.count("scheme"))
  {
    std::string scheme_name = vm["scheme"].as<std::string>();
    if (scheme_name == "explicit")
      TIME_SCHEME = EXPLICIT;
    else if (scheme_name == "cn")
      TIME_SCHEME = CRANK_NICOLSON;
    else
      require(false, "Unknown time scheme : " + scheme_name);
  }

  if (vm.count("ts") && vm.count("nt") && vm.count("tend"))
    require(false, "ts, nt and tend parameters cannot be used together - maximum two of them");

  if (vm.count("tend"))
    TIME_END = vm["tend"].as<double>();
  if (vm.count("ts"))
    TIME_STEP = vm["ts"].as<double>();
  if (vm.count("nt"))
    N_TIME_STEPS = vm["nt"].as<unsigned int>();

  if (!vm.count("nt"))
    N_TIME_STEPS = int((TIME_END - TIME_BEG) / TIME_STEP);
  if (!vm.count("ts"))
    TIME_STEP = (TIME_END - TIME_BEG) / N_TIME_STEPS;
  if (!vm.count("tend"))
    TIME_END = TIME_BEG + N_TIME_STEPS * TIME_STEP;

  require(TIME_BEG >= 0., "time of begin cannot be negative");
  require(TIME_END > TIME_BEG, "time of end is not bigger than time of begin");
  require(TIME_STEP <= (TIME_END - TIME_BEG), "time step cannot be bigger than the whole time");
  require(N_TIME_STEPS > 0, "the number of time steps cannot be less than 1");

  require(fabs((TIME_END - TIME_BEG - N_TIME_STEPS * TIME_STEP) / TIME_END) < 1e-14,
          "time parameters (TIME_END, TIME_STEP and N_TIME_STEPS) do not conform to each other");

  if (vm.count("fe"))
    FE_ORDER = vm["fe"].as<unsigned int>();

  require(FE_ORDER == 1, "This order of basis functions (" + d2s(FE_ORDER) + ") is not implemented");
}



Parameters::~Parameters()
{ }



std::string Parameters::print() const
{
  std::string time_scheme_name[] = { "explicit",
                                     "Crank-Nicolson" };

  std::string str = "list of parameters:\n";
  str += "dim = " + d2s(DIM) + "\n";
  str += "scheme = " + time_scheme_name[TIME_SCHEME] + "\n";
  str += "mesh file name = " + MESH_FILE + "\n";
  str += "mesh cl = " + d2s(CL) + "\n";
  str += "domain = [" + d2s(X_BEG) + ", " + d2s(X_END) + "] x [" + d2s(Z_BEG) + ", " + d2s(Z_END) + "]\n";
  str += "time = from " + d2s(TIME_BEG) + " to " + d2s(TIME_END) + " sec\n";
  str += "time step = " + d2s(TIME_STEP) + "\n";
  str += "number of time steps = " + d2s(N_TIME_STEPS) + "\n";
  str += "fe basis order = " + d2s(FE_ORDER) + "\n";
  return str;
}



void Parameters::generate_paths()
{
  //MESH_FILE = MESH_DIR + "/" + MESH_FILE; // full path to the mesh file

  //RES_DIR

//  RES_DIR = RESULTS_DIR + "
}



void Parameters::check_clean_dirs() const
{
  //using boost::filesystem;
  //path top_res_dir(RESULTS_DIR);
  //if (!exists(top_res_dir)) // if this top directory for results doesn't exist, we create it

}

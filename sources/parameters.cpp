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
}



void Parameters::default_parameters()
{
  MESH_DIR = "/u/artemyev/projects/tat_gmsfem/brandnew/meshes/";
  GEO_DIR = "/u/artemyev/projects/tat_gmsfem/brandnew/geometries/";

  RES_DIR = ""; // should be changed and based on some parameters
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
  SOURCE_FREQUENCY = 20;
  SOURCE_SUPPORT = 10;

  COEF_A_1_FILE = "";
  COEF_A_2_FILE = "";
  COEF_A_1_VALUE = 1.;
  COEF_A_2_VALUE = 1.;

  PRINT_VTU = 1; // print .vtu files by default
  PRINT_SOL = 0; // don't print .sol files by default
  VTU_STEP = 1; // print the .vtu file on each time step
  PRINT_INFO = 1; // print an information to console on each time step
}



void Parameters::read_from_command_line(int argc, char **argv)
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("meshfile", po::value<std::string>(), "name of mesh file")
    ("scheme", po::value<std::string>(), "time scheme")
    ("tend", po::value<double>(), "time ending")
    ("tstep", po::value<double>(), "time step")
    ("nt", po::value<unsigned int>(), "number of time steps")
    ("fe", po::value<unsigned int>(), "order of fe basis functions")
    ("a1file", po::value<std::string>(), "name of the file with coefficient a in main domain distribution")
    ("a2file", po::value<std::string>(), "name of the file with coefficient a in inclusions distribution")
    ("a1val", po::value<double>(), "the value of coefficient a in main homogeneous domain")
    ("a2val", po::value<double>(), "the value of coefficient a in homogeneous inclusions")
    ("vtu", po::value<bool>(), "whether we need to print .vtu files")
    ("sol", po::value<bool>(), "whether we need to print .sol files")
    ("inf", po::value<bool>(), "whether we need to print some info during calculations")
    ("vtu_step", po::value<unsigned int>(), "if we need to print .vtu files then how often. every (vtu_step)-th file will be printed")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    exit(0);
  }

  if (vm.count("meshfile"))
    MESH_FILE = vm["meshfile"].as<std::string>();

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

  require(!(vm.count("tstep") && vm.count("nt") && vm.count("tend")),
          "tstep, nt and tend parameters cannot be used together - maximum two of them");

  if (vm.count("tend"))
    TIME_END = vm["tend"].as<double>();
  if (vm.count("tstep"))
    TIME_STEP = vm["tstep"].as<double>();
  if (vm.count("nt"))
    N_TIME_STEPS = vm["nt"].as<unsigned int>();

  if (!vm.count("nt"))
    N_TIME_STEPS = int((TIME_END - TIME_BEG) / TIME_STEP);
  if (!vm.count("tstep"))
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

  require(!(vm.count("a1file") && vm.count("a1val")), "a1file and a1val can not be used together - only one of them");
  require(!(vm.count("a2file") && vm.count("a2val")), "a2file and a2val can not be used together - only one of them");

  if (vm.count("a1file"))
    COEF_A_1_FILE = vm["a1file"].as<std::string>();
  if (vm.count("a2file"))
    COEF_A_2_FILE = vm["a2file"].as<std::string>();
  if (vm.count("a1val"))
    COEF_A_1_VALUE = vm["a1val"].as<double>();
  if (vm.count("a2val"))
    COEF_A_2_VALUE = vm["a2val"].as<double>();

  if (vm.count("vtu"))
    PRINT_VTU = vm["vtu"].as<bool>();
  if (vm.count("sol"))
    PRINT_SOL = vm["sol"].as<bool>();
  if (vm.count("inf"))
    PRINT_INFO = vm["inf"].as<bool>();
  if (vm.count("vtu_step"))
    VTU_STEP = vm["vtu_step"].as<unsigned int>();
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



void Parameters::establish_environment()
{
  generate_paths(); // generate all necessary path to files and directories
  check_clean_dirs(); // check the existance and clearance of some directories
}



void Parameters::generate_paths()
{
  using namespace boost::filesystem;

  MESH_FILE = MESH_DIR + "/" + MESH_FILE; // full path to the mesh file
  path mesh_file(MESH_FILE);

  std::string coef_a1, coef_a2; // description of coefficients a
  if (COEF_A_1_FILE == "")
    coef_a1 = d2s(COEF_A_1_VALUE);
  else
  {
    path coef_name(COEF_A_1_FILE);
    coef_a1 = coef_name.stem().string();
  }
  if (COEF_A_2_FILE == "")
    coef_a2 = d2s(COEF_A_2_VALUE);
  else
  {
    path coef_name(COEF_A_2_FILE);
    coef_a2 = coef_name.stem().string();
  }
  RES_DIR = RESULTS_DIR + "/" + mesh_file.stem().string() +
            "_T" + d2s(TIME_END) + "_K" + d2s(N_TIME_STEPS) +
            "_f" + d2s(SOURCE_FREQUENCY) + "_P" + d2s(SOURCE_SUPPORT) +
            "_A" + coef_a1 + "_" + coef_a2 + "/";

  VTU_DIR = RES_DIR + "/" + VTU_DIR;
  SOL_DIR = RES_DIR + "/" + SOL_DIR;
  TIME_FILE = RES_DIR + "/" + TIME_FILE;
  INFO_FILE = RES_DIR + "/" + INFO_FILE;
}



void Parameters::check_clean_dirs() const
{
  using namespace boost::filesystem;

  path top_res_dir(RESULTS_DIR); // top directory with results (from config.h, CMakeLists.txt)
  if (!exists(top_res_dir)) // if this top directory for results doesn't exist, we create it
    create_directory(top_res_dir);

  require(RES_DIR != "", "Directory for results has no name");
  path cur_res_dir(RES_DIR); // current directory with results
  if (exists(cur_res_dir) && is_directory(cur_res_dir)) // if this directory exists, we need to clean it up
    remove_all(cur_res_dir); // remove all contents of the directory and the directory itself
  create_directory(cur_res_dir); // now create empty directory

  if (PRINT_VTU)
    create_directory(VTU_DIR);

  if (PRINT_SOL)
    create_directory(SOL_DIR);
}

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
  SOL_DIR = "dat/"; // should be added to RES_DIR after generating of the latter
  TIME_FILE = "time.txt"; // should be added to RES_DIR after generating of the latter
  INFO_FILE = "info.txt"; // should be added to RES_DIR after generating of the latter

  MESH_FILE = "mesh.msh";  // should be added to MESH_DIR after establishing of the latter (means that MESH_DIR can be changed from parameter file of command line)

  TIME_SCHEME = EXPLICIT;
  X_BEG = Y_BEG = 0.;
  X_END = Y_END = 1.;
  CL = 0.1;
  TIME_BEG = 0.;
  TIME_END = 1.;
  TIME_STEP = 1.;
  N_TIME_STEPS = 1;
  FE_ORDER = 1;
  //QUAD_ORDER = 3;
  SOURCE_FREQUENCY = 20;
  SOURCE_SUPPORT = 10;
  SOURCE_CENTER_X = 0;
  SOURCE_CENTER_Y = 0;

  COEF_BETA_1_FILE = "";
  COEF_BETA_2_FILE = "";
  COEF_BETA_1_VALUE = 1.;
  COEF_BETA_2_VALUE = 1.;
  COEF_ALPHA_1_VALUE = 1.;
  COEF_ALPHA_2_VALUE = 1.;

  INCL_CENTER_X = 1000;
  INCL_CENTER_Y = 1500;
  INCL_RADIUS = 0; // once the radius is determined, the program starts to calculate taking this inclusion into account

  PRINT_VTU = 0; // print .vtu files by default (the last time step will be printed in any case)
  VTU_STEP = 1; // print the .vtu file on each time step
  PRINT_INFO = 0; // print an information to console on each time step
  SAVE_SOL = 0; // save solution or not (the last time step will be saved in any case)
  SOL_STEP = 1; // save solution on every time step
}



void Parameters::read_from_command_line(int argc, char **argv)
{
  std::string time_scheme = (TIME_SCHEME == EXPLICIT ? "explicit" : "crank-nicolson");

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("meshfile",  po::value<std::string>(),   std::string("name of mesh file (" + MESH_FILE + ")").c_str())
    ("scheme",    po::value<std::string>(),   std::string("time scheme (" + time_scheme + ")").c_str())
    ("tend",      po::value<double>(),        std::string("time ending (" + d2s(TIME_END) + ")").c_str())
    ("tstep",     po::value<double>(),        std::string("time step (" + d2s(TIME_STEP) + ")").c_str())
    ("nt",        po::value<unsigned int>(),  std::string("number of time steps (" + d2s(N_TIME_STEPS) + ")").c_str())
    ("fe",        po::value<unsigned int>(),  std::string("order of fe basis functions (" + d2s(FE_ORDER) + ")").c_str())
    ("a1val",     po::value<double>(),        std::string("the value of coefficient alpha in main homogeneous domain (" + d2s(COEF_ALPHA_2_VALUE) + ")").c_str())
    ("a2val",     po::value<double>(),        std::string("the value of coefficient alpha in homogeneous inclusions (" + d2s(COEF_ALPHA_2_VALUE) + ")").c_str())
    ("b1val",     po::value<double>(),        std::string("the value of coefficient beta in main homogeneous domain (" + d2s(COEF_BETA_1_VALUE) + ")").c_str())
    ("b2val",     po::value<double>(),        std::string("the value of coefficient beta in homogeneous inclusions (" + d2s(COEF_BETA_2_VALUE) + ")").c_str())
    ("b1file",    po::value<std::string>(),   std::string("name of the file with coefficient a in main domain distribution (" + COEF_BETA_1_FILE + ")").c_str())
    ("b2file",    po::value<std::string>(),   std::string("name of the file with coefficient a in inclusions distribution (" + COEF_BETA_2_FILE + ")").c_str())
    ("vtu",       po::value<bool>(),          std::string("whether we need to print .vtu files (" + d2s(PRINT_VTU) + ")").c_str())
    ("sol",       po::value<bool>(),          std::string("whether we need to save solutions in .dat files (" + d2s(SAVE_SOL) + ")").c_str())
    ("inf",       po::value<bool>(),          std::string("whether we need to print some info during calculations (" + d2s(PRINT_INFO) + ")").c_str())
    ("vtu_step",  po::value<unsigned int>(),  std::string("if we need to print .vtu files then how often. every (vtu_step)-th file will be printed (" + d2s(VTU_STEP) + ")").c_str())
    ("sol_step",  po::value<unsigned int>(),  std::string("if we need to save .dat files then how often. every (sol_step)-th file will be saved (" + d2s(SOL_STEP) + ")").c_str())
    ("f0",        po::value<double>(),        std::string("source frequency (" + d2s(SOURCE_FREQUENCY) + ")").c_str())
    ("p",         po::value<double>(),        std::string("source support (" + d2s(SOURCE_SUPPORT) + ")").c_str())
    ("xcen",      po::value<double>(),        std::string("source center, x-coordinate (" + d2s(SOURCE_CENTER_X) + ")").c_str())
    ("ycen",      po::value<double>(),        std::string("source center, y-coordinate (" + d2s(SOURCE_CENTER_Y) + ")").c_str())
    ("x1",        po::value<double>(),        std::string("x-coordinate of the limit (" + d2s(X_END) + ")").c_str())
    ("y1",        po::value<double>(),        std::string("y-coordinate of the limit (" + d2s(Y_END) + ")").c_str())
    ("inclx",     po::value<double>(),        std::string("x-coordinate of the center of the inclusion (" + d2s(INCL_CENTER_X) + ")").c_str())
    ("incly",     po::value<double>(),        std::string("y-coordinate of the center of the inclusion (" + d2s(INCL_CENTER_Y) + ")").c_str())
    ("inclr",     po::value<double>(),        std::string("radius of the circular inclusion (" + d2s(INCL_RADIUS) + ")").c_str())
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

  require(!(vm.count("b1file") && vm.count("b1val")), "b1file and b1val can not be used together - only one of them");
  require(!(vm.count("b2file") && vm.count("b2val")), "b2file and b2val can not be used together - only one of them");

  if (vm.count("a1val"))
    COEF_ALPHA_1_VALUE = vm["a1val"].as<double>();
  if (vm.count("a2val"))
    COEF_ALPHA_2_VALUE = vm["a2val"].as<double>();
  if (vm.count("b1val"))
    COEF_BETA_1_VALUE = vm["b1val"].as<double>();
  if (vm.count("b2val"))
    COEF_BETA_2_VALUE = vm["b2val"].as<double>();
  if (vm.count("b1file"))
    COEF_BETA_1_FILE = vm["b1file"].as<std::string>();
  if (vm.count("b2file"))
    COEF_BETA_2_FILE = vm["b2file"].as<std::string>();

  if (vm.count("vtu"))
    PRINT_VTU = vm["vtu"].as<bool>();
  if (vm.count("sol"))
    SAVE_SOL = vm["sol"].as<bool>();
  if (vm.count("inf"))
    PRINT_INFO = vm["inf"].as<bool>();
  if (vm.count("vtu_step"))
    VTU_STEP = vm["vtu_step"].as<unsigned int>();
  if (vm.count("sol_step"))
    SOL_STEP = vm["sol_step"].as<unsigned int>();

  if (vm.count("f0"))
    SOURCE_FREQUENCY = vm["f0"].as<double>();
  if (vm.count("p"))
    SOURCE_SUPPORT = vm["p"].as<double>();
  if (vm.count("xcen"))
    SOURCE_CENTER_X = vm["xcen"].as<double>();
  if (vm.count("ycen"))
    SOURCE_CENTER_Y = vm["ycen"].as<double>();

  if (vm.count("x1"))
    X_END = vm["x1"].as<double>();
  if (vm.count("y1"))
    Y_END = vm["y1"].as<double>();

  if (vm.count("inclx"))
    INCL_CENTER_X = vm["inclx"].as<double>();
  if (vm.count("incly"))
    INCL_CENTER_Y = vm["incly"].as<double>();
  if (vm.count("inclr"))
    INCL_RADIUS = vm["inclr"].as<double>();
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
  //str += "mesh cl = " + d2s(CL) + "\n";
  //str += "domain = [" + d2s(X_BEG) + ", " + d2s(X_END) + "] x [" + d2s(Y_BEG) + ", " + d2s(Y_END) + "]\n";
  str += "time = from " + d2s(TIME_BEG) + " to " + d2s(TIME_END) + " sec\n";
  str += "time step = " + d2s(TIME_STEP) + "\n";
  str += "number of time steps = " + d2s(N_TIME_STEPS) + "\n";
  //str += "fe basis order = " + d2s(FE_ORDER) + "\n";
  str += "f0 = " + d2s(SOURCE_FREQUENCY) + "\n";
  str += "P = " + d2s(SOURCE_SUPPORT) + "\n";
  str += "x_cen = " + d2s(SOURCE_CENTER_X) + "\n";
  str += "y_cen = " + d2s(SOURCE_CENTER_Y) + "\n";
  str += "a1 = " + d2s(COEF_ALPHA_1_VALUE) + "\n";
  str += "a2 = " + d2s(COEF_ALPHA_2_VALUE) + "\n";
  str += "b1 = " + d2s(COEF_BETA_1_VALUE) + "\n";
  str += "b2 = " + d2s(COEF_BETA_2_VALUE) + "\n";
  str += "b1file = " + COEF_BETA_1_FILE + "\n";
  str += "b2file = " + COEF_BETA_2_FILE + "\n";
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

  std::string coef_b1, coef_b2; // description of coefficients beta
  if (COEF_BETA_1_FILE == "")
    coef_b1 = d2s(COEF_BETA_1_VALUE);
  else
  {
    path coef_name(COEF_BETA_1_FILE);
    coef_b1 = coef_name.stem().string();
  }
  if (COEF_BETA_2_FILE == "")
    coef_b2 = d2s(COEF_BETA_2_VALUE);
  else
  {
    path coef_name(COEF_BETA_2_FILE);
    coef_b2 = coef_name.stem().string();
  }
  RES_DIR = RESULTS_DIR + "/" + mesh_file.stem().string() +
            "_T" + d2s(TIME_END) + "_K" + d2s(N_TIME_STEPS) +
            "_f" + d2s(SOURCE_FREQUENCY) + "_P" + d2s(SOURCE_SUPPORT) +
            "_A" + d2s(COEF_ALPHA_1_VALUE) + "_" + d2s(COEF_ALPHA_2_VALUE) +
            "_B" + coef_b1 + "_" + coef_b2 + "/";

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

  //if (PRINT_VTU)
    create_directory(VTU_DIR); // always create a VTU directory, since we always print the results of simulation at the last time step

  //if (SAVE_SOL)
    create_directory(SOL_DIR); // always create a SOL directory, since we always save the results of simulation at the last time step
}

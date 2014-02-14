#include "config.h"
#include "parameters.h"
#include "acoustic2d.h"
#include "testing.h"
#include <iostream>
#include <boost/timer/timer.hpp>
#include <gtest/gtest.h>
#include "petscsys.h"


int main(int argc, char **argv)
{
  PetscInitialize(&argc, &argv, NULL, NULL);

  Parameters param(argc, argv);
  param.establish_environment();

  // time measurement
  boost::timer::auto_cpu_timer boost_timer;

#if defined(TESTING)
  std::cout << "\n\nTESTING\n";
  ::testing::InitGoogleTest(&argc, argv);
  int test_ret = RUN_ALL_TESTS();
  std::cout << "\nTesting procedures finished (" << test_ret << " is returned)\n\n";
#endif

//#if defined(DEBUG)
  std::cout << param.print() << std::endl;
//#endif

  Acoustic2D problem(&param);
  problem.solve();

  PetscFinalize();

  return 0;
}


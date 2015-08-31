#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "kernel.h"
#include "gridgeneration.h"
#include <cmath>
#include <random>
#include "MeshMovement.h"

#include "cuda.h"
#include <cuda_runtime.h>

//Test cases
#include "test_list.h"

//Main program
int main(int argc, char *argv[]) {			
  TestCases1D::TestCase1D_SodShockTube test;
  test.Init(argc, argv);
  test.Run(200);
	return 0;
};

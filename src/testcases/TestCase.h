#ifndef TURBO_TestCases_TestCase
#define TURBO_TestCases_TestCase

#include "kernel.h"

//Base class for all automated tests
class TestCase {
protected:
	Kernel* _kernel; //Computational kernel object	
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
public:	
	virtual void RunTest(int* argc, char** argv[]) = 0; //Run test with program arguments
	//virtual void RunTest(Kernel* kernel) = 0; //Main interface function for running test case code
};


#endif
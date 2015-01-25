#ifndef TURBO_TestCases_TestCase
#define TURBO_TestCases_TestCase

#include "kernel.h"

//Base class for all automated tests
class TestCase {
protected:
	Kernel _kernel; //Computational kernel object
public:
	void RunTest(Kernel* kernel); //Main interface function for running test case code
};


#endif
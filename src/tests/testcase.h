#ifndef TURBO_TESTS_TESTCASE
#define TURBO_TESTS_TESTCASE

#include "kernel.h"

//Base class for test cases
class TestCase {
public:
	//Main variables
	int argc; 
	char** argv;

	//Test description info
	std::string TestInfo;

	//Test results report
	std::string TestReport;

	//Destructor
	virtual ~TestCase() { /* Empty implementation */ };

	//Main test code function
	virtual bool RunTest() = 0; 
};

#endif
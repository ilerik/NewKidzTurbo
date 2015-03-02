#ifndef TURBO_TestCases_TestCase
#define TURBO_TestCases_TestCase

#include "kernel.h"

struct TestCaseResultInfo {
	//Test status
	enum class TestStatusType {
		Passed,
		Failed
	} testStatus;
};

//Base class for all automated tests
class TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
public:	
	virtual void RunTestWithKernel(Kernel* kernel) = 0; //Main interface function for running test case code

	virtual TestCaseResultInfo GetTestCaseResultInfo() = 0; //Get results of test run 

	//Attach computational kernel
	void SetKernel(Kernel& kernel) {
		_kernel = &kernel;
	};
};


#endif
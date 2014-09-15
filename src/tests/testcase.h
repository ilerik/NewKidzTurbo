#ifndef TURBO_TESTS_TESTCASE
#define TURBO_TESTS_TESTCASE

#include "kernel.h"

class TestCase {
public:
	//main variables
	int argc; char** argv;

	//constructor
	/*TestCase(int _argc, char *_argv[]) {
		argc = _argc;
		argv = _argv;
	};*/

	//destructor
	~TestCase();

	virtual bool RunTest() {
	return false;
	};
};

#endif
#ifndef TURBO_TestCases_RiemannProblemTests_TestCase1D_SodShockTube
#define TURBO_TestCases_RiemannProblemTests_TestCase1D_SodShockTube

#include "TestCase.h"
#include "gengrid1D.h"

namespace TestCases1D {

class TestCase1D_SodShockTube : public InitialConditions::InitialConditions {
private:
	//Data
	std::unique_ptr<Kernel> _kernel;
	std::unique_ptr<Grid> _grid;
	Configuration _configuration;

	//
	void _prepareGrid(int nCells) {
		//_grid = GenGrid1D();
	};
public:
	
	void PrepareConfiguration();	

	void Run(int nCells) {

	};
};

};

#endif
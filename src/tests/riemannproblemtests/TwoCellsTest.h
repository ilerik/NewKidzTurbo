#ifndef TURBO_TESTS_RIEMANNPROBLEMTEST_TWOCELLSTESTS
#define TURBO_TESTS_RIEMANNPROBLEMTEST_TWOCELLSTESTS

#include "testcase.h"
#include "perfect_gas_roe_3d.h"

class TwoCellsTest : public TestCase
{
public:
	Roe3DSolverPerfectGas RiemanSolver;

	//constructor by default
	TwoCellsTest() {};

	//run test function
	bool RunTest() {
		//Test case description
		TestInfo = "Test for computing Riemann solvers between two cells";
		GasModel::ConservativeVariables UL, UR;

		//set conservative variables
		UL.ro = 1.0;
		UL.rou = 1.0;
		UL.rov = 1.0;
		UL.row = 1.0;
		UL.roE = 1.0;
		UR.ro = 1.0;
		UR.rou = 1.0;
		UR.rov = 1.0;
		UR.row = 1.0;
		UR.roE = 1.0;

		//create a face
		Face f;
		f.FaceNormal = Vector(1, 0, 0);
		f.FaceSquare = 0.5;
		
		//compute flux
		std::vector<double> res = RiemanSolver.ComputeFlux(UL, UR, f);

		//TO DO results of other methods
				
		return true;
	};
};




#endif
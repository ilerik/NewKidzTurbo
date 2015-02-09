#ifndef TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest1
#define TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest1

#include "ToroTest.h"

namespace ToroTests {

class AutomaticTest1 {
public:
	int NCells;		//size of grid
	double total_time;	//time for comparison of solutions

	//constructor
	AutomaticTest1() {
		NCells = 200;
		total_time = 0.2;
	};
	AutomaticTest1(int _NCells, double _total_time) : NCells(_NCells), total_time(_total_time) { };

	//classes of different solvers for Test 1 (E.Toro's book // parag. 4.3.3, pp. 129 - 133)
	//Roe
	class ToroTest1Roe : public ToroTest {
	public:
	
		//constructor
		ToroTest1Roe(ToroTestsInit &params) {
			this->SetParams(params);
			_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Roe;
			solution_name = "test1_Roe";
		};
	}; //ToroTest1

	//HLLC
	class ToroTest1HLLC : public ToroTest {
	public:
	
		//constructor
		ToroTest1HLLC(ToroTestsInit &params) {
			this->SetParams(params);
			_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Roe;
			solution_name = "test1_HLLC";
		};
	}; //ToroTest1

	//Godunov
	class ToroTest1Godunov : public ToroTest {
	public:
	
		//constructor
		ToroTest1Godunov(ToroTestsInit &params) {
			this->SetParams(params);
			_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Roe;
			solution_name = "test1_Godunov";
		};
	}; //ToroTest1

	//TEST 1 (automatic)
	void RunTest(int* argc, char** argv[]) {

		ToroTestsInit test1;
		//fill task settings for Test 1
		test1.discontinuity_pos = 0.5;
		test1.Lx = 1.0;
		test1.nCells = NCells;
		test1.TimeMax = total_time;
		test1.roL = 1.0;
		test1.uL = 0.0;
		test1.pL = 1.0;
		test1.roR = 0.125;
		test1.uR = 0.0;
		test1.pR = 0.1;

		//create tests for all solvers
		Kernel kernel;
		kernel.Initilize(argc, argv);

		ToroTest1Roe test_roe(test1);
		ToroTest1HLLC test_hllc(test1);
		ToroTest1Godunov test_godunov(test1);

		//Run all tests
		std::cout << "Test 1 for Roe solver is running!\n";
		test_roe.RunTestWithKernel(&kernel);
		std::cout << "Test 1 for HLLC solver is running!\n";
		test_hllc.RunTestWithKernel(&kernel);
		std::cout << "Test 1 for Godunov solver is running!\n";
		test_godunov.RunTestWithKernel(&kernel);

		kernel.Finalize();
	};

};

}; //namespace


#endif
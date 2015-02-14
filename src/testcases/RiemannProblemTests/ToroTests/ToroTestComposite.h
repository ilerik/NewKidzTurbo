#ifndef TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest1
#define TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest1

#include "ToroTest.h"

namespace ToroTests {

class ToroTestComposite : public ToroTest {
private:
	std::vector<TestCaseResultInfo::TestStatusType> _statuses;
	std::vector<double> _L2Errors;
	std::vector<double> _LInfErrors;
public:	
	//classes of different solvers for Test 1 (E.Toro's book // parag. 4.3.3, pp. 129 - 133)
	//Roe
	class ToroTestRoe : public ToroTest {
	public:	
		//Constructor overload
		ToroTestRoe(ToroTestsInit &params) : ToroTest(params) {			
			_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Roe;
			_solutionName = "Roe";
		};		

	}; //ToroTest1Roe

	//HLLC
	class ToroTestHLLC : public ToroTest {
	public:	
		//Constructor overload
		ToroTestHLLC(ToroTestsInit &params) : ToroTest(params) {			
			_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;
			_solutionName = "HLLC";
		};
	}; //ToroTest1HLLC

	//Godunov
	class ToroTestGodunov : public ToroTest {
	public:	
		//Constructor overload
		ToroTestGodunov(ToroTestsInit &params) : ToroTest(params) {			
			_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Godunov;
			_solutionName = "Godunov";
		};
	}; //ToroTest1Godunov

	//Constructor
	ToroTestComposite(ToroTestsInit& initSettings) : ToroTest(initSettings) {		
	};	

	//Test1 (automatic)
	void RunTest(int* argc, char** argv[]) {		
		//Initialize kernel
		_kernel = new Kernel;
		_kernel->Initilize(argc, argv);				

		//Prepare grid
		PrepareGrid();
		_kernel->SaveGrid("result.cgns");
		_kernel->VerboseOff();

		//Run all tests	
		_statuses.clear();		
		ToroTest* currentTest = nullptr;
		TestCaseResultInfo info;

		//Run Roe method for this test
		currentTest = new ToroTestRoe(_initSettings);
		std::cout << "Roe solver is running!\n";
		currentTest->RunTestWithKernel(_kernel);
		std::cout << "Roe solver finished!\n";

		//Check Roe test results
		info = currentTest->GetTestCaseResultInfo();
		_statuses.push_back(info.testStatus);
		if (info.testStatus == TestCaseResultInfo::TestStatusType::Passed) {
			std::cout<<"PASSED\n";
		} else {
			std::cout<<"FAILED\n";
		};		

		//Run HLLC method for this test
		currentTest = new ToroTestHLLC(_initSettings);
		std::cout << "HLLC solver is running!\n";
		currentTest->RunTestWithKernel(_kernel);
		std::cout << "HLLC solver finished!\n";

		//Check HLLC test results
		info = currentTest->GetTestCaseResultInfo();
		_statuses.push_back(info.testStatus);
		if (info.testStatus == TestCaseResultInfo::TestStatusType::Passed) {
			std::cout<<"PASSED\n";
		} else {
			std::cout<<"FAILED\n";
		};		

		//Run Godunov method for this test
		currentTest = new ToroTestGodunov(_initSettings);
		std::cout << "Godunov solver is running!\n";
		currentTest->RunTestWithKernel(_kernel);
		std::cout << "Godunov solver finished!\n";

		//Check Roe test results
		info = currentTest->GetTestCaseResultInfo();
		_statuses.push_back(info.testStatus);
		if (info.testStatus == TestCaseResultInfo::TestStatusType::Passed) {
			std::cout<<"PASSED\n";
		} else {
			std::cout<<"FAILED\n";
		};				
		
		//Save analitical solution
		WriteExactSolution("result.cgns", "Analytical");
		std::cout<<"Analytical solution written.\n";		

		//Finilazie kernel
		_kernel->Finalize();
	};

	//Get results of test run 
	virtual TestCaseResultInfo GetTestCaseResultInfo() override 
	{		
		//Compute composite test results
		TestCaseResultInfo result;
		result.testStatus = TestCaseResultInfo::TestStatusType::Passed;
		for (TestCaseResultInfo::TestStatusType status : _statuses) {
			
			//If any of tests failed
			if (status != TestCaseResultInfo::TestStatusType::Passed) {
				result.testStatus = TestCaseResultInfo::TestStatusType::Failed;
			};
		};
		return result;
	};

};

}; //namespace


#endif
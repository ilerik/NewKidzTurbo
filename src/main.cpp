#include "tests.h"

//Main program
int main(int argc, char *argv[]) {
	//RunSAFlatPlate();
	//RunGAWCalculation();
	//RunPoiseuilleTest();
	//RunRiemannProblemTest(ToroTestInitDistribution1, 100, 0.15, 2);
	//RunRiemannProblemTest(SODTestInitDistribution, 400, 0.2, 2);
	//RunShearFlowTest();
	//RunBlasiusTest();		//test has bad grid (not from ANSYS)
	//RunIncompressibleBlasius();
	//RunBumpFlow();
	//CellGradientTest();
	//RunSteadyShock();
	//RunVoronka();
	//ConvertGrid("SimpleCircle.dat");
	//RunBlasiusFlowAnsysGridTest();
	RunBlasiusFlowAnsysGridTest(150.0, 40);

	return 0;

};
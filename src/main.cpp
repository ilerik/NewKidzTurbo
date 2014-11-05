#include "tests.h"
#include "releigh-taylor_comp.h"

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
	//RunBiffFlatPlane();
	//RunGLSFlatPlane();
	//RunBumpFlow();
	//CellGradientTest();
	//RunSteadyShock();
	//RunVoronka();
	//ConvertGrid("SimpleCircle.dat");
	//RunBlasiusFlowAnsysGridTest();
	//RunBlasiusFlowAnsysGridTest(150.0, 40);
	//RunGLSFlatPlane();
	//RunReleighTaylor2D(30, 90);
	//RunReleighTaylor3D(50, 50, 50);
	//RunCollisionTest2DF(100, 100);
	RunCollisionTest2D(100, 100, 1.0e-6, 2.0e-5);
	//RunCollisionTest2DGausseAcceleration(20, 200, 1.0e-7, 1.01e-6);
	//RunReleighTaylor2Dvar2(100, 100);
	//RunReleighTaylor3Dvar1(50, 50, 50);

	//RunFePuCollision1DTest(1000, 5.0e-5, 0.25e-6);

	return 0;

};
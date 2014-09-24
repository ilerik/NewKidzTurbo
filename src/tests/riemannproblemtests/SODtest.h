#ifndef TURBO_TESTS_RIEMANNPROBLEMTEST_SODTEST
#define TURBO_TESTS_RIEMANNPROBLEMTEST_SODTEST

#include "testcase.h"

//SOD Test new version
class SodTestInitialConditions : public InitialConditions::InitialConditions
{
public:
	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;
			//http://en.wikipedia.org/wiki/Sod_shock_tube
			//Grid sizes
			double Lx = 0.02; // SI 2 cm				

			//Other velocities
			double v = 0;
			double w = 0;

			//Left state
			double roL = 1.0;
			double PL = 1.0;
			double uL = 0;
			
			//Right state
			double roR = 0.125;
			double PR = 0.1;
			double uR = 0;

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
				
			//Values
			double u = 0;
			double ro = 0;
			double roE = 0;
			if (x <= 0.5) {
				ro = roL;
				u = uL;
				roE = PL/(_gasModel->Gamma - 1.0);
			} else {
				ro = roR;
				u = uR;
				roE = PR/(_gasModel->Gamma - 1.0);
			};
			
			//Convert to conservative variables
			initValues.resize(_gasModel->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0);

			return initValues;		
	};
};

class SodTest : public TestCase
{
public:
	//Constructor
	//void SodTest(int _arg, char *_argv[]) : TestCase(int _arg, char *_argv[]) {};
	SodTest(int _argc, char *_argv[]) {
		argc = _argc;
		argv = _argv;
	};

	//Destructor
	~SodTest() {};

	//run test function
	bool RunTest() {		
		//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\Mixer\\Mixer.cgns");	
		//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns");	
		//_kernel.ReadInitialConditions("FlowSolution.E:1");
		//_kernel.ReadConfiguration("");

		//Test case description
		TestInfo = "Sod schock tube classical test case.";

		//Test case configuration
		Configuration conf;
		conf.InputCGNSFile = "";
		conf.OutputCGNSFile = "result.cgns";

		//Ideal gas model parameters
		conf.GasModel = CaloricallyPerfect;
		conf.IdealGasConstant = 8.3144621;
		conf.SpecificHeatRatio = 1.4;
		conf.SpecificHeatVolume = 1006.43 / 1.4;
		conf.SpecificHeatPressure = 1006.43;

		//Solver settings
		conf.CFL = 0.5;
		conf.RungeKuttaOrder = 1;

		//Simulation parameters
		conf.MaxIteration = 2000;
		conf.MaxTime = 0.2;

		//Boundary conditions		
		conf.BoundaryConditions["top"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		conf.BoundaryConditions["bottom"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		conf.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		conf.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCSymmetryPlane;

		//Initialize kernel
		Kernel _kernel;	
		_kernel.Initilize(&argc, &argv);

		//Generate grid
		const int nCells = 100;
		const double xMin = 0.0;
		const double xMax = 1.0;
		Grid _grid = GenGrid2D(_kernel.getParallelHelper(), nCells, 1, xMin, xMax, 0.0, 1.0, 1.0, 1.0);
		_kernel.BindGrid(_grid);

		//Set configuration
		_kernel.SetConfiguration(conf);

		//Init calculation
		_kernel.InitCalculation();

		//Set initial conditions 
		SodTestInitialConditions ic;
		_kernel.GenerateInitialConditions(ic);	

		//Run calculation
		_kernel.RunCalculation();

		//Finalize
		_kernel.FinalizeCalculation();	

		//Check result and generate report
		//_kernel.SaveGrid("result.cgns");
		_kernel.SaveSolution("result.cgns", "Solution");
		//_kernel.SaveSolution();
		//_kernel.SaveGrid("grid.cgns");
		_kernel.Finalize();	

		return true;
	};
};




#endif
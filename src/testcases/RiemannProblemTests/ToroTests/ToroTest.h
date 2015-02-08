#ifndef TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest
#define TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest

#include "TestCase.h"
#include "gengrid1D.h"


namespace ToroTests {

//Base class for all test from Toro book (the paragraph 4.3.3 pp. 129 - 133)
class ToroTest : public TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
	RiemannSolverConfiguration::RiemannSolverType _riemannSolverType;		//type of Riemann Solver Problem
	GasModel::ConservativeVariables U_an;			//analitical solution

public:
	//Test parameters
	int nCells;
	double Lx;					//solve problem in [0; Lx] segment
	double discontinuity_pos;		//posiyion of initial discontinuity
	double DiscontunityPos;		//position of discontinuity
	double TimeMax;				//time of solution

	//Left state density, pressure and velocity
	double roL;
	double pL;
	double uL;

	//Right state density, pressure and velocity
	double roR;
	double pR;
	double uR;

	//set parameters for test 1
	void SetParams() {
		//Test constant's
		nCells = 1000;
		Lx = 1.0;
		discontinuity_pos = 0.5;
		TimeMax = 0.2;

		//Left state density, pressure and velocity
		roL = 1.0;
		pL = 1.0;
		uL = 0.0;

		//Right state density, pressure and velocity
		roR = 0.1;
		pR = 0.125;
		uR = 0.0;
	};

	//Prepare computational grid
	void PrepareGrid() {		
		_grid = GenGrid1D(_kernel->getParallelHelper(), nCells, 0.0, Lx, false);
		_kernel->BindGrid(&_grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Godunov;

		//Availible gas models
		_configuration.AddGasModel("Air");			

		//Air (ideal gas)
		_configuration.GasModelsConfiguration["Air"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatPressure", 1006.43);			

		//Boundary conditions				
		_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["left"].MaterialName = "Air";

		_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["right"].MaterialName = "Air";

		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.1;
		_configuration.RungeKuttaOrder = 1;		

		//ALE settings
		_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = 0.1*TimeMax;			

		_kernel->BindConfiguration(_configuration);	

		return _configuration;
	};

	//Main interface function for running test case code
	void RunTestWithKernel(Kernel* kernel) {		
		_kernel = kernel; //Save reference to kernel
		
		//Get grid
		PrepareGrid();

		//Set parameters
		PrepareConfiguration();
		
		//Initialize calculation
		_kernel->InitCalculation();

		//Initial conditions
		_kernel->GenerateInitialConditions(new TestCaseInitialConditions(uL, roL, pL, 0, uR, roR, pR, 0, Lx, discontinuity_pos));	

		//Run computational cycle
		_kernel->RunCalculation();

		//Finilize
		_kernel->FinalizeCalculation();

		//Check results
		//Output result
		_kernel->SaveGrid("result.cgns");
		_kernel->SaveSolution("result.cgns", "Solution");
	};

	//Run test with program arguments
	void RunTest(int* argc, char** argv[]) {
		_kernel = new Kernel();
		_kernel->Initilize(argc, argv);

		//Call main function
		RunTestWithKernel(_kernel);

		_kernel->Finalize();
	};

	//Specify initial conditions
	//Sod's shock tube
	class TestCaseInitialConditions : public InitialConditions::InitialConditions
	{
	public:
		//Left state density, pressure, velocity and material index
		double _uLeft;
		double _roLeft;
		double _pLeft;
		int _nmatLeft;

		//Right state density, pressure, velocity and material index
		double _uRight;
		double _roRight;
		double _pRight;
		int _nmatRight;

		//length of domain and initial discontinuity position
		double _lx;
		double _x0;

		//constructor
		TestCaseInitialConditions(double uLeft, double roLeft, double pLeft, int nmatLeft, double uRight, double roRight, double pRight, int nmatRight, double lx, double x0) :
		_uLeft(uLeft),
		_roLeft(roLeft),
		_pLeft(pLeft),
		_nmatLeft(nmatLeft),
		_uRight(uRight),
		_roRight(roRight),
		_pRight(pRight),
		_nmatRight(nmatRight),
		_lx(lx),
		_x0(x0)
		{ };

		virtual int getInitialGasModelIndex(const Cell& cell) {
			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
		
			return 0; 
		};

		virtual std::vector<double> getInitialValues(const Cell& cell) {
			int nmat = getInitialGasModelIndex(cell); //get material index

			std::vector<double> initValues;				
			//Other velocities
			double v = 0;
			double w = 0;

			//Internal energy
			double e = 0;		

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
			
			//Values
			double u = 0;
			double ro = 0;
			double roE = 0;				
			if (x <= _x0) {
				ro = _roLeft;
				u = _uLeft;			
				e = _pLeft / (0.4 * ro);
			} else {
				ro = _roRight;
				u = _uRight;						
				e = _pRight / (0.4 * ro);
			};
			
			//Convert to conservative variables						
			roE = ro*(e + (u*u + v*v + w*w) / 2.0);
			initValues.resize(_gasModels[nmat]->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;

			return initValues;
		};
	}; //Initial conditions

}; //TestCase

}; //namespace


#endif
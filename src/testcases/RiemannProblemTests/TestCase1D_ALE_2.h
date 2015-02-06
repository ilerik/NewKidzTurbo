#ifndef TURBO_TestCases_TestCases1D_TestCaseALE2
#define TURBO_TestCases_TestCases1D_TestCaseALE2

#include "TestCase.h"
#include "gengrid1D.h"

namespace TestCases1D {

//Numerical test for ALE method from Luo et al 2004
//The well-known Sod shocktube problem is considered in this test case, whose solution contains 
//simultaneously a shock wave, a contact discontinuity, and an expansion fan.
class TestCaseALE2 : public TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
public:
	//Test parameters
	static const int nCells;
	static const double Lx;
	static const double TimeMax;

	//Left state density, pressure and velocity
	static const double roL;
	static const double pL;
	static const double uL;
	static const double gammaL;

	//Right state density, pressure and velocity
	static const double roR;
	static const double pR;
	static const double uR;
	static const double gammaR;

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
		_configuration.RiemannSolverConfiguration.RiemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;

		//Availible gas models		

		//Air (ideal gas)
		_configuration.AddGasModel("First");
		_configuration.GasModelsConfiguration["First"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["First"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["First"].SetPropertyValue("SpecificHeatRatio", gammaL);
		_configuration.GasModelsConfiguration["First"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["First"].SetPropertyValue("SpecificHeatPressure", 1006.43);

		//Air (ideal gas)
		_configuration.AddGasModel("Second");
		_configuration.GasModelsConfiguration["Second"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Second"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Second"].SetPropertyValue("SpecificHeatRatio", gammaR);
		_configuration.GasModelsConfiguration["Second"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Second"].SetPropertyValue("SpecificHeatPressure", 1006.43);

		//Boundary conditions				
		_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCOutflowSupersonic;
		_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::Fixed;
		//_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::FreeSurface;
		_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCOutflowSupersonic;
		_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::Fixed;
		//_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::FreeSurface;
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.5;
		_configuration.RungeKuttaOrder = 1;		

		//ALE settings
		//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";		
		_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";		
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = TimeMax / 10;			

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
		_kernel->GenerateInitialConditions(new TestCaseInitialConditions());	

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
	private:		
	public:			

		virtual int getInitialGasModelIndex(const Cell& cell) {
			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Split in 2 halfs
			if (x <= Lx/2.0) {
				return 0; //First
			} else {
				return 1; //Second
			};				
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
			if (x <= Lx/2.0) {
				ro = roL;
				u = uL;			
				e = pL / ((gammaL - 1.0) * ro);
			} else {
				ro = roR;
				u = uR;						
				e = pR / ((gammaR - 1.0) * ro);
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

//Test constant's
const int TestCaseALE2::nCells = 100;
const double TestCaseALE2::Lx = 1.0;
const double TestCaseALE2::TimeMax = 0.2;

//Left state density, pressure and velocity
const double TestCaseALE2::roL = 1.0;
const double TestCaseALE2::pL = 1.0;
const double TestCaseALE2::uL = 0.0;
const double TestCaseALE2::gammaL = 1.4;

//Right state density, pressure and velocity
const double TestCaseALE2::roR = 0.125;
const double TestCaseALE2::pR = 0.1;
const double TestCaseALE2::uR = 0.0;
const double TestCaseALE2::gammaR = 1.4;

}; //namespace


#endif
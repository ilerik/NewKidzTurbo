#ifndef TURBO_TestCases_MetalsCollision_MetalsCollision1DSteelVSSteel
#define TURBO_TestCases_MetalsCollision_MetalsCollision1DSteelVSSteel

#include "TestCase.h"
#include "gengrid1D.h"

namespace TestCasesMetalsImpact {

//Numerical test for ALE method from Luo et al 2004
//The well-known Sod shocktube problem is considered in this test case, whose solution contains 
//simultaneously a shock wave, a contact discontinuity, and an expansion fan.
class TestCaseMetalsImpact_1D_SteelVSSteel: public TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
public:
	//Test parameters
	static const int nCells;
	static const double LLeft;
	static const double LRight;
	static const double TimeMax;

	//Left state density, material index for EOS5 and velocity
	static const double roSteel;
	static const int nmet;

	//Impact relative speed
	double uImpact;

	//Parametrized constructor
	TestCaseMetalsImpact_1D_SteelVSSteel(double _uImpact) {
		uImpact = _uImpact;
	};

	//Prepare computational grid
	void PrepareGrid() {		
		_grid = GenGrid1D(_kernel->getParallelHelper(), nCells, -2 * LLeft, 2 * LRight, false);
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

		//Stainless steel
		_configuration.AddGasModel("StainlessSteelLeft");
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].GasModelName = "LomonosovFortovGasModel";
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].SetPropertyValue("MaterialIndex", nmet);
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].SetPropertyValue("SpecificHeatVolume", 510); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
		_configuration.GasModelsConfiguration["StainlessSteelLeft"].SetPropertyValue("MeltingTemperature", 1800); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI 

		_configuration.AddGasModel("StainlessSteelRight");
		_configuration.GasModelsConfiguration["StainlessSteelRight"].GasModelName = "LomonosovFortovGasModel";
		_configuration.GasModelsConfiguration["StainlessSteelRight"].SetPropertyValue("MaterialIndex", nmet);
		_configuration.GasModelsConfiguration["StainlessSteelRight"].SetPropertyValue("SpecificHeatVolume", 510); //From http://www.diracdelta.co.uk/science/source/s/p/specific%20heat%20capacity/source.html#.VMr8MP6sXQI
		_configuration.GasModelsConfiguration["StainlessSteelRight"].SetPropertyValue("MeltingTemperature", 1800); //From http://www.diracdelta.co.uk/science/source/m/e/melting%20point/source.html#.VMr8x_6sXQI 

		//Air (ideal gas)
		_configuration.AddGasModel("Air");
		_configuration.GasModelsConfiguration["Air"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatPressure", 1006.43);

		//Boundary conditions				
		_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::Fixed;
		//_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::FreeSurface;
		_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::Fixed;
		//_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::FreeSurface;
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.2;
		_configuration.RungeKuttaOrder = 1;		

		//ALE settings
		//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";		
		_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";		
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = TimeMax / 100;			

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
		_kernel->GenerateInitialConditions(new TestCaseInitialConditions(uImpact));	

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
		double uImpact;

		TestCaseInitialConditions(double _uImpact) {
			uImpact = _uImpact;
		};

		virtual int getInitialGasModelIndex(const Cell& cell) {
			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Outside is air
			if (x <= - LLeft) return 2;
			if (x >= LRight) return 2;

			//Split in 2 halfs
			if (x <= 0) {
				return 0; //First
			} else {
				return 1; //Second
			};		

			return 0;
		};

		virtual std::vector<double> getInitialValues(const Cell& cell) {
			int nmat = getInitialGasModelIndex(cell); //get material index

			std::vector<double> initValues;		
			double uCenterMass = 0; //uImpact * (LLeft * roSteel) / (LLeft * roSteel + LRight * roSteel);
			double uL = uImpact - uCenterMass;
			double uR = -uCenterMass;

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
			double roE = 0;		
			double ro = roSteel;
			if (x <= 0) {
				u = uL;			
			} else {
				u = uR;						
			};

			//Outside is air
			double atm = 1e5;
			double gamma = 1.4;
			if ((x <= - LLeft) || (x >= LRight)) {
				ro = 1.225;
				u = 0;
				e = atm / ((gamma - 1.0) * ro);
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
const int TestCaseMetalsImpact_1D_SteelVSSteel::nCells = 1000;
const double TestCaseMetalsImpact_1D_SteelVSSteel::LLeft = 15e-2; // 15 cm;
const double TestCaseMetalsImpact_1D_SteelVSSteel::LRight = 15e-2; // 15 cm;
const double TestCaseMetalsImpact_1D_SteelVSSteel::TimeMax = 100e-6; // 100 mks

//Density and material index for EOS5 and velocity
const double TestCaseMetalsImpact_1D_SteelVSSteel::roSteel = 1000 * 1.0 / 0.127; //SI	for stainless steel;
const int TestCaseMetalsImpact_1D_SteelVSSteel::nmet = 0;

}; //namespace


#endif
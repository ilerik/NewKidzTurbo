#ifndef TURBO_TestCases_RMInstabilityTests_TestCase1
#define TURBO_TestCases_RMInstabilityTests_TestCase1

#include "TestCase.h"
#include "gengrid2D.h"

namespace TestCasesRMI {

//Base class for all automated tests
class TestCase1 : public TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
public:
	//Test static constant's 
	static const int nCellsX;
	static const int nCellsY;
	static const double xMax;
	static const double xMin;
	static const double yMax;
	static const double yMin;
	static const double TimeMax;

	//Density, pressure and gamma for mixing fluids
	static const double atm;
	static const double roHeavy;
	static const double roLight;
	static const double pressure;
	static const double gamma;

	//Initial pertrubation parameters
	static const int ModesNumber;
	static const double a0;
	static const double lambdaX;

	//Shock wave parameters (in lighter fluid)
	static const double SoundSpeedLightFluid;
	static const double MachNumber;
	static const double pShock;
	static const double roShock;
	static const double uShock;

	//Prepare computational grid
	void PrepareGrid() {		
		_grid = GenGrid2D(_kernel->getParallelHelper(), nCellsX, nCellsY, xMin, xMax, yMin, yMax, 1.0, 1.0, true, false);
		_kernel->BindGrid(&_grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.RiemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Roe;

		//Availible gas models
		_configuration.AddGasModel("Air");			

		//Air (ideal gas)
		_configuration.GasModelsConfiguration["Air"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", gamma);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatPressure", 1006.43);			

		//Boundary conditions				
		_configuration.BoundaryConditions["top"].BoundaryConditionType = BCType_t::BCOutflowSupersonic;
		_configuration.BoundaryConditions["bottom"].BoundaryConditionType = BCType_t::BCOutflowSupersonic;
		
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
		_configuration.SaveSolutionSnapshotTime = TimeMax / 50.0;			

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
		
			return 0; 
		};

		virtual std::vector<double> getInitialValues(const Cell& cell) {
			int nmat = getInitialGasModelIndex(cell); //get material index

			std::vector<double> initValues;				
			//Other velocities
			double v = 0;
			double w = 0;
			double u = 0;

			//Internal energy
			double e = 0;		

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Position of inteface (y = 0)
			double yInterface = 0 + a0 * std::cos((2*PI / lambdaX) * x);

			//Position of shock wave (y = 2*a0)
			double yShockWave = std::abs(a0);
				
			//Values			
			double ro = 0;
			double roE = 0;
			if (y <= yInterface) {
				//Light fluid
				ro = roLight;
				v = 0;			
				e = pressure / ((gamma-1) * ro);
			};
			if ((y > yInterface) && (y <= yShockWave)) {
				//Heavy fluid at rest
				ro = roHeavy;
				v = 0;						
				e = pressure / ((gamma-1) * ro);
			};
			if (y > yShockWave) {
				//Shock wave in heavy fluid
				ro = roShock;
				v = -uShock;
				e = pShock / ((gamma-1) * ro);
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
const int TestCase1::nCellsX = 200;
const int TestCase1::nCellsY = 800;
const double TestCase1::xMax = TestCase1::ModesNumber * (TestCase1::lambdaX * 0.5);
const double TestCase1::xMin = TestCase1::ModesNumber * (-TestCase1::lambdaX * 0.5);
const double TestCase1::yMax = 3 * 1e-2; //[cm]
const double TestCase1::yMin = -4 * 1e-2; //[cm]
const double TestCase1::TimeMax = 0.05 * 1e-3; // [ms]

//Density, pressure and gamma for mixing fluids
const double TestCase1::atm = 101.325e3; // std atmosphere pressure [Pa]
const double TestCase1::roHeavy = 2.74; //xenon [kg/m^3]
const double TestCase1::roLight = 0.084; //helium [kg/m^3]
const double TestCase1::pressure = 0.5*atm; //[atm]
const double TestCase1::gamma = 1.67;

//Initial pertrubation parameters
const int TestCase1::ModesNumber = 1; //Number of modes
const double TestCase1::a0 = -0.25 * 1e-2; //Pertrubation amplitude [cm] 
const double TestCase1::lambdaX = 0.8 * 1e-2; //Wave number [cm]

//Shock wave parameters (in lighter fluid)
const double TestCase1::MachNumber = 2.5;
const double TestCase1::SoundSpeedLightFluid = std::sqrt(TestCase1::gamma * TestCase1::pressure / TestCase1::roLight);
const double TestCase1::pShock = 1.34848*atm; //[atm]
const double TestCase1::roShock = 20.7347; //[kg/m^3]
const double TestCase1::uShock = 262.259; //[m/s]
//const double TestCase1::pShock = 1.34848*atm; //[atm]
//const double TestCase1::roShock = 0.635663; //[kg/m^3]
//const double TestCase1::uShock = 1497.84; //[m/s]

}; //namespace


#endif
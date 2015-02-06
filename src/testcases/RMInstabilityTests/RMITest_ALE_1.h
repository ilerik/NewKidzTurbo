#ifndef TURBO_TestCases_RMInstabilityTests_TestCaseALE1
#define TURBO_TestCases_RMInstabilityTests_TestCaseALE1

#include "TestCase.h"
#include "gengrid2D.h"

namespace TestCasesRMI {

//Base class for all automated tests
class TestCaseALE1 : public TestCase {
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
		_kernel->PartitionGrid(_grid);
		_kernel->GenerateGridGeometry(_grid);

		MeshMovement moveHelper;
		int nDeformationSteps = 2;

		//Modify grid for initial disturbances
		double A = 1e-4;
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-A,A);
		std::vector<int> nodes;
		std::vector<Vector> displacements;

		//Determine coordinate of nodes closest to shock front
		double yInterfacePosition = 0;
		double yInterfaceNodes = yMin;
		double yShockWave = 1.5 * std::abs(a0);
		double yShockWaveNodes = 0;
		for (Node& n : _grid.localNodes) {
			if (std::abs(n.P.y - yShockWave) < std::abs(yShockWave - yShockWaveNodes)) yShockWaveNodes = n.P.y;
			if (std::abs(n.P.y - yInterfacePosition) < std::abs(yInterfaceNodes - yInterfacePosition)) yInterfaceNodes = n.P.y;
		};

		for (Node& n : _grid.localNodes) {
			//Unmovable borders
			//if ((n.P.x == xMin) || (n.P.x == xMax) || (n.P.y== yMin) || (n.P.y == yMax)) {
			if ((n.P.y== yMin) || (std::abs(n.P.y - yMax) < 1e-10)) {
				nodes.push_back(n.GlobalIndex);
				displacements.push_back(Vector(0,0,0));
			} else {
				//Random distortion
				/*if (n.P.x == 0) {
					double delta = distribution(generator);
					Vector dr(delta, 0, 0);
					nodes.push_back(n.GlobalIndex);
					displacements.push_back(dr);
				};*/
				// Harmonic
				if ((std::abs(n.P.y - yInterfaceNodes) < 1e-10)) {
					double yInterface = yInterfacePosition + a0 * std::cos((2*PI / lambdaX) * n.P.x);
					double delta = yInterface;
					Vector dr(0, delta, 0);
					nodes.push_back(n.GlobalIndex);
					displacements.push_back(dr / nDeformationSteps);
				};
				
			};			
		};

		for (int i = 0; i<nDeformationSteps; i++) {
			moveHelper.MoveNodes(_grid, nodes, displacements);
		};

		//Fixed straight shock front
		for (Node& n : _grid.localNodes) {
			//Unmovable borders
			//if ((n.P.x == xMin) || (n.P.x == xMax) || (n.P.y== yMin) || (n.P.y == yMax)) {
			if ((n.P.y== yMin) || (std::abs(n.P.y - yMax) < 1e-10)) {
				nodes.push_back(n.GlobalIndex);
				displacements.push_back(Vector(0,0,0));
			} else {
				if ((std::abs(n.P.y - yShockWaveNodes) < 1e-10)) {
					double delta = yShockWaveNodes - yShockWave;
					Vector dr(0, delta, 0);
					nodes.push_back(n.GlobalIndex);
					displacements.push_back(dr / nDeformationSteps);
				};
			};
		};

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
			
		//He (ideal gas)
		_configuration.AddGasModel("He");
		_configuration.GasModelsConfiguration["He"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["He"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["He"].SetPropertyValue("SpecificHeatRatio", gamma);
		_configuration.GasModelsConfiguration["He"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["He"].SetPropertyValue("SpecificHeatPressure", 1006.43);		

		//Xe (ideal gas)
		_configuration.AddGasModel("Xe");
		_configuration.GasModelsConfiguration["Xe"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Xe"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Xe"].SetPropertyValue("SpecificHeatRatio", gamma);
		_configuration.GasModelsConfiguration["Xe"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Xe"].SetPropertyValue("SpecificHeatPressure", 1006.43);	

		//Boundary conditions				
		_configuration.BoundaryConditions["top"].BoundaryConditionType = BCType_t::BCOutflowSupersonic;
		_configuration.BoundaryConditions["top"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["top"].MaterialName = "He";
		_configuration.BoundaryConditions["bottom"].BoundaryConditionType = BCType_t::BCOutflowSupersonic;
		_configuration.BoundaryConditions["bottom"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["bottom"].MaterialName = "Xe";
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.5;
		_configuration.RungeKuttaOrder = 4;		

		//ALE settings
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";
		_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";

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

		_kernel->SaveGrid("grid.cgns");

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

			//Position of inteface (y = 0)
			double yInterface = 0 + a0 * std::cos((2*PI / lambdaX) * x);

			if (y <= yInterface) {
				return 1;
			} else {
				return 0;
			};
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
			double yShockWave = 1.5 * std::abs(a0);
				
			//Values			
			double ro = 0;
			double roE = 0;
			if (y <= yInterface) {
				//Heavy fluid
				ro = roHeavy;
				v = 0;			
				e = pressure / ((gamma-1) * ro);
			};
			if ((y > yInterface) && (y <= yShockWave)) {
				//Light fluid at rest
				ro = roLight;
				v = 0;						
				e = pressure / ((gamma-1) * ro);
			};
			if (y > yShockWave) {
				//Shock wave in light fluid
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
const int TestCaseALE1::nCellsX = 40;
const int TestCaseALE1::nCellsY = 175;
const double TestCaseALE1::xMax = TestCaseALE1::ModesNumber * (TestCaseALE1::lambdaX * 0.5);
const double TestCaseALE1::xMin = TestCaseALE1::ModesNumber * (-TestCaseALE1::lambdaX * 0.5);
const double TestCaseALE1::yMax = 3 * 1e-2; //[cm]
const double TestCaseALE1::yMin = -4 * 1e-2; //[cm]
const double TestCaseALE1::TimeMax = 0.05 * 1e-3; // [ms]

//Density, pressure and gamma for mixing fluids
const double TestCaseALE1::atm = 101.325e3; // std atmosphere pressure [Pa]
const double TestCaseALE1::roHeavy = 2.74; //xenon [kg/m^3]
const double TestCaseALE1::roLight = 0.084; //helium [kg/m^3]
const double TestCaseALE1::pressure = 0.5*atm; //[atm]
const double TestCaseALE1::gamma = 1.67;

//Initial pertrubation parameters
const int TestCaseALE1::ModesNumber = 1; //Number of modes
const double TestCaseALE1::a0 = 0.25 * 1e-2; //Pertrubation amplitude [cm] 
const double TestCaseALE1::lambdaX = 0.8 * 1e-2; //Wave number [cm]

//Shock wave parameters (in lighter fluid)
const double TestCaseALE1::MachNumber = 2.5;
const double TestCaseALE1::SoundSpeedLightFluid = std::sqrt(TestCaseALE1::gamma * TestCaseALE1::pressure / TestCaseALE1::roLight);
const double TestCaseALE1::pShock = 1.34848*atm; //[atm]
const double TestCaseALE1::roShock = 0.635663; //[kg/m^3]
const double TestCaseALE1::uShock = 1497.84; //[m/s]

}; //namespace


#endif
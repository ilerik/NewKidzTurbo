#ifndef TURBO_TestCases_MetalsCollision_MetalsImpact1DTestCase
#define TURBO_TestCases_MetalsCollision_MetalsImpact1DTestCase

#include "MetalsImpactTestCase.h"
#include "DSUClusteredSet.h"
#include "gengrid1D.h"

namespace TestCasesMetalsImpact {

//Base class for 1D impact test setting
class MetalsImpact1DTestCase: public MetalsImpactTestCase {	
protected:
	//Grid parameters
	int _nCells;

	//Pertrubation position
	double _pFunction(Vector r) {
		return r.x;
	};

public:
	//Parametrized constructor
	MetalsImpact1DTestCase(int nCells, int nSnapshots, double widthLeft, double widthRight, double TimeMax, MetalType metalLeft, double uLeft, MetalType metalRight, double uRight) {
		std::function<double(Vector)> ptr = std::bind(&TestCasesMetalsImpact::MetalsImpact1DTestCase::_pFunction, this, std::placeholders::_1);
		this->MetalsImpactTestCase::MetalsImpactTestCase(ptr, nSnapshots, widthLeft, widthRight, TimeMax, metalLeft, uLeft, metalRight, uRight);
		_nCells = nCells;
	};

	//Prepare computational grid
	void PrepareGrid() override  {		
		_grid = GenGrid1D(_kernel->getParallelHelper(), _nCells, -_widthLeft, _widthRight, false);
		_kernel->PartitionGrid(_grid);
		_kernel->GenerateGridGeometry(_grid);

		//Determine nodes that a closest to zero
		int interfaceNodeIndex;
		double interfaceNodeX = std::max(_widthLeft, _widthRight); 
		for (Node& node : _grid.localNodes) {
			if (std::abs(node.P.x) <  interfaceNodeX) {
				interfaceNodeX = std::abs(node.P.x);
				interfaceNodeIndex = node.GlobalIndex;
			};
		};

		//Compute displacements
		std::vector<int> nodes;
		std::vector<Vector> displacements;
		nodes.push_back(interfaceNodeIndex);
		displacements.push_back(Vector(-interfaceNodeX, 0, 0));

		//Move mesh
		MeshMovement moveHelper;
		moveHelper.meshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::Linear1D;
		moveHelper.MoveNodes(_grid, nodes, displacements);
		
		_kernel->BindGrid(&_grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() override  {
		//Hardcode configuration for now
		_configuration.WorkingDirectory = "./results";
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;

		//Availible gas models		
		//Left metal
		_nmatLeft = 0;
		_configuration.AddGasModel("MetalLeft");
		_configuration.GasModelsConfiguration["MetalLeft"] = GetGasModelConfiguration("LomonosovFortovGasModel", _metalLeft, _roLeft);

		//Right metal
		_nmatRight = 1;
		_configuration.AddGasModel("MetalRight");
		_configuration.GasModelsConfiguration["MetalRight"] = GetGasModelConfiguration("LomonosovFortovGasModel", _metalRight, _roRight);

		//Boundary conditions
		//Air (ideal gas)
		double gamma = 1.4;
		_nmatBoundary = 2;
		_configuration.AddGasModel("Air");
		_configuration.GasModelsConfiguration["Air"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatPressure", 1006.43);

		//Boundary conditions
		_configuration.BoundaryConditions["left"] = GetBoundaryConditionConfiguration("MetalLeft", BoundaryConditionType::Wall);
		_configuration.BoundaryConditions["right"] = GetBoundaryConditionConfiguration("MetalRight", BoundaryConditionType::FreeSurface);		
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.5;
		_configuration.RungeKuttaOrder = 4;		

		//ALE settings
		_configuration.ALEConfiguration.MeshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::Linear1D;
		//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";		
		_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";		
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = _TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = _TimeMax / _nSnapshots;

		_kernel->VerboseOn();
		_kernel->BindConfiguration(_configuration);	

		return _configuration;
	};

	//Get results of test run 
	virtual TestCaseResultInfo GetTestCaseResultInfo() override { 
		throw new Exception("Not implemented");
		return TestCaseResultInfo();
	};
	
	//Prepare kernel
	void PrepareKernel(Kernel* kernel) {				
		_kernel = kernel; //Save reference to kernel		
		
		//Get grid
		PrepareGrid();

		//Set parameters
		PrepareConfiguration();
	};

	//Main interface function for running test case code
	virtual void RunTestWithKernel(Kernel* kernel) override {
		//Set history logger
		_kernel->setStepHistoryLogger(new TestCaseHistoryLogger());
		
		//Initialize calculation
		_kernel->InitCalculation();

		//Initial conditions
		_kernel->GenerateInitialConditions(new TestCaseInitialConditions( getPertrubationPositionFunction,
			_uLeft, 
			_roLeft,
			_nmatLeft,
			_uRight,
			_roRight,
			_nmatRight,
			_uBoundary,
			_roBoundary,
			_nmatBoundary
			));	

		//Run computational cycle
		_kernel->RunCalculation();

		//Finilize
		_kernel->FinalizeCalculation();

		//Output result
		_kernel->SaveGrid("result.cgns");
		_kernel->SaveSolution("result.cgns", "Solution");
	};

	//Run test with program arguments
	virtual void RunTest(int* argc, char** argv[]) {
		_kernel = new Kernel();
		_kernel->Initilize(argc, argv);
		
		//Prepare grid and configuration
		MetalsImpact1DTestCase::PrepareGrid();
		MetalsImpact1DTestCase::PrepareConfiguration();

		//Call main function
		MetalsImpact1DTestCase::RunTestWithKernel(_kernel);

		_kernel->Finalize();
	}; 

}; //TestCase

}; //namespace


#endif
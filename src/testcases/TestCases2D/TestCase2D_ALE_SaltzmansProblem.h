#ifndef TURBO_TestCases_MetalsCollision_MetalsImpact2DTestCase
#define TURBO_TestCases_MetalsCollision_MetalsImpact2DTestCase

#include "TestCase.h"
#include "DSUClusteredSet.h"
#include "gengrid2D.h"
#include <cmath>

namespace TestCases2D {

//Base class for 2D impact test setting
class TestCase2D_ALE_SaltzmanProblem : public TestCase {	
protected:
	//Grid parameters

	//Test parameters
	int _nSnapshots;
	double _TimeMax;

	//Pertrubation parameters

public:
	//Constructor
	TestCase2D_ALE_SaltzmanProblem() {
		_TimeMax = 0.6;
		_nSnapshots = 100;
	};

	//Prepare computational grid
	void PrepareGrid() {		
		double xMin = 0;
		double xMax = 1.0;
		double yMin = 0;
		double yMax = 0.1;
		_grid = GenGrid2D(_kernel->getParallelHelper(), 100, 10, xMin, xMax, yMin, yMax, 1.0, 1.0, false, false);
		_kernel->PartitionGrid(_grid);
		_kernel->GenerateGridGeometry(_grid);	

		//Compute displacements
		std::vector<int> nodes;
		std::vector<Vector> displacements;
		nodes.clear();
		displacements.clear();
		double dx = 0.01;
		double dy = 0.01;
		for (Node& node : _grid.localNodes) {
			int i = std::round(node.P.x / dx) + 1;
			int j = std::round(node.P.y / dy) + 1;
			double dY = (11 - j)*std::sin(PI * (i-1) / 100.0) * dy;
			Vector dr = Vector(0, dY, 0);
			nodes.push_back(node.GlobalIndex);
			displacements.push_back(dr);
		};

		//Move mesh
		MeshMovement moveHelper;
		moveHelper.meshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDWnoRotation;
		moveHelper.MoveNodes(_grid, nodes, displacements);
		
		_kernel->BindGrid(&_grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.WorkingDirectory = "./results";
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.RiemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;

		//Boundary conditions
		//Ideal gas
		int _nmatBoundary = 0;
		std::string boundaryMaterialName = "Gas";
		_configuration.AddGasModel(boundaryMaterialName);
		_configuration.GasModelsConfiguration[boundaryMaterialName] = GetBoundaryGasModelConfiguration();				
		_configuration.BoundaryConditions["left"] = GetBoundaryConditionConfiguration(boundaryMaterialName);
		_configuration.BoundaryConditions["right"] = GetBoundaryConditionConfiguration(boundaryMaterialName);		
		_configuration.BoundaryConditions["top"] = GetBoundaryConditionConfiguration(boundaryMaterialName);		
		_configuration.BoundaryConditions["bottom"] = GetBoundaryConditionConfiguration(boundaryMaterialName);				
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.5;
		_configuration.RungeKuttaOrder = 4;		

		//ALE settings
		_configuration.ALEConfiguration.MeshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDW;
		_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";		
		//_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";		
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = _TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = _TimeMax / _nSnapshots;

		_kernel->BindConfiguration(_configuration);	

		return _configuration;
	};	
	
	//Prepare kernel
	void PrepareKernel(Kernel* kernel) {				
		_kernel = kernel; //Save reference to kernel
		
		//Get grid
		PrepareGrid();

		//Set parameters
		PrepareConfiguration();
	};

	//Run kernel
	void RunKernel(Kernel *kernel) {
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
	virtual void RunTest(int* argc, char** argv[]) override {
		_kernel = new Kernel();
		_kernel->Initilize(argc, argv);
		
		//Prepare grid and configuration
		PrepareGrid();
		PrepareConfiguration();

		//Call main function
		RunKernel(_kernel);

		_kernel->Finalize();
	}; 

}; //TestCase

}; //namespace


#endif
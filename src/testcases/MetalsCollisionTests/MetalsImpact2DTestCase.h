#ifndef TURBO_TestCases_MetalsCollision_MetalsImpact2DTestCase
#define TURBO_TestCases_MetalsCollision_MetalsImpact2DTestCase

#include "MetalsImpactTestCase.h"
#include "DSUClusteredSet.h"
#include "gengrid2D.h"

namespace TestCasesMetalsImpact {

//Base class for 2D impact test setting
class MetalsImpact2DTestCase: public MetalsImpactTestCase {	
protected:
	//Grid parameters
	int _nCellsX;
	int _nCellsY;
	double _widthY;

	//Pertrubation parameters
	double _D;
	double _A;
	double _n;


	//Pertrubation parameters
	//double lamda = 1e-3; // 1mm
	//d

	//Pertrubation position
	double _pFunction(Vector r) {	
		if (std::abs(r.y) > _D / 2.0) return r.x;
		double x = _A * std::cos(2.0 * PI * r.y / (_n * _D));
		return r.x - x;
	};

public:
	//Parametrized constructor
	MetalsImpact2DTestCase(int nCellsX, int nCellsY, int nSnapshots, double widthY, double widthLeft, double widthRight, double TimeMax, MetalType metalLeft, double uLeft, MetalType metalRight, double uRight) {
		std::function<double(Vector)> ptr = std::bind(&TestCasesMetalsImpact::MetalsImpact2DTestCase::_pFunction, this, std::placeholders::_1);
		this->MetalsImpactTestCase::MetalsImpactTestCase(ptr, nSnapshots, widthLeft, widthRight, TimeMax, metalLeft, uLeft, metalRight, uRight);
		_widthY = widthY;
		_nCellsX = nCellsX;
		_nCellsY = nCellsY;

		//Pertrubation parameters
		_D = 0.2e-3;
		_A = 0.1e-3;
		_n = 1;
	};

	//Prepare computational grid
	void PrepareGrid() override  {		
		double xMin = -_widthLeft;
		double xMax = _widthRight;
		double yMin = - (_widthY / 2.0);
		double yMax = + (_widthY / 2.0);
		_grid = GenGrid2D(_kernel->getParallelHelper(), _nCellsX, _nCellsY, xMin, xMax, yMin, yMax, 1.0, 1.0, false, false);
		_kernel->PartitionGrid(_grid);
		_kernel->GenerateGridGeometry(_grid);

		//Determine nodes that a closest to zero
		std::map<double, int> interfaceNodeIndexes;
		std::map<double, double> interfaceNodeXs;
		for (Node& node : _grid.localNodes) {
			double x = node.P.x;
			double y = node.P.y;
			if (interfaceNodeXs.find(y) == std::end(interfaceNodeXs)) {
				interfaceNodeXs[y] = x;
				continue;
			};
			Vector leader = Vector(interfaceNodeXs[y], y, 0);
			if (std::abs(_pFunction(node.P)) <  std::abs(_pFunction(leader))) {
				interfaceNodeXs[y] = x;
				interfaceNodeIndexes[y] = node.GlobalIndex;
			};
		};		

		//Compute displacements
		std::vector<int> nodes;
		std::vector<Vector> displacements;
		nodes.clear();
		displacements.clear();
		for (std::pair<double, int> pair : interfaceNodeIndexes) {
			double y = pair.first;
			double dX = _pFunction(Vector(interfaceNodeXs[y], y, 0));						
			displacements.push_back(Vector(-dX, 0, 0));
			int nodeIndex = pair.second;
			nodes.push_back(nodeIndex);
		};

		//Move mesh
		MeshMovement moveHelper;
		moveHelper.meshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDWnoRotation;
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
		_configuration.GasModelsConfiguration["MetalLeft"] = GetGasModelConfiguration(_metalLeft, _roLeft);

		//Right metal
		_nmatRight = 1;
		_configuration.AddGasModel("MetalRight");
		_configuration.GasModelsConfiguration["MetalRight"] = GetGasModelConfiguration(_metalRight, _roRight);

		//Boundary conditions
		//Air (ideal gas)	
		_nmatBoundary = 2;
		std::string boundaryMaterialName = "Air";
		_configuration.AddGasModel(boundaryMaterialName);
		_configuration.GasModelsConfiguration[boundaryMaterialName] = GetBoundaryGasModelConfiguration();				
		_configuration.BoundaryConditions["left"] = GetBoundaryConditionConfiguration(boundaryMaterialName, BoundaryConditionType::Wall);
		_configuration.BoundaryConditions["right"] = GetBoundaryConditionConfiguration(boundaryMaterialName, BoundaryConditionType::FreeSurface);		
		_configuration.BoundaryConditions["top"] = GetBoundaryConditionConfiguration(boundaryMaterialName, BoundaryConditionType::FreeSurface);		
		_configuration.BoundaryConditions["bottom"] = GetBoundaryConditionConfiguration(boundaryMaterialName, BoundaryConditionType::FreeSurface);				
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.5;
		_configuration.RungeKuttaOrder = 4;		

		//ALE settings
		_configuration.ALEConfiguration.MeshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDWnoRotation;
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

	//Run kernel
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
		MetalsImpact2DTestCase::PrepareGrid();
		MetalsImpact2DTestCase::PrepareConfiguration();

		//Call main function
		MetalsImpact2DTestCase::RunTestWithKernel(_kernel);

		_kernel->Finalize();
	}; 

}; //TestCase

}; //namespace


#endif
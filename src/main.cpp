#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "kernel.h"
#include "gridgeneration.h"
#include <cmath>
#include <random>
#include "MeshMovement.h"

#include "cuda.h"
#include <cuda_runtime.h>

//Test cases
#include "test_list.h"

//Main program
int main(int argc, char *argv[]) {			
	RTInstabilityTests::TestSettings settings;

	//Pertrubation
	settings.geometrySettings.lambda = 0.5;
	settings.geometrySettings.A = 0.0; //0.01; //0.1e-2;	
	settings.geometrySettings.IsUserDefinedDisturbance = false;

	//Mesh
	int nWaves = 1;
	settings.geometrySettings.nCellX = 20;
	settings.geometrySettings.nCellY = 60;
	settings.geometrySettings.xMin = -0.5 * (nWaves * settings.geometrySettings.lambda);
	settings.geometrySettings.xMax = +0.5 * (nWaves * settings.geometrySettings.lambda);
	settings.geometrySettings.yMin = -0.75;
	settings.geometrySettings.yMax = +0.75;
	settings.geometrySettings.yInterface = 0.0;	

	//Boundaries
	settings.bcTop = RTInstabilityTests::BoundaryConditionType::Wall;
	settings.bcBottom = RTInstabilityTests::BoundaryConditionType::Wall;

	//Materials
	settings.materialSettings.gasModelHeavy = RTInstabilityTests::GasModelType::PerfectGas;
	settings.materialSettings.gasModelLight = RTInstabilityTests::GasModelType::PerfectGas;
	settings.materialSettings.roHeavy = 2.0;
	settings.materialSettings.roLight = 1.0;
	settings.Pinterface = 2.5;
	//settings.gravity = Vector(0.0, 0.0, 0.0);
	settings.gravity = Vector(0.0, -0.25, 0.0);

	//Method
	settings.methodSettings.meshMotionType = "ALEMaterialInterfaces";
	//settings.methodSettings.meshMotionType = "Lagrangian";			
	//settings.methodSettings.meshMotionType = "Eulerian";		
	settings.methodSettings.spatialReconstruction = SpatialDiscretisationType::ENO;

	//Run parameters
	int nSnapshots = 100;
	settings.MaxTime = 10.0;
	settings.SaveSolutionSnapshotTime = settings.MaxTime / nSnapshots;
	settings.MaxIteration = 1000000;
	settings.SaveSolutionSnapshotIterations = 0;

	//Create test object
	ParallelManager MPIManager(argc, argv);
	RTInstabilityTests::RTInstabilityTestCase test(MPIManager, settings);
	test.RunTest(&argc, &argv);

	//ParallelManager MPIManager(argc, argv);
	//GenGrid2DInterfacePertrubation(MPIManager, 

	return 0;
};

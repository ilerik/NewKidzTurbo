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
	settings.geometrySettings.lambda = 2.0e-2;
	settings.geometrySettings.A = 0.0; //0.1e-2;
	settings.geometrySettings.IsUserDefinedDisturbance = false;

	//Mesh
	int nWaves = 5;
	settings.geometrySettings.nCellX = 50;
	settings.geometrySettings.nCellY = 50;
	settings.geometrySettings.xMin = -0.5 * (nWaves * settings.geometrySettings.lambda);
	settings.geometrySettings.xMax = +0.5 * (nWaves * settings.geometrySettings.lambda);
	settings.geometrySettings.yMin =  0.0e-2;
	settings.geometrySettings.yMax = 10.0e-2;
	settings.geometrySettings.yInterface = 5.0e-2;

	settings.geometrySettings.lambda *= nWaves; // Make one wave only

	//Boundaries
	settings.bcTop = RTInstabilityTests::BoundaryConditionType::Wall;
	settings.bcBottom = RTInstabilityTests::BoundaryConditionType::Wall;

	//Materials
	settings.materialSettings.gasModelHeavy = RTInstabilityTests::GasModelType::PerfectGas;
	settings.materialSettings.gasModelLight = RTInstabilityTests::GasModelType::PerfectGas;
	settings.materialSettings.roHeavy = 2.0;
	settings.materialSettings.roLight = 1.0;
	settings.Pinterface = 1e10;
	settings.gravity = Vector(0.0, 0.0, 0.0);
	//settings.gravity = Vector(0.0, -1.0, 0.0);

	//Method
	settings.methodSettings.meshMotionType = "ALEMaterialInterfaces";
	//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";			
	//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		
	settings.methodSettings.spatialReconstruction = RTInstabilityTests::MethodSettings::SpatialDiscretisationType::Constant;

	//Run parameters
	int nSnapshots = 100;
	settings.MaxTime = 10e-6; //10 mks
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

#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "GridLoading\cgnsload.h"
#include "grid.h"
#include "gridgeneration.h"
#include "model.h"
#include "riemannsolvers.h"
#include "garbarukmodel.h"
#include "SAmodel.h"
#include "biffplatemodel.h"
#include "bicgstab.h"
#include "gridloading\triangleload.h"
#include "mpi.h"
#include "kernel.h"
#include "SODtest.h"

//Main program ))
int main(int argc, char *argv[]) {		
	//MPI_Init(&argc, &argv);
	//std::cout<<sizeof(MPI_LONG_DOUBLE)<<"\n";
	//int size;
	//MPI_Type_size( MPI_DOUBLE, &size );
	//std::cout<<size<<"\n";		
	//std::cout<<sizeof(double)<<"\n";
	//MPI_Finalize();

	//test SODtest class
	SodTest test(argc, argv);
	bool testResult = test.RunTest();
	return 0;

	//SimpleChannelTest(); return 0;
	//ImplicitSolverTestOneCell(); return 0;
	//BlasiusTest(); return 0;
	//BumpFlowTestExplicit(); return 0;
	//ImplicitSolverTest(); return 0;
	//LinearSolverTests(); return 0;
	//GodunovTests();
	//RunSAFlatPlate();
	//RunGAWCalculation();
	//RunPoiseuilleTest();
	//RunSODTest();
	//RunBiffFlatPlate();
    
	//RunSODTest();
	//return 0;

	//runSodTest(argc, argv);
	return 0;
	
	Kernel _kernel;	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 30, 90, -0.25, 0.25, -0.75, 0.75, 1.0, 1.0);
	_kernel.BindGrid(_grid);
	//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\Mixer\\Mixer.cgns");	
	//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns");	
	//_kernel.ReadInitialConditions("FlowSolution.E:1");
	_kernel.ReadConfiguration("");
	_kernel.InitCalculation();
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();	
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	//_kernel.SaveSolution();
	//_kernel.SaveGrid("grid.cgns");
	_kernel.Finalize();	
	//_kernel.LoadGridTopologyAndInfo("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\SODShockTube\\SodShockTube.cgns");
	//std::cout<<argc<<"\n"<<argv[1]<<"\n";
	
	//_kernel.LoadGridTopology("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns");
	//_kernel.PartitionGrid();
	//_kernel.GenerateGridGeometry();
	

	return 0;	
};
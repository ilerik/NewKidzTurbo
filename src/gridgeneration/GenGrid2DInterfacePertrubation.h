#ifndef NewKidzTurbo_gridgeneration_GenGrid2DInterfacePertrubation
#define NewKidzTurbo_gridgeneration_GenGrid2DInterfacePertrubation

#include "grid.h"
#include "cgnslib.h"
#include "GenGrid2DUniform.h"

#include <cassert>
#include <functional>

Grid GenGrid2DInterfacePertrubation(std::shared_ptr<ParallelManager> MPIManager, int N, int M, 
									double xMin, double xMax, 
									double yMin, double yMax , 
									double q_x, double q_y, 
									bool periodicX, bool periodicY,
									Vector interfaceCenter, Vector interfaceNormal,
									std::function<double(Vector r)> disturbanceDistance) 
{
	Grid grid = GenGrid2DUniform(MPIManager, N, M, xMin, xMax, yMin, yMax, 1.0, 1.0, periodicX, periodicY);		

	//Determine nodes that a closest to interface
	std::map<double, int> interfaceNodesIndexes;	
	std::map<double, double> interfaceNodesY;
	double interfaceY = interfaceCenter.y;
	for (Node& node : grid.localNodes) {
		double x = node.P.x;
		double y = node.P.y;
		if (interfaceNodesY.find(x) == std::end(interfaceNodesY)) {
			interfaceNodesY[x] = y;			
			continue;
		};
		Vector leader = Vector(x, interfaceNodesY[x], 0);
		double d = disturbanceDistance(node.P);
		double dLeader = disturbanceDistance(leader);		
		if ( std::abs(y - interfaceY) <  std::abs(interfaceNodesY[x] - interfaceY)) {
		//if ( std::abs(d) <  std::abs(dLeader)) {
			interfaceNodesY[x] = y;
			interfaceNodesIndexes[x] = node.GlobalIndex;
		};
	};		

	//Compute displacements
	std::vector<int> nodes;
	std::vector<Vector> displacements;
	nodes.clear();
	displacements.clear();
	for (std::pair<double, int> pair : interfaceNodesIndexes) {
		double x = pair.first;
		double dY = disturbanceDistance(Vector(x, interfaceNodesY[x], 0));								
		int nodeIndex = pair.second;

		nodes.push_back(nodeIndex);		
		displacements.push_back(Vector(0, -dY, 0));
	};
	
	//Add unmovable nodes
	for (Node& node : grid.localNodes) {
		double x = node.P.x;
		double y = node.P.y;

		//if ((x == -_widthLeft) || (x == _widthRight) || (std::abs(y) == _widthY / 2.0)) {				
		if ((y == yMin) || (y == yMax)) {				
			displacements.push_back(Vector(0, 0, 0));
			nodes.push_back(node.GlobalIndex);
			continue;
		};		
	};

	grid.PartitionGrid(MPIManager);
	grid.GenerateLocalCells(MPIManager->rank(), grid.cellsPartitioning);
	grid.GenerateLocalFaces(MPIManager->rank());
	grid.UpdateGeometricProperties();
	//grid.GenerateLocalGrid();

	//Move mesh
	MeshMovement moveHelper;
	moveHelper.meshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDWnoRotation;
	moveHelper.MoveNodes(grid, nodes, displacements);

	grid.PartitionGrid(MPIManager);
	grid.GenerateLocalCells(MPIManager->rank(), grid.cellsPartitioning);
	grid.GenerateLocalFaces(MPIManager->rank());
	grid.UpdateGeometricProperties();

	return grid;
};

#endif
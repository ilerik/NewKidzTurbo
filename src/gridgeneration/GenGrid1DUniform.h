#ifndef NewKidzTurbo_gridgeneration_GenGrid1DUniform
#define NewKidzTurbo_gridgeneration_GenGrid1DUniform

#include "grid.h"
#include "cgnslib.h"

#include <cassert>
#include <functional>
#include <memory>

std::unique_ptr<Grid> GenGrid1DUniform(std::shared_ptr<ParallelManager> MPIManager, int N, double xMin, double xMax, bool IsPeriodic)
{
	std::unique_ptr<Grid> grid = std::unique_ptr<Grid>(new Grid(MPIManager));

	double size_x = xMax - xMin;	
	int rank = grid->gridStructure.rank = MPIManager->rank();
	int nProc = grid->gridStructure.np = MPIManager->np();

	//Generate grid topology
	const double q_x = 1.0;
	N++;	
	std::vector<double> x_p(N);	
	grid->periodicNodesIdentityList.clear();
	double h_x = 1.0 * size_x * (1.0 - q_x) / (1.0 - pow(q_x, N-1));	
	for (int i = 0; i<N; i++) {
		if (q_x == 1.0) {
			x_p[i] = (1.0*size_x*i) / (N-1);
		} else {
			x_p[i] = h_x * (1.0 - pow(q_x, i)) / (1.0-q_x);
		};
	};
	x_p[0] = 0;
	x_p[N-1] = size_x;
	
	//Create nodes
	grid->localNodes.clear();	
	for (int i = 0; i<N; i++) {
		Node new_node;
		new_node.GlobalIndex = i;
		new_node.P.x = x_p[i] + xMin;
		new_node.P.y = 0;
		new_node.P.z = 0;
		grid->localNodes.push_back(new_node); 
	};
	
	//Create cells
	grid->nCells = (N-1);		
	for (int i = 0; i<(N-1); i++) {
		Cell c;
		c.GlobalIndex = i;
		c.Nodes.resize(2);
		c.Nodes[0] = i;
		c.Nodes[1] = i + 1;		
		c.CellCenter = Vector(0,0,0);
		for (int k = 0; k<c.Nodes.size(); k++) {
			c.CellCenter =  c.CellCenter + grid->localNodes[c.Nodes[k]].P;
		};
		c.CellCenter = (1.0/c.Nodes.size()) * c.CellCenter;		
		c.CellVolume = (grid->localNodes[c.Nodes[1]].P-grid->localNodes[c.Nodes[0]].P).mod();	
		c.CellHSize = 0;
		c.CGNSType = BAR_2;
		c.IsDummy = false;
		grid->Cells.push_back(c);
	};

	grid->nProperCells = grid->nCells;		

	//Dummy cells and connectivity info		
	for (int i = 0; i<(N-1); i++) {
		int cIndex = i;
		//Left neighbour
		int neighbour = -1;
		if ((i == 0) && (!IsPeriodic)) {				
			//Create dummy cell
			Cell dummyCell;
			dummyCell.GlobalIndex = grid->nCells++;
			dummyCell.NeigbourCells.push_back(cIndex);
			dummyCell.IsDummy = true;
			dummyCell.CGNSType = NODE;				
			dummyCell.Nodes.clear();
			dummyCell.Nodes.push_back(i);			
			dummyCell.BCMarker = 1;				
			grid->Cells.push_back(dummyCell);
			neighbour = dummyCell.GlobalIndex;
		};
		if ((i == 0) && (IsPeriodic)) {
			neighbour = N-2;
		};
		if (i != 0) neighbour = i-1;
		grid->Cells[cIndex].NeigbourCells.push_back(neighbour);

		//Right neighbour
		if ((i == N-2) && (!IsPeriodic)) {				
			//Create dummy cell
			Cell dummyCell;
			dummyCell.GlobalIndex = grid->nCells++;
			dummyCell.NeigbourCells.push_back(cIndex);
			dummyCell.IsDummy = true;
			dummyCell.CGNSType = NODE;				
			dummyCell.Nodes.clear();
			dummyCell.Nodes.push_back(i+1);			
			dummyCell.BCMarker = 2;				
			grid->Cells.push_back(dummyCell);
			neighbour = dummyCell.GlobalIndex;
		};
		if ((i == N-2) && (IsPeriodic)) neighbour = 0;			
		if (i != N-2) neighbour = i+1;
		grid->Cells[cIndex].NeigbourCells.push_back(neighbour);		
	};
	
	grid->nDummyCells = grid->Cells.size() - grid->nProperCells;

	//Add periodic boundary nodes identity information
	if (IsPeriodic) {		
		grid->periodicNodesIdentityList[ 0 ].insert( N-1 );
		grid->periodicNodesIdentityList[ N-1 ].insert( 0 );		
	};	

	//Fill in connectivity info
	grid->vdist.clear();
	grid->xadj.clear();
	grid->adjncy.clear();

	//TO DO generalize
	grid->gridStructure.cellsPart.resize(grid->nProperCells);
	grid->cellsPartitioning.resize(grid->nProperCells);
	int nProcessors = nProc;
	int vProc = grid->nProperCells / nProcessors;
	int vLeft = grid->nProperCells % nProcessors;
	grid->vdist.resize(nProcessors + 1);
	grid->vdist[0] = 0;
	for (int i = 0; i<nProcessors; i++) {
		grid->vdist[i+1] = grid->vdist[i] + vProc;
		if (vLeft > 0) {
			grid->vdist[i+1]++;
			vLeft--;
		};

		//Set initial partitioning
		for (int j = grid->vdist[i]; j < grid->vdist[i + 1]; j++) {
			grid->gridStructure.cellsPart[j] = i;
			grid->cellsPartitioning[j] = i;
		};
	};
	grid->vdist[nProcessors] = grid->nProperCells;

	int currentInd = 0;
	grid->xadj.push_back(currentInd);
	for (int i = grid->vdist[rank]; i < grid->vdist[rank + 1]; i++)
	{			
		Cell& cell = grid->Cells[i];				
		for (int j = 0; j<cell.NeigbourCells.size(); j++) {
			int nIndex = cell.NeigbourCells[j];
			if (nIndex >= grid->nProperCells) continue; //Skip dummy
			grid->adjncy.push_back(nIndex);
			currentInd++;
		};
		grid->xadj.push_back(currentInd); 
	};

	//Eliminate non local cells		
	grid->nCellsLocal = -grid->vdist[rank] + grid->vdist[rank + 1];		

	//Fill in grid info
	grid->gridInfo.CellDimensions = 1;	
	
	//Add boundary names
	if (!IsPeriodic) {
		grid->addPatch("left", 1);
		grid->addPatch("right", 2);
	};	

	grid->PartitionGrid(MPIManager);
	grid->GenerateLocalCells(MPIManager->rank(), grid->cellsPartitioning);
	grid->GenerateLocalFaces(MPIManager->rank());
	grid->UpdateGeometricProperties();
	MPIManager->Barrier();
	
	return grid;
};

#endif
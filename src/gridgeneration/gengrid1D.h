#ifndef TURBO_GRIDGENERATION_GENGRID1D
#define TURBO_GRIDGENERATION_GENGRID1D

#include "grid.h"
#include "cgnslib.h"
#include "parallelHelper.h"
#include <assert.h>

Grid GenGrid1D(ParallelHelper* pHelper, int N, double lBegin, double lEnd, Vector direction, bool IsPeriodic = false)
{
	Grid g;
	int rank = pHelper->getRank();
	int nProc = pHelper->getProcessorNumber();

	//Normalize direction
	direction = direction / direction.mod();

	//Grid parameters	
	g.gridInfo.GridDimensions = 1;
	g.gridInfo.CellDimensions = 1;
	g.gridInfo.nCoords = 3;

	//Generate grid topology	
	N++;
	std::vector<double> x_p(N);
	std::vector<double> y_p(N);
	std::vector<double> z_p(N);
	g.periodicNodesIdentityList.clear();
	double L = lEnd - lBegin;	
	for (int i = 0; i<N; i++) {
		double curL = lBegin + L * i / (N-1);
		Vector r = curL * direction;
		x_p[i] = r.x;
		y_p[i] = r.y;
		z_p[i] = r.z;
	};
		
	//Fill in nodes
	g.localNodes.clear();
	for (int i = 0; i<N; i++) {		
		Node new_node;
		new_node.GlobalIndex = i;
		new_node.P.x = x_p[i];
		new_node.P.y = y_p[i];
		new_node.P.z = z_p[i];		
		//Add node
		g.localNodes.push_back(new_node);  		
	};		

	//Create cells
	g.nCells = (N-1);		
	for (int i = 0; i<(N-1); i++) {
		Cell c;
		c.GlobalIndex = i;
		c.Nodes.resize(2);
		c.Nodes[0] = i;
		c.Nodes[1] = i + 1;			
		c.CGNSType = BAR_2;
		c.IsDummy = false;
		g.Cells.push_back(c);	
	};
	g.nProperCells = g.nCells;		

	//Dummy cells and connectivity info		
	for (int i = 0; i<(N-1); i++) {
		int cIndex = i;
		//Left neighbour
		int neighbour = -1;
		if ((i == 0) && (!IsPeriodic)) {				
			//Create dummy cell
			Cell dummyCell;
			dummyCell.GlobalIndex = g.nCells++;
			dummyCell.NeigbourCells.push_back(cIndex);
			dummyCell.IsDummy = true;
			dummyCell.CGNSType = NODE;				
			dummyCell.Nodes.clear();
			dummyCell.Nodes.push_back(i);			
			dummyCell.BCMarker = 1;				
			g.Cells.push_back(dummyCell);
			neighbour = dummyCell.GlobalIndex;
		};
		if ((i == 0) && (IsPeriodic)) {
			neighbour = N-2;
		};
		if (i != 0) neighbour = i-1;
		g.Cells[cIndex].NeigbourCells.push_back(neighbour);

		//Right neighbour
		if ((i == N-2) && (!IsPeriodic)) {				
			//Create dummy cell
			Cell dummyCell;
			dummyCell.GlobalIndex = g.nCells++;
			dummyCell.NeigbourCells.push_back(cIndex);
			dummyCell.IsDummy = true;
			dummyCell.CGNSType = NODE;				
			dummyCell.Nodes.clear();
			dummyCell.Nodes.push_back(i+1);			
			dummyCell.BCMarker = 2;				
			g.Cells.push_back(dummyCell);
			neighbour = dummyCell.GlobalIndex;
		};
		if ((i == N-2) && (IsPeriodic)) neighbour = 0;			
		if (i != N-2) neighbour = i+1;
		g.Cells[cIndex].NeigbourCells.push_back(neighbour);		
	};
	g.nDummyCells = g.Cells.size() - g.nProperCells;

	//Add periodic boundary nodes identity information
	if (IsPeriodic) {		
		g.periodicNodesIdentityList[ 0 ].insert( N-1 );
		g.periodicNodesIdentityList[ N-1 ].insert( 0 );		
	};	

	//Fill in connectivity info
	g.vdist.clear();
	g.xadj.clear();
	g.adjncy.clear();

	//TO DO generalize
	int nProcessors = pHelper->getProcessorNumber();
	int vProc = g.nProperCells / nProcessors;
	int vLeft = g.nProperCells % nProcessors;
	g.vdist.resize(nProcessors + 1);
	g.vdist[0] = 0;
	for (int i = 0; i<nProcessors; i++) {
		g.vdist[i+1] = g.vdist[i] + vProc;
		if (vLeft > 0) {
			g.vdist[i+1]++;
			vLeft--;
		};
	};
	g.vdist[nProcessors] = g.nProperCells;

	/*assert(nProc == 1);
	g.vdist.push_back(0);		
	g.vdist.push_back(g.nProperCells);*/

	int currentInd = 0;
	g.xadj.push_back(currentInd);
	for (int i = g.vdist[rank]; i < g.vdist[rank + 1]; i++)
	{			
		Cell& cell = g.Cells[i];				
		for (int j = 0; j<cell.NeigbourCells.size(); j++) {
			int nIndex = cell.NeigbourCells[j];
			if (nIndex >= g.nProperCells) continue; //Skip dummy
			g.adjncy.push_back(nIndex);
			currentInd++;
		};
		g.xadj.push_back(currentInd);
	};

	//Eliminate non local cells		
	g.nCellsLocal = -g.vdist[rank] + g.vdist[rank + 1];				
	
	//Add boundary names
	if (!IsPeriodic) {
		g.addPatch("left", 1);
		g.addPatch("right", 2);
	};	

	g.ConstructAndCheckPatches();
	pHelper->Barrier();
	
	return g;
};

#endif
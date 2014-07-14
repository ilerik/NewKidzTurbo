#ifndef TURBO_GRIDGENERATION_GENGRID2D
#define TURBO_GRIDGENERATION_GENGRID2D

#include "grid.h"
#include "cgnslib.h"
#include "parallelHelper.h"
#include <assert.h>

Grid GenGrid2D(ParallelHelper* pHelper, int N, int M, double size_x, double size_y, double q_x, double q_y, bool periodicX = false, bool periodicY = false)
{
	int rank = pHelper->getRank();
	int nProc = pHelper->getProcessorNumber();
	int cartI;
	int cartJ;

	//Generate grid topology
	Grid g;
	N++;
	M++;
	std::vector<double> x_p(N);
	std::vector<double> y_p(M);
	double h_x = 1.0 * size_x * (1.0 -q_x) / (1.0 - pow(q_x, N-1));
	double h_y = 1.0 * size_y * (1.0 -q_y) / (1.0 - pow(q_y, M-1));
	for (int i = 0; i<N; i++) {
		if (q_x == 1.0) {
			x_p[i] = (1.0*size_x*i) / (N-1);
		} else {
			x_p[i] = h_x * (1.0 - pow(q_x, i)) / (1.0-q_x);
		};
	};
	x_p[0] = 0;
	x_p[N-1] = size_x;
	for (int i = 0; i<M; i++) {
		if (q_y == 1.0) {
			y_p[i] = (1.0 *size_y*i) / (M-1);
		} else {
			y_p[i] = 1.0 * h_y * (1.0 - pow(q_y, i)) / (1.0-q_y);
		};
	};
	y_p[0] = 0;
	y_p[M-1] = size_y;
	g.localNodes.clear();
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<M; j++) {
			Node new_node;
			new_node.GlobalIndex = i + j * N;
			new_node.P.x = x_p[i];//(size_x*i) / (N-1);
			new_node.P.y = y_p[j];//(size_y*j) / (M-1);
			new_node.P.z = 0;
			g.localNodes.push_back(new_node); 
		};
	};	

	//Create cells
	g.nCells = (N-1)*(M-1);
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			Cell c;
			c.GlobalIndex = i + j*(N-1);
			c.Nodes.resize(4);
			c.Nodes[0] = i + j*N;
			c.Nodes[1] = i + (j+1)*N;
			c.Nodes[2] = (i+1) + (j+1)*N;
			c.Nodes[3] = (i+1) + j*N;			
			c.CellCenter = Vector(0,0,0);
			for (int k = 0; k<c.Nodes.size(); k++) {
				c.CellCenter =  c.CellCenter + g.localNodes[c.Nodes[k]].P;
			};
			c.CellCenter = (1.0/c.Nodes.size()) * c.CellCenter;		
			c.CellVolume = ((g.localNodes[c.Nodes[1]].P-g.localNodes[c.Nodes[0]].P).mod()) * ((g.localNodes[c.Nodes[2]].P-g.localNodes[c.Nodes[1]].P).mod());			
			c.CellHSize = 0;
			c.CGNSType = QUAD_4;
			c.IsDummy = false;
			g.Cells.push_back(c);
		};
	};
	g.nProperCells = g.nCells;	

	//Dummy cells and connectivity info
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			int cIndex = i + j*(N-1);
			//Left neighbour
			int neighbour = -1;
			if ((i == 0) && (!periodicX)) {				
				//Create dummy cell
				Cell dummyCell;
				dummyCell.GlobalIndex = g.nCells++;
				dummyCell.NeigbourCells.push_back(cIndex);
				dummyCell.IsDummy = true;
				dummyCell.BCMarker = 1;				
				g.Cells.push_back(dummyCell);
				neighbour = dummyCell.GlobalIndex;
			};
			if ((i == 0) && (periodicX)) neighbour = N-2 + j*(N-1);			
			if (i != 0) neighbour = i-1 + j*(N-1);
			g.Cells[cIndex].NeigbourCells.push_back(neighbour);

			//Right neighbour
			if ((i == N-2) && (!periodicX)) {				
				//Create dummy cell
				Cell dummyCell;
				dummyCell.GlobalIndex = g.nCells++;
				dummyCell.NeigbourCells.push_back(cIndex);
				dummyCell.IsDummy = true;
				dummyCell.BCMarker = 2;				
				g.Cells.push_back(dummyCell);
				neighbour = dummyCell.GlobalIndex;
			};
			if ((i == N-2) && (periodicX)) neighbour = 0 + j*(N-1);			
			if (i != N-2) neighbour = i+1 + j*(N-1);
			g.Cells[cIndex].NeigbourCells.push_back(neighbour);

			//Bottom neighbour
			if ((j == 0) && (!periodicY)) {				
				//Create dummy cell
				Cell dummyCell;
				dummyCell.GlobalIndex = g.nCells++;
				dummyCell.NeigbourCells.push_back(cIndex);
				dummyCell.IsDummy = true;
				dummyCell.BCMarker = 3;				
				g.Cells.push_back(dummyCell);
				neighbour = dummyCell.GlobalIndex;
			};
			if ((j == 0) && (periodicY)) neighbour = i + (M-2)*(N-1);			
			if (j != 0) neighbour = i + (j-1)*(N-1);
			g.Cells[cIndex].NeigbourCells.push_back(neighbour);

			//Top neighbour
			if ((j == M-2) && (!periodicY)) {				
				//Create dummy cell
				Cell dummyCell;
				dummyCell.GlobalIndex = g.nCells++;
				dummyCell.NeigbourCells.push_back(cIndex);
				dummyCell.IsDummy = true;
				dummyCell.BCMarker = 4;				
				g.Cells.push_back(dummyCell);
				neighbour = dummyCell.GlobalIndex;
			};
			if ((j == M-2) && (periodicY)) neighbour = i + 0*(N-1);			
			if (j != M-2) neighbour = i+1 + (j+1)*(N-1);
			g.Cells[cIndex].NeigbourCells.push_back(neighbour);
		}
	}
	

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
	//g.Cells.erase(g.Cells.begin() + g.vdist[rank + 1], g.Cells.end());
	//g.Cells.erase(g.Cells.begin(), g.Cells.begin() + g.vdist[rank]);	

	//Fill in grid info
	g.gridInfo.CellDimensions = 2;	
	
	//Add boundary names
	if (!periodicX) {
		g.addPatch("left", 1);
		g.addPatch("right", 2);
	};
	if (!periodicY) {
		g.addPatch("bottom", 3);
		g.addPatch("top", 4);
	};
	//g.ConstructAndCheckPatches();
	pHelper->Barrier();
	
	return g;
};

#endif
#ifndef TURBO_GRIDGENERATION_GENGRID3D
#define TURBO_GRIDGENERATION_GENGRID3D

#include "grid.h"
#include "cgnslib.h"
#include "parallelHelper.h"
#include <assert.h>

inline int GetCartIndex(int i, int j, int k, int N, int M, int K) {	
	i = i % N;
	j = j % M;
	k = k % K;
	if (i < 0) i += N;
	if (j < 0) j += M;
	if (k < 0) k += K;
	return i + j * N + k * N * M;
};

Grid GenGrid3D(ParallelHelper* pHelper, int N, int M, int K, 
			   double x_min, double x_max, //X axis range 
			   double y_min, double y_max, //Y axis range
			   double z_min, double z_max, //Z axis range
			   double q_x, double q_y, double q_z,
			   bool periodicX = false, bool periodicY = false, bool periodicZ = false
			   )
{
	double size_x = x_max - x_min;
	double size_y = y_max - y_min;
	double size_z = z_max - z_min;
	int rank = pHelper->getRank();
	int nProc = pHelper->getProcessorNumber();
	int cartI;
	int cartJ;

	//Generate grid topology
	Grid g;
	N++;
	M++;
	K++;
	std::vector<double> x_p(N);
	std::vector<double> y_p(M);
	std::vector<double> z_p(K);
	g.periodicNodesIdentityList.clear();
	double h_x = 1.0 * size_x * (1.0 -q_x) / (1.0 - pow(q_x, N-1));
	double h_y = 1.0 * size_y * (1.0 -q_y) / (1.0 - pow(q_y, M-1));
	double h_z = 1.0 * size_z * (1.0 -q_z) / (1.0 - pow(q_z, K-1));
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
	for (int i = 0; i<K; i++) {
		if (q_z == 1.0) {
			z_p[i] = (1.0 *size_z*i) / (M-1);
		} else {
			z_p[i] = 1.0 * h_z * (1.0 - pow(q_z, i)) / (1.0-q_z);
		};
	};
	z_p[0] = 0;
	z_p[M-1] = size_y;
	g.localNodes.clear();	
	for (int k = 0; k<K; k++) {
		for (int j = 0; j<M; j++) {
			for (int i = 0; i<N; i++) {
				Node new_node;
				new_node.GlobalIndex = i + j * N + k * N * M;
				new_node.P.x = x_p[i] + x_min;//(size_x*i) / (N-1);
				new_node.P.y = y_p[j] + y_min;//(size_y*j) / (M-1);
				new_node.P.z = z_p[k] + z_min;
				g.localNodes.push_back(new_node); 
			};
		};	
	};

	//Create cells
	g.nCells = (N-1)*(M-1)*(K-1);	
	for (int k = 0; k<(K-1); k++) {
		for (int j = 0; j<(M-1); j++) {
			for (int i = 0; i<(N-1); i++) {
				Cell c;
				c.GlobalIndex = GetCartIndex(i,j,k, N-1, M-1, K-1);
				c.Nodes.resize(8);
				c.Nodes[0] = i + j*N + k*N*M;
				c.Nodes[1] = i + (j+1)*N + k*N*M;
				c.Nodes[2] = (i+1) + (j+1)*N + k*N*M;
				c.Nodes[3] = (i+1) + j*N + k*N*M;
				c.Nodes[4] = i + j*N + (k+1)*N*M;
				c.Nodes[5] = i + (j+1)*N + (k+1)*N*M;
				c.Nodes[6] = (i+1) + (j+1)*N + (k+1)*N*M;
				c.Nodes[7] = (i+1) + j*N + (k+1)*N*M;
				/*c.CellCenter = Vector(0,0,0);
				for (int k = 0; k<c.Nodes.size(); k++) {
					c.CellCenter =  c.CellCenter + g.localNodes[c.Nodes[k]].P;
				};*/
				//c.CellCenter = (1.0/c.Nodes.size()) * c.CellCenter;		
				//c.CellVolume = ((g.localNodes[c.Nodes[1]].P-g.localNodes[c.Nodes[0]].P).mod()) * ((g.localNodes[c.Nodes[2]].P-g.localNodes[c.Nodes[1]].P).mod());			
				//c.CellHSize = 0;
				c.CGNSType = HEXA_8;
				c.IsDummy = false;
				g.Cells.push_back(c);
			};
		};
	};
	g.nProperCells = g.nCells;		

	//Dummy cells and connectivity info	
	for (int k = 0; k<(K-1); k++) {
		for (int j = 0; j<(M-1); j++) {
			for (int i = 0; i<(N-1); i++) {
				int cIndex = i + j*(N-1)+k*(N-1)*(M-1);
				//Left neighbour
				int neighbour = -1;
				if ((i == 0) && (!periodicX)) {				
					//Create dummy cell
					Cell dummyCell;
					dummyCell.GlobalIndex = g.nCells++;
					dummyCell.NeigbourCells.push_back(cIndex);
					dummyCell.IsDummy = true;
					dummyCell.CGNSType = QUAD_4;				
					dummyCell.Nodes.clear();
					dummyCell.Nodes.push_back(GetCartIndex(i,j,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j+1,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j+1,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j,k+1, N, M, K));
					dummyCell.BCMarker = 1;				
					g.Cells.push_back(dummyCell);
					neighbour = dummyCell.GlobalIndex;
				} else {
					neighbour = GetCartIndex(i-1,j,k, N-1, M-1, K-1);
				};
				g.Cells[cIndex].NeigbourCells.push_back(neighbour);

				//Right neighbour
				if ((i == N-2) && (!periodicX)) {				
					//Create dummy cell
					Cell dummyCell;
					dummyCell.GlobalIndex = g.nCells++;
					dummyCell.NeigbourCells.push_back(cIndex);
					dummyCell.IsDummy = true;
					dummyCell.CGNSType = QUAD_4;				
					dummyCell.Nodes.clear();
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j+1,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j+1,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j,k+1, N, M, K));
					dummyCell.BCMarker = 2;				
					g.Cells.push_back(dummyCell);
					neighbour = dummyCell.GlobalIndex;
				} else {
					neighbour = GetCartIndex(i+1,j,k, N-1, M-1, K-1);
				};
				g.Cells[cIndex].NeigbourCells.push_back(neighbour);

				//Top neighbour
				if ((j == 0) && (!periodicY)) {				
					//Create dummy cell
					Cell dummyCell;
					dummyCell.GlobalIndex = g.nCells++;
					dummyCell.NeigbourCells.push_back(cIndex);
					dummyCell.IsDummy = true;					
					dummyCell.CGNSType = QUAD_4;				
					dummyCell.Nodes.clear();
					dummyCell.Nodes.push_back(GetCartIndex(i,j,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j,k+1, N, M, K));
					dummyCell.BCMarker = 3;				
					g.Cells.push_back(dummyCell);
					neighbour = dummyCell.GlobalIndex;
				} else {
					neighbour = GetCartIndex(i,j-1,k, N-1, M-1, K-1);
				};
				g.Cells[cIndex].NeigbourCells.push_back(neighbour);

				//Bottom neighbour
				if ((j == M-2) && (!periodicY)) {				
					//Create dummy cell
					Cell dummyCell;
					dummyCell.GlobalIndex = g.nCells++;
					dummyCell.NeigbourCells.push_back(cIndex);
					dummyCell.IsDummy = true;
					dummyCell.CGNSType = QUAD_4;				
					dummyCell.Nodes.clear();
					dummyCell.Nodes.push_back(GetCartIndex(i,j+1,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j+1,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j+1,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j+1,k+1, N, M, K));		
					dummyCell.BCMarker = 4;				
					g.Cells.push_back(dummyCell);
					neighbour = dummyCell.GlobalIndex;
				}else {
					neighbour = GetCartIndex(i,j+1,k, N-1, M-1, K-1);
				};
				g.Cells[cIndex].NeigbourCells.push_back(neighbour);

				//Rear neighbour
				if ((k == 0) && (!periodicZ)) {				
					//Create dummy cell
					Cell dummyCell;
					dummyCell.GlobalIndex = g.nCells++;
					dummyCell.NeigbourCells.push_back(cIndex);
					dummyCell.IsDummy = true;					
					dummyCell.CGNSType = QUAD_4;				
					dummyCell.Nodes.clear();
					dummyCell.Nodes.push_back(GetCartIndex(i,j,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j+1,k, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j+1,k, N, M, K));
					dummyCell.BCMarker = 5;				
					g.Cells.push_back(dummyCell);
					neighbour = dummyCell.GlobalIndex;
				} else {
					neighbour = GetCartIndex(i,j,k-1, N-1, M-1, K-1);
				};
				g.Cells[cIndex].NeigbourCells.push_back(neighbour);

				//Front neighbour
				if ((k == K-2) && (!periodicZ)) {				
					//Create dummy cell
					Cell dummyCell;
					dummyCell.GlobalIndex = g.nCells++;
					dummyCell.NeigbourCells.push_back(cIndex);
					dummyCell.IsDummy = true;
					dummyCell.CGNSType = QUAD_4;				
					dummyCell.Nodes.clear();
					dummyCell.Nodes.push_back(GetCartIndex(i,j,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i,j+1,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j+1,k+1, N, M, K));
					dummyCell.Nodes.push_back(GetCartIndex(i+1,j,k+1, N, M, K));		
					dummyCell.BCMarker = 6;				
					g.Cells.push_back(dummyCell);
					neighbour = dummyCell.GlobalIndex;
				}else {
					neighbour = GetCartIndex(i,j,k+1, N-1, M-1, K-1);
				};
				g.Cells[cIndex].NeigbourCells.push_back(neighbour);
			}
		}
	};
	g.nDummyCells = g.Cells.size() - g.nProperCells;

	//Add periodic boundary nodes identity information
	if (periodicX) {
		for (int k = 0; k<K; k++) {
			for (int j = 0; j<M; j++) {
				g.periodicNodesIdentityList[ GetCartIndex(0, j, k, N, M, K) ].insert( GetCartIndex(-1, j, k, N, M, K) );
				g.periodicNodesIdentityList[ GetCartIndex(-1, j, k, N, M, K) ].insert( GetCartIndex(0, j, k, N, M, K) );
			};
		};
	};
	if (periodicY) {
		for (int k = 0; k<K; k++) {
			for (int i = 0; i<N; i++) {
				g.periodicNodesIdentityList[ GetCartIndex(i, 0, k, N, M, K) ].insert( GetCartIndex(i, -1, k, N, M, K) );
				g.periodicNodesIdentityList[ GetCartIndex(i, -1, k, N, M, K) ].insert( GetCartIndex(i, 0, k, N, M, K) );
			};
		};
	};
	if (periodicZ) {
		for (int j = 0; j<M; j++) {
			for (int i = 0; i<N; i++) {
				g.periodicNodesIdentityList[ GetCartIndex(i, j, 0, N, M, K) ].insert( GetCartIndex(i, j, -1, N, M, K) );
				g.periodicNodesIdentityList[ GetCartIndex(i, j, -1, N, M, K) ].insert( GetCartIndex(i, j, 0, N, M, K) );
			};
		};
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
	//g.Cells.erase(g.Cells.begin() + g.vdist[rank + 1], g.Cells.end());
	//g.Cells.erase(g.Cells.begin(), g.Cells.begin() + g.vdist[rank]);	

	//Fill in grid info	
	g.gridInfo.CellDimensions = 3;	
	
	//Add boundary names
	if (!periodicX) {
		g.addPatch("left", 1);
		g.addPatch("right", 2);
	};
	if (!periodicY) {
		g.addPatch("bottom", 3);
		g.addPatch("top", 4);
	};
	if (!periodicZ) {
		g.addPatch("rear", 5);
		g.addPatch("front", 6);
	};

	//g.ConstructAndCheckPatches();
	pHelper->Barrier();
	
	return g;
};

#endif
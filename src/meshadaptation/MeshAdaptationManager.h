#ifndef TURBO_MeshAdaptation_MeshAdaptation
#define TURBO_MeshAdaptation_MeshAdaptation

#include "grid.h"
#include "basetypes.h"
#include "meshquality.h"
#include "geomfunctions.h"
#include <cassert>
#include <unordered_map>
#include <unordered_set>

//Impelementation of mesh movement algorithms
class MeshAdaptationManager {
	std::shared_ptr<Grid> _gridPtr;		//Reference to underlying grid
	std::unordered_map<int, std::set<int>> _balls;	//Balls of elements for each node
	std::unordered_map<int, double> _cellsQuality;	//Quality for all grid non dummy cells

	//Utility functions
	std::unordered_map<int, std::set<int>>& ComputeBallsOfElements() {
		for (Cell* cellPtr : _gridPtr->localCells) {
			for (int node : cellPtr->Nodes) _balls[node].insert(cellPtr->GlobalIndex);
		};
		return _balls;
	};

	//Compute given cell quality
	double ComputeCellQuality(Cell& cell) {
		double Q = MeshQuality::Anisotropy(_gridPtr, cell);
		return Q;
	};

	//Compute cells quality
	std::unordered_map<int, double>& ComputeCellsQuality() {
		for (Cell* cellPtr : _gridPtr->localCells) {
			_cellsQuality[cellPtr->GlobalIndex] = ComputeCellQuality(*cellPtr);
		};
		return _cellsQuality;
	};
public:
	//Constructor
	MeshAdaptationManager(std::shared_ptr<Grid>& gridPtr) : _gridPtr(gridPtr) {
	};	

	//Decompose cell in all possible ways
	std::vector<std::vector<Cell> > DecomposeCell(const Cell& cell) {
		cell.Nodes;
	};

	//Compose cell from several cells
	Cell ComposeCell(std::vector<const Cell&> cells) {
		Cell unitedCell;
		//Build cells ajacency graph

		//Unite cells one by one

		//Compute convex hull
	};

	//Find optimal position for a given node
	Vector FindOptimalPosition(int node) {
		std::set<int> BOE = _balls[node];

		//Compute minimum quality for neighbour cells given node position
		auto computeMinQuality = [&](Vector P) {
			
		};
		
	};

	//Apply vertex smoothing algorithm
	void ApplySmoothing(std::unordered_set<int> freeNodes) {
		//For each free node define it's ball of elements
		ComputeBallsOfElements();

		//Compute quality for all cells
		
		//For each node try to improve minimum quality
		for (int node : freeNodes) {
			Vector optimalP = FindOptimalPosition(node, nodePositions);
		};
	};
	
	//Apply topology changing adaptation
	void ApplyRecombination(std::unordered_map<int, std::vector<double>> cellValues) {
		//For QUAD_4 cells
		
	};
};

#endif
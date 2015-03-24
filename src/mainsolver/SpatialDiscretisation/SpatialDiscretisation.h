#ifndef NewKidzTurbo_MainSolver_SpatialDiscretisation_SpatialDiscretisation
#define NewKidzTurbo_MainSolver_SpatialDiscretisation_SpatialDiscretisation

#include <vector>

#include "grid.h"

//Satial discretisation type
enum class SpatialDiscretisationType {
	PiecewiseConstant,
	PiecewiseLinear,
	WENO
};

//Cell spatial reconstruction info
class CellSpatialDiscretisation {
	//std::shared_ptr<Cell> _cellPtr;
	Cell* _cellPtr;
	std::vector<int> _stencilIndexes;
	bool _initialized;
	int _localIndex;	

	//Linear reconstruction for now TO DO downshift to child class
	std::vector<double> _meanValues;
	/*Vector gradRo;
	Vector gradRoU;
	Vector gradRoV;
	Vector gradRoW;
	Vector gradRoE;*/
public:
	//Constructor
	CellSpatialDiscretisation(Cell* cell) : _initialized(false), _cellPtr(cell) {		
	};

	//Public properties
	inline const bool initialized() { return _initialized; };
	inline const std::vector<int> stencil() { return _stencilIndexes; };
	inline const int localIndex() { return _localIndex; };
	inline const int globalIndex() { return _cellPtr->GlobalIndex; };
	
	//Interface functions		
	virtual std::vector<int> CalculateStencil() {		
		_stencilIndexes.clear();
		_stencilIndexes.push_back(_cellPtr->GlobalIndex); //Cell itself
		return _stencilIndexes;
	};

	//Reconstruct cell solution given stencil cells and values
	virtual void ReconstructSolution(std::vector<std::vector<double> >& values) {
		_meanValues = values[0];
	};

	//Interpolate face values
	virtual std::vector<double> GetSolutionAtFace(Face& face) {
		return _meanValues;
	};
};

//Spatial discretisation implementation
class SpatialDiscretisation {
	std::shared_ptr<Grid> _gridPtr;
	SpatialDiscretisationType _discretisationType;
public:
	//Public constructor
	SpatialDiscretisation(std::shared_ptr<Grid> gridPtr) : _gridPtr(gridPtr) {

	};

	//Public properties
	inline const SpatialDiscretisationType type() { return _discretisationType; };
	
	//Interface functions

	//Obtain stencil for a given cell
	std::vector<int> GetStencilIndexes(int localIndex) {
		CellSpatialDiscretisation cellInfo(_gridPtr->localCells[localIndex]);
		return cellInfo.CalculateStencil();
	};

	//Reconstruct solution in cell
	CellSpatialDiscretisation ReconstructCellSolution(Cell* cell, std::vector<std::vector<double> >& values) {
		CellSpatialDiscretisation cellInfo(cell);
		cellInfo.CalculateStencil();		
		cellInfo.ReconstructSolution(values);
		return cellInfo;
	};
};

#endif
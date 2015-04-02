#ifndef NewKidzTurbo_MainSolver_SpatialDiscretisation_SpatialDiscretisation
#define NewKidzTurbo_MainSolver_SpatialDiscretisation_SpatialDiscretisation

#include <vector>

#include "grid.h"

//Satial discretisation type
enum class SpatialDiscretisationType {
	PiecewiseConstant,	
	WENO,
	ENO	
};

//Stencil structure
struct Stencil {
	//Stencil graph structure with root vertex pointing at center
	int root;
	std::set<int> cells;
	std::map<int, std::vector<int> > neighbours;
	inline const int& size() { return cells.size(); };
	inline const bool& has_cell(int ind) { return (cells.find(ind) != cells.end()); };
	
	//Construct empy stencil
	Stencil() {
		root = -1;
		cells.clear();
	};

	//Construct one cell stencil
	Stencil(int root_) {
		root = root_;
		neighbours[root] = std::vector<int>();
		cells.insert(root);
	};	

	//Construct of disconnected cells
	Stencil(int root_, std::set<int> cells_) {
		root = root_;
		cells = cells_;
		neighbours.clear();
	};

	//Add cell to stencil and make connections
	bool AddCell(int cell, std::vector<int> nCells) {
		//Return false if this cell is in stencil
		if (cells.find(cell) != cells.end()) return false; 

		//Add cell
		cells.insert(cell);
		neighbours[cell] = nCells;
		for (int nCell : nCells) neighbours[nCell].push_back(cell);
		return true;
	};	
};


struct ReconstructionSettings {
	SpatialDiscretisationType type; //Reconstruction procedure type
	int nDims;
	int nVars;
	int minStencilDepth;
	int stencilCells;	
};

class ReconstructionData {
public:
	//Stencils used
	Vector center;
	Stencil stencil;

	//Oscilation indicator values	
	std::vector<double> oscilation;
	std::vector<double> weight;	
	double distance;	

	//Reconstruction data	
	std::vector<double> meanValues;
	std::vector<Vector> gradients;

	//Get interpolated values at point
	std::vector<double> InterpolateValues(SpatialDiscretisationType type, Vector r) {
		std::vector<double> result = meanValues;
		if ((type == SpatialDiscretisationType::WENO) || (type == SpatialDiscretisationType::ENO)) {
			Vector dr = r - center;
			result += gradients * dr;
			return result;
		};

		//Default piecewise constant reconstruction
		return result;
	};

	//Get interpolated value at point
	double InterpolateValue(int nValue, SpatialDiscretisationType type, Vector r) {
		double result = meanValues[nValue];
		if ((type == SpatialDiscretisationType::WENO) || (type == SpatialDiscretisationType::ENO)) {
			Vector dr = r - center;
			result += gradients[nValue] * dr;
			return result;
		};

		//Default piecewise constant reconstruction
		return result;
	};
};


//Cell spatial reconstruction info
class CellSpatialDiscretisation {
	//Base data
	std::shared_ptr<Grid> _gridPtr;		//Grid object reference		
	Stencil _stencil;					//Current cell global stencil	
	ReconstructionSettings _settings;	//Settings

	std::vector<ReconstructionData> _reconstructions;	//Reconstruction data for each substencil
	std::vector<std::vector<double>> _weights;			//Weight of each substencil
	std::vector<double> _mino; //ENO stencilt oscilation for every variable
	std::vector<int> _ENOStencilIndex;	//ENO stencil index for every variable
	
	//Temp
	std::vector<int> _stencilIndexes;
	bool _initialized;
	int _localIndex;		
public:
	//Constructor
	CellSpatialDiscretisation(int cell, std::shared_ptr<Grid>& grid, SpatialDiscretisationType type) : 
		_initialized(false),
		_stencil(cell),
		_gridPtr(grid)
	{
		_settings.nDims = _gridPtr->gridInfo.CellDimensions;
		_settings.type = type;
		_settings.nVars = 5;
		if (type == SpatialDiscretisationType::PiecewiseConstant) {
			_settings.minStencilDepth = 0;
			_settings.stencilCells = 1;
		};
		if ((type == SpatialDiscretisationType::WENO) || (type == SpatialDiscretisationType::ENO))
		{
			_settings.minStencilDepth = 1;
			_settings.stencilCells = 1 + _settings.nDims;
		};		
	};

	//Public properties
	inline const bool initialized() { return _initialized; };
	inline const Stencil& stencil() { return _stencil; };
	inline const std::vector<int>& stencilIndexes() { return _stencilIndexes; };
	inline const int localIndex() { return _localIndex; };
	inline const int globalIndex() { return _stencil.root; };
	
	//Interface functions		
	const Stencil& CalculateStencil( std::function<bool(Cell&)> isGoodCell ) {		
		//Required number of neighbours
		std::deque<int> BFSQueue;
		std::set<int> visited;
		std::map<int, int> distance;

		//Prepare
		BFSQueue.push_back(_stencil.root);
		distance[_stencil.root] = 0;
		visited.clear();
		
		//Prepare output data structures		
		_stencilIndexes.clear();		

		//Breadth first search for neighbours
		while ( !BFSQueue.empty() ) {
			//Get current cell index from queue
			int currentCellIndex = BFSQueue.at(0);
			BFSQueue.pop_front();
			visited.insert(currentCellIndex);
			Cell& currentCell = _gridPtr->Cells[currentCellIndex];

			//If we excided predefined depth and have enough cells to get all needed info
			if ((distance[currentCellIndex] > _settings.minStencilDepth) && (_stencilIndexes.size() >= _settings.stencilCells)) break;

			//Add eligible cells to stencil
			_stencilIndexes.push_back(currentCellIndex);
			_stencil.AddCell(currentCellIndex, currentCell.NeigbourCells);
			
			for (int& nIndex : currentCell.NeigbourCells) {
				//Skip if already visited
				if (visited.find(nIndex) != visited.end())  continue;

				//Obtain neighbour cell
				Cell& nCell = _gridPtr->Cells[nIndex];

				//Skip dummy
				if (!isGoodCell(nCell)) continue;

				//Add neighbour
				BFSQueue.push_back(nIndex);
				distance[nIndex] = distance[currentCellIndex] + 1;
			};
		};

		//Remove all edges that have cells outside of stencil as vertex		
		std::set<int> cellsToRemove;
		for (auto& pair : _stencil.neighbours) {
			int cell = pair.first;
			if (!_stencil.has_cell(cell)) {
				cellsToRemove.insert(cell);		
			};			

			std::vector<int> neighboursLeft;
			for (int ncell : pair.second) if (_stencil.has_cell(ncell)) {
				neighboursLeft.push_back(ncell);
			};
			pair.second = neighboursLeft;
		};

		for (auto cell : cellsToRemove) {			
			_stencil.neighbours.erase(cell);						
		};

		return _stencil; 
	};

	//Obtain all substencils of given size
	std::vector<Stencil> ObtainSubstencils(int size) {	
		//Store cells taken
		std::set<int> stencilCells;
		std::set<std::set<int> > subStencilsCells;
		std::vector<Stencil> subStencils;

		//Obtain all possible substencils
		std::function<void()> StencilSearch = [&]() {
			//If we gathered enough
			if (size == stencilCells.size()) {						
				subStencilsCells.insert(stencilCells);
				return;
			};

			//Find all cells neighbours
			std::unordered_set<int> neighbourCells;
			for (int cell : stencilCells) {
				auto ncells = _stencil.neighbours[cell];
				for (int ncell : ncells) neighbourCells.insert(ncell);
			};

			//Remove cells that already in stencil
			for (int cell : stencilCells) {
				if (neighbourCells.find(cell) != std::end(neighbourCells)) neighbourCells.erase(cell);
			};

			//Make step to all possible next choises
			for (int ncell : neighbourCells) {
				stencilCells.insert(ncell);
				StencilSearch();
				stencilCells.erase(ncell);
			};

			return;
		};

		//Launch stencil search	
		stencilCells.insert(_stencil.root);
		StencilSearch();

		//Generate stencils
		for (std::set<int> cells : subStencilsCells) {
			subStencils.push_back(Stencil(_stencil.root, cells));
		};

		return subStencils; 
	};

	//Make reconstruction on given stencil
	ReconstructionData ReconstructStencilSolution(const Stencil& stencil, std::map<int, std::shared_ptr<std::vector<double>>>& values) {
		ReconstructionData data;
		int nv = _settings.nVars;
		int rootIndex = stencil.root;
		Cell& rootCell = _gridPtr->Cells[rootIndex];
		Vector center = rootCell.CellCenter;
		data.stencil = stencil;
		data.center = center;
		data.meanValues = *values[rootIndex];

		//Cell characteristic sizes
		double hx = 0;
		double hy = 0;
		double hz = 0;
		for (int nI = 0; nI < rootCell.Nodes.size(); nI++) {
			Vector dr = _gridPtr->localNodes[nI].P - center;
			hx += std::abs(dr.x);
			hy += std::abs(dr.y);
			hz += std::abs(dr.z);
		};

		if (_settings.type == SpatialDiscretisationType::PiecewiseConstant) {
			//Null gradients
			data.gradients = std::vector<Vector>(nv, Vector());
			data.oscilation.resize(nv, 1.0);
			data.weight.resize(nv, 1.0);
		};
		if ((_settings.type == SpatialDiscretisationType::WENO) || (_settings.type == SpatialDiscretisationType::ENO)) {
			bool isGradientFailed = false;
			bool status = true;
			data.meanValues = *values[rootIndex];
			data.gradients = std::vector<Vector>(nv, Vector());
			for (int i = 0; i < nv; i++) {
				//Solve LLS problem to find gradient
				std::vector<double> vs(0);
				std::vector<Vector> points(0);
				for (int cellIndex : stencil.cells) {
					if (cellIndex == rootIndex) continue; //Skip center
					Cell& cell = _gridPtr->Cells[cellIndex];
					points.push_back(cell.CellCenter);
					vs.push_back((*values[cellIndex])[i]);					
				};

				//Actual LLS call
				Vector grad = ComputeGradientByPoints(_settings.nDims, center, data.meanValues[i], points, vs, status);
				data.gradients[i] = grad;

				//If reconstruction failed finish
				if (!status) {
					isGradientFailed = true;
				};
			};

			//Compute oscilation indicator value based on obtained reconstruction		
			data.oscilation.resize(nv,0);
			data.weight.resize(nv, 0);

			if (isGradientFailed) {
				//If no reconstruction obtained
				for (int i = 0; i < nv; i++) {
					data.oscilation[i] = std::numeric_limits<double>::max();
					data.weight[i] = 0;
				};
			} else {
				//If we have obtained correct reconstruction				

				//Compute smoothness (oscilation)			
				for (int cellIndex : stencil.cells) {
					if (cellIndex == rootIndex) continue; //Skip center
					Cell& cell = _gridPtr->Cells[cellIndex]; 					
					for (int i = 0; i < nv; i++) {
						double dux = (data.gradients[i].x) * hx;
						double duy = (data.gradients[i].y) * hy;
						double duz = (data.gradients[i].z) * hz;
						data.oscilation[i] += (dux * dux) + (duy * duy) + (duz * duz);						
					};
				};						

				//Compute weight
				const double r = -4;
				const double epsilon = std::numeric_limits<double>::epsilon();
				for (int i = 0; i < nv; i++) {
					data.weight[i] = std::pow(epsilon + data.oscilation[i], r);
				};
			};

			
		};
		values.clear();

		return data;
	};

	//Reconstruct cell solution given stencil cells and values
	virtual void ReconstructSolution(std::vector< int >& materials, std::vector<std::vector<double> >& values) {		
		//Now obtain all availible stencils
		int nv = _settings.nVars;
		int size = _settings.stencilCells;
		std::vector<Stencil> subStencils = ObtainSubstencils(size);		

		//Make map from cell index to values
		std::map<int, std::shared_ptr<std::vector<double>> > valuesMap;
		for (int i = 0; i < _stencilIndexes.size(); i++) {
			int cell = _stencilIndexes[i];			
			valuesMap[cell] = std::shared_ptr<std::vector<double>>(new std::vector<double>(values[i]));
		};

		//Initialize equal weights
		int nstencils = subStencils.size();
		_weights.clear();//.resize(nstencils); 
		/*for (int i = 0; i < nstencils; i++) {
			_weights[i].resize(nv, 1.0 / subStencils.size());
		};*/

		//Now compute smoothness of every stencil and make reconstruction
		_ENOStencilIndex.resize(nv, -1);
		_mino.resize(nv, std::numeric_limits<double>::max());
		std::map<int, std::shared_ptr<std::vector<double>> > vs;
		for (int sI = 0; sI < subStencils.size(); sI++) {
			Stencil& s = subStencils[sI];

			//Gather values
			for (int cell : s.cells) {
				vs[cell] = valuesMap[cell];
			};	

			//Compute reconstruction
			ReconstructionData data = ReconstructStencilSolution(s, vs);
			_reconstructions.push_back(data);
			_weights.push_back(data.weight);

			//Find ENO stencil
			for (int i = 0; i < nv; i++) {
				if (data.oscilation[i] < _mino[i]) {
					_mino[i] = data.oscilation[i];
					_ENOStencilIndex[i] = sI;
				};
			};
		};

		//Normilize weights		
		for (int i = 0; i < nv; i++) {
			//Component wise
			double wSum = 0;
			for (int j = 0; j < nstencils; j++) {						
				wSum += _weights[j][i];
			};
			for (int j = 0; j < nstencils; j++) {						
				_weights[j][i] /= wSum;
			};
		};
	
	};

	//Interpolate face values
	virtual std::vector<double> GetSolutionAtFace(Face& face) {		
		if ((_settings.type == SpatialDiscretisationType::ENO)) {
			//ENO
			//Obtain interpolated values at face center
			std::vector<double> U(_settings.nVars);
			for (int i = 0; i < _settings.nVars; i++) {
				U[i] = _reconstructions[_ENOStencilIndex[i]].InterpolateValue(i, _settings.type, face.FaceCenter);
			};
			return U;
		};

		//WENO general reconstruction
		int nStencils = _weights.size();
		std::vector<double> result(_settings.nVars, 0.0);
		std::vector<double> wSum(_settings.nVars); //Sum of distance corrected weights at face
		for (int si = 0; si < nStencils; si++) {						
			//Obtain interpolated values at face center
			for (int i = 0; i < _settings.nVars; i++) {				
				double u =  _reconstructions[si].InterpolateValue(i, _settings.type, face.FaceCenter);
				double w = _weights[si][i];				
				result[i]	+= w * u;
				wSum[i]		+= w;
			};

			//Compute distance coefficient for current stencil (as inverse average distance to face center)
			double d = 0.0;						
			for (int cellIndex : _reconstructions[si].stencil.cells) {	
				Cell& cell = _gridPtr->Cells[cellIndex];
				Vector dr = cell.CellCenter - face.FaceCenter;
				d += dr.mod();
			};
			d = 1.0 / std::pow(d, 5);						
		};

		//Normilize		
		for (int i = 0; i < _settings.nVars; i++) {								
			result[i] /=  wSum[i];
		};		

		return result;
	};
};

#endif
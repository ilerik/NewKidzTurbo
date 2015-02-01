#ifndef TURBO_Utility_MeshQuality
#define TURBO_Utility_MeshQuality

//Various mesh elements quality measures

namespace MeshQuality {

//Normilized anisotropy
inline double Anisotropy(const Grid& grid, const Cell& cell) {
	//Meaningless for 1D
	if (grid.gridInfo.CellDimensions == 1) {
		return 1.0;
	};

	double sumFaceSquare = 0; //Sum of face squares
	double normRatio = 0; //Ratio of perimeter to volume for perfect polygon	
	double exponent = 0;
	
	switch (cell.CGNSType) {
	case (TRI_3) : {
		//Triangle
		normRatio = 12 * std::sqrt(3);
		exponent = 2; //2D case
		break;
		};
	case (QUAD_4) : {
		//Quadriliteral
		normRatio = 16;
		exponent = 2; //2D case
		break;
		};
	default : {
		//
		throw new Exception("Unknown type of element");
		};
	};

	//Compute sum of face normals	
	for (int faceInd : cell.Faces) {
		double faceS = grid.localFaces[faceInd].FaceSquare;
		sumFaceSquare += faceS;
	};

	//Compute normalized quality
	if (grid.gridInfo.CellDimensions == 2) {
		double q = std::pow(sumFaceSquare, 2.0) / cell.CellVolume;
		q /= normRatio;
		return q;
	};
	if (grid.gridInfo.CellDimensions == 3) {
		double q = std::pow(sumFaceSquare, 3.0 / 2.0) / cell.CellVolume;
		q /= normRatio;
		return q;
	};

	throw new Exception("Unsupported number of dimensions");
	return 0;
};


//Cell linear size
inline double LinearSize(const Grid& grid, const Cell& cell) {
	if (grid.gridInfo.CellDimensions>0) {
		double linearSize = std::pow(cell.CellVolume, 1.0/grid.gridInfo.CellDimensions);
		return linearSize;
	};

	throw new Exception("Unsupported number of dimensions");
	return 0;
};

//Cell linear size ratio
inline double LinearSizeRatio(const Grid& grid, const Cell& cell) {
	if (grid.gridInfo.CellDimensions<=0) throw new Exception("Unsupported number of dimensions");
	
	//Compute ratio of linear sizes of given cell with every neighbour
	double maxRatio = 1.0;
	std::vector<double> ratios(0);
	double linearSize = std::pow(cell.CellVolume, 1.0/grid.gridInfo.CellDimensions);
	
	for (int nGlobalIndex : cell.NeigbourCells) {
		int nLocalIndex = grid.cellsGlobalToLocal.find(nGlobalIndex)->second;
		double nVolume = grid.localCells[nLocalIndex]->CellVolume;
		double nLinearSize = std::pow(nVolume, 1.0/grid.gridInfo.CellDimensions);
		double ratio = nLinearSize / linearSize;
		if (ratio < 1.0) ratio = 1.0 / ratio; //If less then one take inverse
		if (maxRatio < ratio) maxRatio = ratio;
	};

	return maxRatio;
};


};

#endif
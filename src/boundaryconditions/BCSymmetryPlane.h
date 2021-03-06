#ifndef TURBO_BoundaryConditions_BCSymmetryPlane
#define TURBO_BoundaryConditions_BCSymmetryPlane

#include "BoundaryCondition.h"

namespace BoundaryConditions {

class BCSymmetryPlane : public BoundaryCondition {
public:	
	virtual std::vector<double> getDummyValues(int nmat, std::vector<double> values, const Cell& dummyCell) {
		//Obtain face
		if (dummyCell.Faces.size() != 1) throw new Exception("Dummy cell has more than one face");
		Face& face = _grid->localFaces[dummyCell.Faces[0]];
		//Obtain neighbour cell
		int nCellIndex = _grid->cellsGlobalToLocal[face.FaceCell_1];
		int nVariables = _gasModels[nmat]->nConservativeVariables;
		//Obtain interrior values
		std::vector<double> inV(&values[nCellIndex * nVariables], &values[nVariables * nCellIndex] + nVariables);

		//Compute dummy values
		std::vector<double> res = inV;	
		Vector roU(inV[1], inV[2], inV[3]);			
		Vector newRoU = roU - 2*(roU * face.FaceNormal)*face.FaceNormal/face.FaceNormal.mod();
		res[1] = newRoU.x;
		res[2] = newRoU.y;
		res[3] = newRoU.z;		
		return res;
	};		

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {			
		//Read configuration
	}; 
};

};

#endif
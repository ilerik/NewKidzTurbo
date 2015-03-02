#ifndef TURBO_BoundaryConditions_BCOutflowSupersonic
#define TURBO_BoundaryConditions_BCOutflowSupersonic

#include "BoundaryCondition.h"

namespace BoundaryConditions {

class BCOutflowSupersonic : public BoundaryCondition {
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
		
		return res;
	};

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {			
		//Read configuration
	}; 
};

};

#endif
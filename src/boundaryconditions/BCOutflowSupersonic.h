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
		double pressure = 1e11;
		/*
		double density = 1.0;
		res[0] = density;
		res[1] = 0;
		res[2] = 0;
		res[3] = 0;
		res[4] = pressure / (0.4);*/
		//res[4] = pressure / (0.4);
		/*res[0] = inV[0];
		res[1] = inV[1];
		res[2] = inV[2];
		res[3] = inV[3];
		res[4] = 0;*/
		Vector roU(inV[1], inV[2], inV[3]);					
		return res;
	};

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {			
		//Read configuration
	}; 
};

};

#endif
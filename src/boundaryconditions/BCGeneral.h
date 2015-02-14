#ifndef TURBO_BoundaryConditions_BCGeneral
#define TURBO_BoundaryConditions_BCGeneral

#include "BoundaryCondition.h"

namespace BoundaryConditions {

enum class BoundaryVariableType {
	Pressure,
	Density,
	VelocityX,
	VelocityY,
	VelocityZ,
	InternalEnergy
};

class CompositeBoundaryConditionInfo {
public:
	double CValue;
	double CGradient;
	double Value;

	//Set values
	void SetValues(double cValue, double cGradient, double value) {
		if (cValue > 0 || cGradient > 0) {
			CValue = cValue;
			CGradient = cGradient;
			Value = value;
			return;
		};
		throw new Exception("One of coefficients must be non-zero");
	};

	//Set dirichlet boundary condition
	void SetDirichletBoundary(double Value) {
		SetValues(1.0, 0.0, Value);
	};

	//Set neuman boundary condition
	void SetNeumanBoundary(double Gradient) {
		SetValues(0.0, 1.0, Gradient);
	};

	//Interpolate
	double GetDummyValue(const double inV, const Face& f, const Vector cellCenter) {
		double dn = -(cellCenter - f.FaceCenter) * f.FaceNormal;
		double a = (Value - CValue * inV) / (CValue + CGradient / dn);
		double b = inV + a;
		double value = a + b;
		return value;
	};
};

class BCGeneral : public BoundaryConditions::BoundaryCondition {
public:
	//Boundary conditions in general form
	std::map<BoundaryVariableType, CompositeBoundaryConditionInfo> boundaryConditions;

	//Get dummy cell values
	std::vector<double> getDummyValues(int nmat, std::vector<double> values, const Cell& dummyCell) {
		//Obtain face
		if (dummyCell.Faces.size() != 1) throw new Exception("Dummy cell has more than one face");
		Face& face = _grid->localFaces[dummyCell.Faces[0]];

		//Obtain neighbour cell
		int nCellIndex = _grid->cellsGlobalToLocal[face.FaceCell_1];
		int nVariables = _gasModels[nmat]->nConservativeVariables;
		Vector center = _grid->localCells[nCellIndex]->CellCenter;

		//Obtain interrior values
		std::vector<double> inV(&values[nCellIndex * nVariables], &values[nVariables * nCellIndex] + nVariables);

		//Compute dummy values
		double ro = inV[0];
		double u = inV[1]/inV[0];
		double v = inV[2]/inV[0];
		double w = inV[3]/inV[0];
		double E = inV[4]/inV[0]; 
		double e = ro*E - ro*(u*u + v*v + w*w)/2.0;

		double roDummy = boundaryConditions[BoundaryVariableType::Density].GetDummyValue(ro, face, center);
		double uDummy = boundaryConditions[BoundaryVariableType::VelocityX].GetDummyValue(u, face, center);
		double vDummy = boundaryConditions[BoundaryVariableType::VelocityY].GetDummyValue(v, face, center);
		double wDummy = boundaryConditions[BoundaryVariableType::VelocityZ].GetDummyValue(w, face, center);
		double eDummy = boundaryConditions[BoundaryVariableType::InternalEnergy].GetDummyValue(e, face, center);

		std::vector<double> res(nVariables);	
		res[0] = roDummy;
		res[1] = roDummy * uDummy;
		res[2] = roDummy * vDummy;
		res[3] = roDummy * wDummy;
		res[4] = eDummy + roDummy*(uDummy*uDummy + vDummy*vDummy + wDummy*wDummy)/2.0;
		return res;
	};

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {			
		//Read configuration
		//std::pair<double, bool> propertyValue;
		//propertyValue = bcConfig.GetPropertyValue("Density");
		//Density = propertyValue.first;		

		//TO DO figure out configuration
		boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::InternalEnergy] = CompositeBoundaryConditionInfo();

		boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::InternalEnergy].SetDirichletBoundary(0);
	}; 

};

};

#endif
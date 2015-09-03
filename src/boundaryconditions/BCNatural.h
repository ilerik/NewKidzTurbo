#ifndef TURBO_BoundaryConditions_BCNatural
#define TURBO_BoundaryConditions_BCNatural

#include "BCGeneral.h"

namespace BoundaryConditions {

class BCNatural : public BoundaryConditions::BCGeneral {
public:

	virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) override {			
		//Read configuration
		//std::pair<double, bool> propertyValue;
		//propertyValue = bcConfig.GetPropertyValue("Density");
		//Density = propertyValue.first;		

		//TO DO figure out configuration
		boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

		boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);
	}; 

};

};

#endif
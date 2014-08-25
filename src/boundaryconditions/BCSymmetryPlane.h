#ifndef TURBO_BoundaryConditions_BCSymmetryPlane
#define TURBO_BoundaryConditions_BCSymmetryPlane

#include "BoundaryCondition.h"

namespace BoundaryConditions {

class BCSymmetryPlane : public BoundaryCondition {
public:	
	std::vector<double> getDummyValues(std::vector<double> inV, const Face& face) {
		std::vector<double> res = inV;				
		return res;
	};		

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {			
		//Read configuration
	}; 
};

};

#endif
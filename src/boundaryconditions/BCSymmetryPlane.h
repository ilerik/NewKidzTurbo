#ifndef TURBO_BoundaryConditions_BCSymmetryPlane
#define TURBO_BoundaryConditions_BCSymmetryPlane

#include "BoundaryCondition.h"

namespace BoundaryConditions {

class BCSymmetryPlane : public BoundaryCondition {
public:	
	virtual std::vector<double> getDummyValues(std::vector<double> inV, const Face& face) {
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
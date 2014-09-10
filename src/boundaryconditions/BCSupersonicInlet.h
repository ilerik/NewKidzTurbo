#ifndef TURBO_BoundaryConditions_BCSupersonicInlet
#define TURBO_BoundaryConditions_BCSupersonicInlet

#include "BoundaryCondition.h"

class BCSupersonicInlet : public BoundaryConditions::BoundaryCondition {
public:
	Vector Velocity;

	std::vector<double> getDummyValues(std::vector<double> inV, , const Cell& dummyCell) {
		std::vector<double> res = inV;
		/*Vector V_in(inV.rou, inV.rov, inV.row);
		V_in = (1.0/inV.ro)*V_in;

		Vector V_dummy = 2.0*V - V_in;
		res.rou = inV.ro*V_dummy.x;
		res.rov = inV.ro*V_dummy.y;
		res.row = inV.ro*V_dummy.z;
		res.roE = inV.roE - 0.5*inV.ro*V_in.mod()*V_in.mod() + 0.5*res.ro*V_dummy.mod()*V_dummy.mod();*/
		return res;
	};		

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {			
		//Read configuration

	}; 

};

#endif
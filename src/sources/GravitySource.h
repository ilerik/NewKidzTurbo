#ifndef NewKidzTurbo_Sources_GravitySource
#define NewKidzTurbo_Sources_GravitySource

#include "ISource.h"

//Source term base class
class GravitySource : public ISource {	
	Vector _g; // free fall acceleration
public:		
	//Constructor	
	GravitySource(Vector g) {
		_g = g;
	};

	virtual std::vector<double> GetResidual(std::vector<double>& values) override {		
		double ro = values[0];
		double u = values[1] / values[0];
		double v = values[2] / values[0];
		double w = values[3] / values[0];
		Vector velocity(u,v,w);
		
		//Compute residual
		std::vector<double> R(5, 0);
		R[0] = 0;
		R[1] = ro * _g.x;
		R[2] = ro * _g.y;
		R[3] = ro * _g.z;
		R[4] = ro * velocity * _g;

		return R;
	};
};

#endif

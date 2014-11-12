#ifndef TURBO_GasModels_PerfectGasModel
#define TURBO_GasModels_PerfectGasModel

#include "GasModel.h"

//Gas model base class
class PerfectGasModel : public GasModel {
public:	
	//Medium properties
	double Gamma;
	double Cv;
	double Cp;	
	double R;		

	//Obtain medium pressure
	virtual double GetPressure(GasModel::ConservativeVariables U) {
		double ro = U.ro;
		double u = U.rou / ro;
		double v = U.rov / ro;
		double w = U.row / ro;		
		double P = (U.roE - ro*(u*u + v*v + w*w) / 2.0) * (Gamma - 1.0);
		return P;
	};

	//Obtain medium pressure, soundspeed, Gruneisen coefficient and adiabatic exponent values
	void GetPressureAndSoundSpeed(GasModel::ConservativeVariables U, double& pressure, double& soundspeed, double& gruneisen) {
		double ro = U.ro;
		double u = U.rou / ro;
		double v = U.rov / ro;
		double w = U.row / ro;	
		double e = (U.roE / ro) - (u*u + v*v + w*w) / 2.0;
		pressure = (ro * e) * (Gamma - 1.0);
		soundspeed = sqrt((Gamma - 1.0) * (e + pressure/ro));
		gruneisen = Gamma - 1.0;		
	};		


	//Read configuration
	void loadConfiguration(Configuration _configuration) {
		nConservativeVariables = 5;
		R = _configuration.IdealGasConstant;
		Gamma = _configuration.SpecificHeatRatio;
		Cp = _configuration.SpecificHeatPressure;
		Cv = _configuration.SpecificHeatVolume;		
	};
};

#endif
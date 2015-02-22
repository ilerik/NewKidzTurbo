#ifndef TURBO_GasModels_PerfectGasModel
#define TURBO_GasModels_PerfectGasModel

#include "GasModel.h"

//Gas model base class
class PerfectGasModel : public GasModel {
public:	
	//Constructor
	PerfectGasModel(Logger& logger) : GasModel(&logger) {
		//Model information
		GasModelName = "PerfectGasModel";
		GasModelType = ModelType_t::CaloricallyPerfect;		
	};

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

	//Obtain medium temperature
	virtual double GetTemperature(GasModel::ConservativeVariables U) {
		double ro = U.ro;
		double u = U.rou / ro;
		double v = U.rov / ro;
		double w = U.row / ro;	
		double e = (U.roE / ro) - (u*u + v*v + w*w) / 2.0;
		double T = e / Cv;
		return T;
	};

	//Obtain information about phase
	virtual GasModel::MediumPhase GetPhase(GasModel::ConservativeVariables U) {
		return GasModel::MediumPhase::BelowMeltingPoint;
	};

	//Given density and pressure find internal energy
	virtual double FindInternalEnergy(double ro, double p) {
		double e = p / (ro * (Gamma - 1.0));
		return e;
	};


	//Read configuration
	void loadConfiguration(GasModelConfiguration conf) {
		nConservativeVariables = 5;
		std::pair<double, bool> res;
		//Load properties
		res = conf.GetPropertyValue("IdealGasConstant");
		R = res.first;
		res = conf.GetPropertyValue("SpecificHeatRatio");
		Gamma = res.first;

		//Specific heats
		res = conf.GetPropertyValue("SpecificHeatPressure");
		bool isCpSet = res.second;
		Cp = res.first;
		res = conf.GetPropertyValue("SpecificHeatVolume");		
		bool isCvSet = res.second;
		Cv = res.first;
		
		//Check if everything is fine
		if ((!isCvSet) && (!isCpSet)) {
			//Error
			//TO DO
			throw 1;
		};		
	};
};

#endif
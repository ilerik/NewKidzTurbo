#ifndef TURBO_GasModels_BarotropicGasModel
#define TURBO_GasModels_BarotropicGasModel

#include "GasModel.h"

//Gas model base class
class BarotropicGasModel : public GasModel {
public:	
	//Constructor
	BarotropicGasModel(Logger& logger) : GasModel(&logger) {
		//Model information
		GasModelName = "BarotropicGasModel";
		GasModelType = ModelType_t::ModelTypeUserDefined;		
	};

	//Medium properties
	double E;
	double ReferencePressure;
	double ReferenceDensity;		

	//Obtain medium pressure
	virtual double GetPressure(GasModel::ConservativeVariables U) {
		double ro = U.ro;			
		double P = ReferencePressure + E * (ro - ReferenceDensity) / ReferenceDensity;
		return P;
	};

	//Obtain medium pressure, soundspeed, Gruneisen coefficient and adiabatic exponent values
	void GetPressureAndSoundSpeed(GasModel::ConservativeVariables U, double& pressure, double& soundspeed, double& gruneisen) {
		double ro = U.ro;		
		pressure = ReferencePressure + E * (ro - ReferenceDensity) / ReferenceDensity;
		soundspeed = std::sqrt(E / ReferenceDensity);
		gruneisen = 0.0;		
	};

	//Obtain medium temperature
	virtual double GetTemperature(GasModel::ConservativeVariables U) {
		double T = 0;
		return T;
	};

	//Obtain information about phase
	virtual GasModel::MediumPhase GetPhase(GasModel::ConservativeVariables U) {
		return GasModel::MediumPhase::BelowMeltingPoint;
	};

	//Given density and pressure find internal energy
	virtual double FindInternalEnergy(double ro, double p) {
		double e = 0;
		return e;
	};

	//Read configuration
	void loadConfiguration(GasModelConfiguration conf) {
		nConservativeVariables = 5;
		std::pair<double, bool> res;
		//Load properties
		res = conf.GetPropertyValue("YoungModulus");
		E = res.first;			
		res = conf.GetPropertyValue("ReferencePressure");
		ReferencePressure = res.first;	
		res = conf.GetPropertyValue("ReferenceDensity");
		ReferenceDensity = res.first;			
	};
};

#endif
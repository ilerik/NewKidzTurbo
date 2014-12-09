#ifndef TURBO_GasModels_GasModel
#define TURBO_GasModels_GasModel

#include "cgnslib.h"
#include "configuration.h"

//Gas model base class
class GasModel {
public:
	//Thermodynamic gas model
	std::string GasModelName;
	ModelType_t GasModelType;

	//Number of variables
	int nConservativeVariables;

	//Conservative variables class
	struct ConservativeVariables {
	public:
		double ro; //Density
		double rou; //Specific momentum X
		double rov; //Specific momentum Y
		double row; //Specific momentum Z
		double roE; //Specific total energy

		inline operator std::vector<double>() const {
			std::vector<double> u(5);
			u[0] = ro;
			u[1] = rou;
			u[2] = rov;
			u[3] = row;
			u[4] = roE;
			return u;
		};

		ConservativeVariables() {
		};

		ConservativeVariables(double* value) {
			ro = value[0];
			rou = value[1];
			rov = value[2];
			row = value[3];
			roE = value[4];
		};

		ConservativeVariables(std::vector<double> value) {
			ro = value[0];
			rou = value[1];
			rov = value[2];
			row = value[3];
			roE = value[4];
		};
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
	virtual void GetPressureAndSoundSpeed(GasModel::ConservativeVariables U, double& pressure, double& soundspeed, double& gruneisen) = 0;

	//Read configuration
	virtual void loadConfiguration(GasModelConfiguration _configuration) = 0;

	//Storing / saving service functions
	std::vector<std::string> GetStoredFieldsNames() {
		std::vector<std::string> names(nConservativeVariables);
		names[0] = "Density";
		names[1] = "MomentumX";
		names[2] = "MomentumY";
		names[3] = "MomentumZ";
		names[4] = "EnergyStagnationDensity";
		return names;
	};

	//Get values for storing from conservative variables
	std::vector<double> GetStoredValuesFromConservative(const ConservativeVariables& u) {
		return (std::vector<double>)u;
	};

	//Initialize conservative variables from stored vector
	ConservativeVariables GetConservativeFromStoredValues(const std::vector<double> values) {
		ConservativeVariables u;
		u.ro = values[0];
		u.rou = values[1];
		u.rov = values[2];
		u.row = values[3];
		u.roE = values[4];
		return u;
	};
};

#endif

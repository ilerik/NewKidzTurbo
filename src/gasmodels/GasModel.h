#ifndef TURBO_GasModels_GasModel
#define TURBO_GasModels_GasModel

#include "cgnslib.h"

//Gas model base class
class GasModel {
	//Thermodynamic gas model
	ModelType_t GasModelType;

	//Number of variables
	int nConservativeVariables;

	//Conservative variables class
	struct ConservativeVariables {
		double ro; //Density
		double rou; //Specific momentum X
		double rov; //Specific momentum Y
		double row; //Specific momentum Z
		double roE; //Specific total energy

		//Initialize conservative from vector	
	};

	std::vector<std::string> GetSolutionFieldsNames() {
		std::vector<std::string> names(nConservativeVariables);
		names[0] = "Density";
		names[1] = "MomentumX";
		names[1] = "MomentumY";
		names[1] = "MomentumZ";
		names[4] = "EnergyStagnationDensity";
		return names;
	};

	//Medium properties
	double Cv;
	double Cp;	
	double R;
	
	//Read configuration
	void loadConfiguration() {
	};
};

#endif

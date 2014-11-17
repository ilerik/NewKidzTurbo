#ifndef TURBO_MAINSOLVER_CONFIGURATION_Configuration
#define TURBO_MAINSOLVER_CONFIGURATION_Configuration

#include "BoundaryConditionConfiguration.h"
#include "GasModelConfiguration.h"
#include "grid.h"

//Class that represents configuration file methods and structure
class Configuration {
public:
	//Files locations
	std::string InputCGNSFile;
	std::string OutputCGNSFile;
	std::string InitialConditionsSolution;	

	//Set of equations settings	
	GoverningEquationsType_t GoverningEquations;

	//Thermodynamic Gas Model Structure
	ModelType_t GasModel;
	double IdealGasConstant;
	double SpecificHeatRatio;
	double SpecificHeatVolume;
	double SpecificHeatPressure;

	//Boundary conditions		
	std::map<std::string, BoundaryConditionConfiguration> BoundaryConditions;
};

#endif
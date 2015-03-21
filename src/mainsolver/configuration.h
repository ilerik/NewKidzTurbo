#ifndef TURBO_MAINSOLVER_CONFIGURATION_Configuration
#define TURBO_MAINSOLVER_CONFIGURATION_Configuration

#include "BoundaryConditionConfiguration.h"
#include "GasModelConfiguration.h"
#include "ALEConfiguration.h"
#include "RiemannSolverConfiguration.h"
#include "grid.h"

//Class that represents configuration file methods and structure
class Configuration {
	std::map<std::string, double> _properties;

	//Helpers
	//Checks if property was set
	bool IsPropertySet(std::string name) {
		return _properties.find(name) != _properties.end();
	};

	//Returns property value and set status
	std::pair<double, bool> GetPropertyValue(std::string name) {
		if (!IsPropertySet(name)) return std::pair<double, bool>(0, false);
		return std::pair<double, bool>(_properties[name], true);
	};

	//Sets property value
	void SetPropertyValue(std::string name, double value) {
		_properties[name] = value;
	};
public:
	//Files locations
	std::string WorkingDirectory;
	std::string InputCGNSFile;
	std::string OutputCGNSFile;
	std::string InitialConditionsSolution;	

	//Simulation parameters
	SimulationType_t SimulationType; //Simulation type
	double CFL; //CFL number
	int RungeKuttaOrder; //Runge-Kutta time stepping scheme order

	//Gravity
	Vector g;
	
	//ALE configuration
	ALEConfiguration ALEConfiguration;

	//Riemann solver configuration
	RiemannSolverConfiguration RiemannSolverConfiguration;

	//Run settings
	int MaxIteration; //Maximum number of iterations
	double MaxTime;   //Maximum computational time [s]
	int SaveSolutionSnapshotIterations; //Save solution snapshot after number of iterations
	double SaveSolutionSnapshotTime;	//Save solution snapshot after time interval

	//Set of equations settings	
	GoverningEquationsType_t GoverningEquations;

	//Thermodynamic Gas Model Structure
	std::vector<std::string> GasModelNames;
	std::map<std::string, int> GasModelNameToIndex;
	std::map<std::string, GasModelConfiguration> GasModelsConfiguration;
	
	//Add new gas model to configuration
	void AddGasModel(std::string name) {
		for (std::string& existingName : GasModelNames) {
			if (existingName == name) {
				return; //Already exist
			};
		};
		
		//Add new gas model
		GasModelNames.push_back(name);
		GasModelNameToIndex[name] = GasModelNames.size() - 1;		
	};

	//Boundary conditions		
	std::map<std::string, BoundaryConditionConfiguration> BoundaryConditions;
};

#endif
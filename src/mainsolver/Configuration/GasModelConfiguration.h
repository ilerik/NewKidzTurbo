#ifndef TURBO_MAINSOLVER_CONFIGURATION_GasModelConfiguration
#define TURBO_MAINSOLVER_CONFIGURATION_GasModelConfiguration

#include "cgnslib.h"

class GasModelConfiguration
{
	std::map<std::string, double> _properties;
public:
	//Gas model name (type)
	std::string GasModelName;

	//CGNS style types
	ModelType_t GasModelType;

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

	//Read configuration element from stream
	void ReadFromStream(std::istringstream s) {
		double value;
		std::string token, name;
		s>>token;
		if (token != "{") return;
		for (;;) {
			s>>token;
			if (token == "}") break;		
			double value;
			s>>name>>token>>value;			
			SetPropertyValue(name, value);
		};
	};
};

#endif
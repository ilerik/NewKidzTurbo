#ifndef TURBO_DATATYPES
#define TURBO_DATATYPES

#include <map>
#include <string>

class Dictionary {
private:
	std::map<std::string, double> _values;
public:
	bool isSet(std::string name) {
		return _values.find(name) != _values.end();
	};

	double getValue(std::string name) {
		if (!isSet(name)) throw Exception("Property \""+name+"\" must be set");
		return _values[name];
	};

	double setValue(std::string name, double v) {
		_values[name] = v;
	};
	
};

//Conservative variables
struct ConservativeVariables {
	int GlobalIndex;
	double ro;	//Density
	double rou;	//VelocityX * Density
	double rov;	//VelocityY * Density
	double row;	//VelocityZ * Density
	double roE; //Total energy * Density

	ConservativeVariables() {		 
	};

	ConservativeVariables(int ind) {
		GlobalIndex = ind;
	};

	//TO DO - VEctor ro
	void FillValuesFromPrimitiveVariables(Vector ro, double p, double T) {

	};

	double& operator[](int i) {
		/*double res;
		switch (i) {
		case 0: {
					res = ro;
					break;
				};
		case 1: {
					res = rou;
					break;
				};
		case 2: {
				res = rov;
				break;
			};
		case 3: {
				res = row;
				break;
			};
		case 4: {
				res = roE;
				break;
			};
		};
		return res;*/	
		switch (i) {
		case 0: {
					return ro;
					break;
				};
		case 1: {
					return rou;
					break;
				};
		case 2: {
					return rov;
					break;
				};
		case 3: {
					return row;
					break;
				};
		case 4: {
					return roE;
					break;
				};
		};
	};
};

#endif

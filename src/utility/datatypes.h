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
		GlobalIndex = 0;
		ro = 0;
		rou = 0;
		rov = 0;
		row = 0;
		roE = 0;
	};

	ConservativeVariables(int ind) {
		GlobalIndex = ind;
		ro = 0;
		rou = 0;
		rov = 0;
		row = 0;
		roE = 0;
	};

	double& operator[](int i) {
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

	ConservativeVariables& operator+=(const std::vector<double>& b) {
		ro += b[0];
		rou += b[1];
		rov += b[2];
		row += b[3];
		roE += b[4];
		return *this;
	};

	ConservativeVariables& operator+=(const ConservativeVariables& b) {
		ro += b.ro;
		rou += b.rou;
		rov += b.rov;
		row += b.row;
		roE += b.roE;
		return *this;
	};

	ConservativeVariables& operator*=(const double& c) {
		ro *= c;
		rou *= c;
		rov *= c;
		row *= c;
		roE *= c;
		return *this;
	};
};

//Overload some operators
ConservativeVariables operator+(const ConservativeVariables& a, const ConservativeVariables& b) {
	ConservativeVariables res;
	res.ro = a.ro + b.ro;
	res.rou = a.rou + b.rou;
	res.rov = a.rov + b.rov;
	res.row = a.row + b.row;
	res.roE = a.roE + b.roE;
	return res;
};

ConservativeVariables operator*(const ConservativeVariables& a, const double& c) {
	ConservativeVariables res(a);					
	res *= c;
	return res;
};

#endif

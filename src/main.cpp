#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "kernel.h"
#include "gridgeneration.h"
#include <cmath>

template< typename T >
std::string int_to_hex( T i )
{
  std::stringstream stream;
  stream << "0x" 
         << std::setfill ('0') << std::setw(sizeof(T)*2) 
         << std::hex << i;
  return stream.str();
}

//SOD Test new version
class SodTestInitialConditions : public InitialConditions::InitialConditions
{
public:
	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;
			//http://en.wikipedia.org/wiki/Sod_shock_tube
			//Grid sizes
			double Lx = 0.02; // SI 2 cm				

			//Other velocities
			double v = 0;
			double w = 0;

			//Left state
			double roL = 1.0;
			double PL = 1.0;
			double uL = 0;
			
			//Right state
			double roR = 0.125;
			double PR = 0.1;
			double uR = 0;

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
				
			//Values
			double u = 0;
			double ro = 0;
			double roE = 0;
			if (x <= 0.5) {
				ro = roL;
				u = uL;
				roE = PL/(_gasModel->Gamma - 1.0);
			} else {
				ro = roR;
				u = uR;
				roE = PR/(_gasModel->Gamma - 1.0);
			};
			
			//Convert to conservative variables
			initValues.resize(_gasModel->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0);

			return initValues;		
	};
};

void runSodTest(int argc, char *argv[]) {
	Kernel _kernel;
	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 200, 1, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, true, true);
	_kernel.BindGrid(_grid);
	_kernel.ReadConfiguration("");		
	_kernel.InitCalculation();

	//Initial conditions
	SodTestInitialConditions ic;
	_kernel.GenerateInitialConditions(ic);	

	//Run test
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();

	//Output result
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	_kernel.Finalize();	
};

//Periodic boundary test
class PeriodicTestInitialConditions : public InitialConditions::InitialConditions
{
public:
	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;
			//Other velocities
			double ro = 1.0;
			double u = 1;
			double v = 1;
			double w = 1;
			double P = 1.0;

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
				
			//Values
			double roE = P/(_gasModel->Gamma - 1.0);;
			
			//Convert to conservative variables
			initValues.resize(_gasModel->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0);

			return initValues;		
	};
};

void runPeriodicTest2D(int argc, char *argv[]) {
	Kernel _kernel;
	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 20, 20, 
		0.0, 1.0,
		0.0, 1.0,
		1.0, 1.0,
		true, true);
	_kernel.BindGrid(_grid);
	_kernel.ReadConfiguration("");		
	_kernel.InitCalculation();

	//Initial conditions
	PeriodicTestInitialConditions ic;
	_kernel.GenerateInitialConditions(ic);	

	//Run test
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();

	//Output result
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	_kernel.Finalize();	
};

void runPeriodicTest3D(int argc, char *argv[]) {
	Kernel _kernel;
	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid3D(_kernel.getParallelHelper(), 10, 10, 10, 
		0.0, 1.0, 
		0.0, 1.0, 
		0.0, 1.0,
		1.0, 1.0, 1.0, 
		true, true, true);
	_kernel.BindGrid(_grid);
	_kernel.ReadConfiguration("");		
	_kernel.InitCalculation();

	//Initial conditions
	PeriodicTestInitialConditions ic;
	_kernel.GenerateInitialConditions(ic);	

	//Run test
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();

	//Output result
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	_kernel.Finalize();	
};


//Shear layer computation
class ShearLayerInitialConditions : public InitialConditions::InitialConditions
{
private:
	double _topBound;
	double _bottomBound;	
	double _ro;
	double _topV;
	double _bottomV;
	double _P;	
	double _PI;
public:	

	ShearLayerInitialConditions() {
		_PI = 3.14159265359;
		_topBound = _PI + 0.5;
		_bottomBound = _PI - 0.5;
		_ro = 1.0;
		_P = 1000;
		_topV = 5;
		_bottomV = -5;
	};

	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;
			//Other velocities
			double ro = _ro;
			double u = 0;
			double v = 0;
			double w = 0;
			double P = _P;

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Below bottom bound
			if (z < _bottomBound) {
				u = 0;
				v = _bottomV;
				w = 0;
			};

			//Medium layer
			if ((_bottomBound <= z) && (z <= _topBound)) {
				double alpha = (z - _bottomBound) / (_topBound - _bottomBound);
				u = 0.2*sin(y)*cos(x);
				v = (1.0 - alpha) * _bottomV + alpha * _topV;
				w = 0.2*sin(y)*cos(x);
			};

			//Above top bound
			if (z > _topBound) {
				u = 0;
				v = _topV;
				w = 0;
			};
				
			//Values
			double roE = P/(_gasModel->Gamma - 1.0) + ro*(u*u + v*v + w*w) / 2.0;
			
			//Convert to conservative variables
			initValues.resize(_gasModel->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0);

			return initValues;		
	};
};

void runShearLayer(int argc, char *argv[]) {
	Kernel _kernel;
	const double PI = 3.14159265359;
	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid3D(_kernel.getParallelHelper(), 40, 40, 40,
		0.0, 2*PI,
		0.0, 2*PI,
		0.0, 2*PI,
		1.0, 1.0, 1.0,	
		true, true, false);
	_kernel.BindGrid(_grid);
	_kernel.ReadConfiguration("");		
	_kernel.InitCalculation();

	//Initial conditions
	ShearLayerInitialConditions ic;
	_kernel.GenerateInitialConditions(ic);	

	//Run test
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();

	//Output result
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	_kernel.Finalize();	
};

//Main program ))
 int main(int argc, char *argv[]) {		
	//runPeriodicTest2D(argc, argv);
	//runSodTest(argc, argv);
	runShearLayer(argc, argv);
	return 0;	
};
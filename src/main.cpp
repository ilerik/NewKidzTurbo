#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "kernel.h"
#include "gridgeneration.h"
#include <cmath>
#include "LomonosovFortovGasModel.h"

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
	//Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 200, 1, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, false, true);
	Grid _grid = GenGrid1D(_kernel.getParallelHelper(), 1000, 0, 1.0, false);
	_kernel.BindGrid(_grid);
	//_kernel.LoadGrid("C:\\Users\\Ilya\\Dropbox\\Science\\ValidationCFD\\Mixer\\Mixer.cgns");	
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
	double _Pressure;
	double _topBound;
	double _bottomBound;	
	double _ro;
	double _topV;
	double _bottomV;	
	double _PI;
public:	

	ShearLayerInitialConditions() {
		_PI = 3.14159265359;
		_topBound = _PI + 0.5;
		_bottomBound = _PI - 0.5;
		_ro = 1.0;
		_Pressure = 1000.0;
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
			double P = _Pressure;

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;

			//Below bottom bound
			if (y < _bottomBound) {
				u = _bottomV;
				v = 0;
				w = 0;
			};

			//Medium layer
			if ((_bottomBound <= y) && (y <= _topBound)) {
				double alpha = (y - _bottomBound) / (_topBound - _bottomBound);
				//u = 0.2*sin(y)*cos(x);
				//u = 0.2*sin(8*y)*cos(8*x);
				//u = 0.5*sin(2*ae_pi*y)*cos(ae_pi*x);
				//v = (1.0 - alpha) * _bottomV + alpha * _topV;
				//w = 0.2*sin(y)*cos(x);
				//w = 0.2*sin(8*y)*cos(8*x);				
				//w = 0.5*sin(2*ae_pi*y)*cos(ae_pi*x);
				u = (1.0 - alpha) * _bottomV + alpha * _topV;
				double lambda = 1;
				double A = 1.0;
				double B = 1.0;
				double C = 1.0;
				w += A*cos(lambda * x) + B*sin(lambda*y);
				v += C*cos(lambda * z) + A*sin(lambda*x);
				u += B*cos(lambda * y) + C*sin(lambda*z);
			};

			//Above top bound
			if (y > _topBound) {
				u = _topV;
				v = 0;
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

//Shear layer computation
class ImpactShockInitialConditions : public InitialConditions::InitialConditions
{
private:
	//State of target	
	double _u0;
	double _ro0;
	//Speed of projectile	
	double _V;	

	
public:	

	ImpactShockInitialConditions(double V, double ro0) {
		_V = V;
		_ro0 = ro0;		
	};

	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;				
		//Other velocities
		double v = 0;
		double w = 0;

		//Internal energy
		double e = 10;

		//Left state		
		double roL = _ro0;		
		double uL = 0.5*_V;
			
		//Right state		
		double roR = _ro0;		
		double uR = -0.5*_V;

		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
				
		//Values
		double u = 0;
		double ro = 0;
		double roE = 0;				
		if (x <= 0) {
			ro = roL;
			u = uL;			
		} else {
			ro = roR;
			u = uR;						
		};
			
		//Convert to conservative variables
		roE = ro*(e + (u*u + v*v + w*w) / 2.0);
		initValues.resize(_gasModel->nConservativeVariables);
		initValues[0] = ro;
		initValues[1] = ro * u;
		initValues[2] = ro * v;
		initValues[3] = ro * w;
		initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0); //TO DO check
		return initValues;
	};
};

void runImpactShockTest(int argc, char *argv[]) {
	Kernel _kernel;
	const double PI = 3.14159265359;
	
	_kernel.Initilize(&argc, &argv);
	double L = 15e-2; // 5 cm; 
	Grid _grid = GenGrid1D(_kernel.getParallelHelper(), 1800, -L, L, false); //Change grid size here
	_kernel.BindGrid(_grid);
	//_kernel.LoadGrid("C:\\Users\\Ilya\\Downloads\\cilindr5.11(1).cgns");

	//Fill in configuration
	_kernel.ReadConfiguration(""); //Change run parameters here
	_kernel.InitCalculation();

	//Initial conditions
	//double ro0 = 1;
	double ro0 = 1000 * 1.0 / 0.88200003E-01; // SI	for Pb
	//double ro0 = 1000 * 1.0 / 0.127; // SI	for stainless steel
	double V = 500; //m/s
	ImpactShockInitialConditions ic(V, ro0);
	_kernel.GenerateInitialConditions(ic);	

	//Run test
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();

	//Output result
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	_kernel.Finalize();	
};

//Kelvin-Helmgolz instability
//Shear layer computation
class KelvinHelmholzInitialConditions : public InitialConditions::InitialConditions
{
private:
	double _Pressure;
	double _V;
	double _roUp;
	double _roDown;
	double _LX;
	double _LY;
	double _A;
	int _nModes;
public:	
	//Constructor
	KelvinHelmholzInitialConditions(double Pressure, double V, double roUp, double roDown, double LX, double LY, double A, int nModes) {
		_Pressure = Pressure;
		_V = V;
		_roUp = roUp;	
		_roDown = roDown;
		_LX = LX;
		_LY = LY;
		_A = A;
		_nModes = nModes;
	};

	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;				
		//Other velocities
		double v = 0;
		double w = 0;

		//Internal energy
		double e = 10;

		//Upper state		
		double roUp = _roUp;		
		double uUp = 0.5*_V;
			
		//Right state		
		double roDown = _roDown;		
		double uDown = -0.5*_V;

		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
				
		//Values
		double u = 0;
		double ro = 0;
		double roE = 0;			
		if (y <= 0) {
			ro = roDown;
			u = uDown;			
		} else {
			ro = roUp;
			u = uUp;						
		};
		
			
		//Convert to conservative variables
		roE = ro*(e + (u*u + v*v + w*w) / 2.0);
		initValues.resize(_gasModel->nConservativeVariables);
		initValues[0] = ro;
		initValues[1] = ro * u;
		initValues[2] = ro * v;
		initValues[3] = ro * w;
		initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0); //TO DO check
		return initValues;
	};
};

void runKelvinHelmholzTest(int argc, char *argv[]) {

};

//Main program ))
 int main(int argc, char *argv[]) {
	//runPeriodicTest2D(argc, argv);
	//runSodTest(argc, argv);
	//runShearLayer(argc, argv);
 	//runImpactShockTest(argc, argv);
	//LomonosovFortovGasModel gasModel(1);
	//double P;
	//double c;
	//bool nonPhysical;
	//double ro = 1./0.88200003E-01;
	//double e = 0;	
	//double Gr = 0;
	//gasModel.EOSE5(ro, e, P, c, Gr, nonPhysical);
	//double P_SI = 1e5;
	//double ro_SI = ro * 1000;
	////e = gasModel.FindInternalEnergy(ro_SI, P_SI);
	//std::cout<<"Pressure = "<<P<<"\n";
	////std::cout<<"Energy = "<<e<<"\n";
	return 0;	
};
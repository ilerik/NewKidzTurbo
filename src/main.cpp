#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "kernel.h"
#include "gridgeneration.h"
#include <cmath>
#include <random>
#include "MeshMovement.h"

#include "cuda.h"
#include <cuda_runtime.h>

//Test cases
#include "test_list.h"


template< typename T >
std::string int_to_hex( T i )
{
  std::stringstream stream;
  stream << "0x" 
         << std::setfill ('0') << std::setw(sizeof(T)*2) 
         << std::hex << i;
  return stream.str();
}


//Shear layer computation
class ImpactShockInitialConditions : public InitialConditions::InitialConditions
{
private:		
	double roL_;
	double roR_;
	int nmatLeft_;
	int nmatRight_;
  double uL_;
  double uR_;
  double PL_;
  double PR_;			
public:	
  ImpactShockInitialConditions(int nmatLeft, double roL, double uL, double PL, int nmatRight, double roR, double uR, double PR) :
    nmatLeft_{ nmatLeft },
    nmatRight_{ nmatRight },
    roL_{ roL },
    roR_{ roR },
    uL_{ uL },
    uR_{ uR },
    PL_{ PL },
    PR_{ PR }
  {	};

	virtual int getInitialGasModelIndex(const Cell& cell) {
		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
		
		if (x <= 0) {
			return nmatLeft_; 
		} else {
			return nmatRight_; 
		};
	};

	virtual std::vector<double> getInitialValues(const Cell& cell) {
		int nmat = getInitialGasModelIndex(cell); //get material index
		std::vector<double> initValues;				

		//Other velocities
		double v = 0;
		double w = 0;

		//Internal energy
		double e = 0;

		//Cell center
		double x = cell.CellCenter.x;
		double y = cell.CellCenter.y;
		double z = cell.CellCenter.z;
				
		//Values
		double u = 0;
		double ro = 0;
    double P = 0;
		double roE = 0;				
		if (x <= 0) {
      P = PL_;
			ro = roL_;
			u = uL_;			
		} else {
      P = PR_;
			ro = roR_;
			u = uR_;						
		};
			
		//Convert to conservative variables		    
    double gamma = 5.0 / 3.0;
		e = P / ((gamma - 1.0) * ro);
		roE = ro*(e + (u*u + v*v + w*w) / 2.0);
		initValues.resize(_gasModels[nmat]->nConservativeVariables);
		initValues[0] = ro;
		initValues[1] = ro * u;
		initValues[2] = ro * v;
		initValues[3] = ro * w;
		initValues[4] = roE;
		return initValues;
	};
};

void runImpactShockTest(int argc, char *argv[]) {
	Kernel _kernel;
	const double PI = 3.14159265359;
	
	_kernel.Initilize(&argc, &argv);
	double L = 70e-3; // 70 mm; 
	Grid _grid = GenGrid1D(_kernel.getParallelHelper(), 500, -L, L, false); //Change grid size here
	_kernel.BindGrid(&_grid);
	//_kernel.LoadGrid("C:\\Users\\Ilya\\Downloads\\cilindr5.11(1).cgns");

	//Fill in configuration
	_kernel.ReadConfiguration(""); //Change run parameters here
	_kernel.InitCalculation();

	//Initial conditions
	//double ro0 = 1;
  double gamma = 5.0 / 3.0;
	double roXe = 5.894;  //kg / m^3
  double roAr = 1.784;   //kg / m^3
	double pXe = 50000;    //0.5 bar = 0.05 Pa
  double pAr = 50000;    //0.5 bar = 0.05 Pa	
  double uXe = 0.0;     //m/s
  double uAr = 0.0;     //m/s	
  
  //After shock
  double D = -600;       //shock speed
  double uArAS = 2 * (D*D*roAr - gamma * pAr) / (D * roAr * (1.0 + gamma));
  double roArAS = D * roAr / (D - uArAS);
  double pArAS = pAr + D * roAr * uArAS;

  D = 600;       //shock speed
  double uXeAS = 2 * (D*D*roXe - gamma * pXe) / (D * roXe * (1.0 + gamma));
  double roXeAS = D * roXe / (D - uXeAS);
  double pXeAS = pXe + D * roXe * uXeAS;
  ImpactShockInitialConditions ic(0, roXeAS, uXeAS, pXeAS, 1, roAr, uAr, pAr);
  //ImpactShockInitialConditions ic(0, roXeAS, uXeAS, pXeAS, 1, roXe, uXe, pXe);
  //ImpactShockInitialConditions ic(0, roXe, uXe, pXe, 1, roArAS, uArAS, pArAS);  
	//ImpactShockInitialConditions ic(1, roAr, uAr, pAr, 1, roArAS, uArAS, pArAS);
	_kernel.GenerateInitialConditions(&ic);	

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
  runImpactShockTest(argc, argv);
	return 0;

};

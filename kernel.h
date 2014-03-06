#include "grid.h"

////Calculation kernel
class Kernel {
public:
	//Grid data
	Grid _grid;

	//Task parameters
	int RungeKuttaOrder;	

	//Internal storage	
	double* Uconservative;	//Conservative variables
	double* Residual;		//Residual

	//Compute residual
	void ComputeResidual(double *R, double *U){

	};

	//Compute convective flux through face
	void ComputeConvectiveFlux() {
	};
	

	//Explicit time step
	void ExplicitTimeStep() {
	};

	//Implicit time step
	void ImplicitTimeStep() {
	};
};
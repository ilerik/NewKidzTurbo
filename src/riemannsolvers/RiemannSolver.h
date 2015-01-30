#ifndef TURBO_RIEMANNSOLVERS_RiemannSolver
#define TURBO_RIEMANNSOLVERS_RiemannSolver

//Base class for riemann solvers

#include "datatypes.h"
#include "basetypes.h"
#include "grid.h"
#include "RiemannSolverConfiguration.h"
#include "GasModel.h"

//Solution information structure
class RiemannProblemSolutionResult {
public:
	std::vector<double> Fluxes;
	std::vector<double> FluxesLeft; //Conservative flux to the left state
	std::vector<double> FluxesRight; //Conservative flux to the right state
	double MaxEigenvalue; //Maximal local eigenvalue of Jacobian
	Vector Velocity; //Interface velocity estimate
	double Pressure; //Interface pressure estimate
};

//Base class for all riemann solvers
class RiemannSolver {	
protected:
	std::vector<GasModel*> _gasModels; //Link to gasmodel to use
public:	
	//Check gas models and bind them to solver
	virtual bool BindGasModels(std::vector<GasModel*>& gasModels) = 0;

	//Load settings from configuration object
	virtual bool loadConfiguration(Logger* logger, RiemannSolverConfiguration configuration) = 0;

	//Solve riemann problem	
	virtual RiemannProblemSolutionResult Solve(int nmatL, const GasModel::ConservativeVariables& UL, int nmatR, const GasModel::ConservativeVariables& UR, const Face& f, Vector faceVelocity) = 0;
};

#endif
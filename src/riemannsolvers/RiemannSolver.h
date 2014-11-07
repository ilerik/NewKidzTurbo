#ifndef TURBO_RIEMANNSOLVERS_RiemannSolver
#define TURBO_RIEMANNSOLVERS_RiemannSolver

//Base class for riemann solvers

#include "datatypes.h"
#include "basetypes.h"
#include "grid.h"
#include "GasModel.h"

//Solution information structure
class RiemannProblemSolutionResult {
public:
	std::vector<double> Fluxes; //Conservative flux
	double MaxEigenvalue; //Maximal local eigenvalue of Jacobian
	Vector Velocity; //Interface velocity estimate	
};

//Base class for all riemann solvers
class RiemannSolver {	
	GasModel* _gasModel; //Link to gasmodel to use

	//Private helpers
	inline double takeRoeAverage(double& roLRoot, double roRRoot, double fL, double fR) {
		double roeAvg = (roLRoot * fL + roRRoot * fR) / (roLRoot + roRRoot);
		return roeAvg;
	};
public:	
	RiemannSolver(GasModel* gasModel) : _gasModel(gasModel) { };	

	//Solve riemann problem	
	RiemannProblemSolutionResult Solve(const GasModel::ConservativeVariables& UL, const GasModel::ConservativeVariables& UR, const Face& f) {
		RiemannProblemSolutionResult result;
		result.Fluxes.resize(_gasModel->nConservativeVariables, 0);

		//The Harten, Lax, and van Leer with contact restoration (HLLC) Riemann solver
		//Left and right states
		double roL = UL.ro;
		Vector vL = Vector(UL.rou / roL, UL.rov / roL, UL.row / roL);
		double uL = vL * f.FaceNormal;
		double pL = 0;
		double cL = 0;
		double GrL = 0;		
		_gasModel->GetPressureAndSoundSpeed(UL, pL, cL, GrL);
		double phiL = cL*cL - GrL * pL / roL;		

		double roR = UR.ro;
		Vector vR = Vector(UR.rou / roR, UR.rov / roR, UR.row / roR);
		double uR = vR * f.FaceNormal;
		double pR = 0;
		double cR = 0;
		double GrR = 0;
		double GammaR = 0;		
		_gasModel->GetPressureAndSoundSpeed(UR, pR, cR, GrR);			
		double phiR = cR*cR - GrR * pR / roR;

		//Generalized Roe averages (according to Hu et al)
		double roLRoot = sqrt(roL);
		double roRRoot = sqrt(roR);
		double roRoe = roLRoot * roRRoot; //Roe averaged density
		double uRoe = takeRoeAverage(roLRoot, roRRoot, uL, uR); //Roe averaged velocity (Hu et al (17))
		double pOverRoRoe = takeRoeAverage(roLRoot, roRRoot, pL/roL, pR/roR) + 0.5*pow((uR - uL) / (roRRoot + roLRoot), 2.0); //(Hu et al (18))		
		double phiRoe = takeRoeAverage(roLRoot, roRRoot, phiL, phiR); //Roe averaged phi (Hu et al (25))
		double GrRoe = takeRoeAverage(roLRoot, roRRoot, GrL, GrR); //Roe averaged Gruneisen coefficient (Hu et al (25))
		double cRoe = phiRoe + GrRoe * pOverRoRoe; //Roe averaged sound speed (Hu et al (16))
		cRoe = sqrt(cRoe);

		//Wave speeds for two waves (Hu et al (12))
		double SL = min(uL - cL, uRoe - cRoe); //bl
		double SR = max(uR + cR, uRoe + cRoe); //br

		//The HLLC Flux for Euler equations (Toro 10.4.2, p324)
		double SStar = (pR - pL) + roL*uL*(SL - uL) - roR*uR*(SR - uR); //Toro 10.37
		SStar /=  roL*(SL - uL) - roR*(SR - uR); // Speed of middle wave
		double PStarL = pL + roL * (SL - uL)*(SStar - uL); //Pressure estimate for left star region Toro 10.36
		double PStarR = pR + roR * (SR - uR)*(SStar - uR); //Pressure estimate for right star region Toro 10.36
		
		//Choose flux according to wave speeds pattern Toto 10.26
		//And choose estimate for interface velocity
		Vector velocity;
		std::vector<double> DL(5);
		DL[0] = 0;
		DL[1] = 1;
		DL[2] = 0;
		DL[3] = 0;
		DL[4] = uL;
		std::vector<double> FL = uL*UL + pL*DL;
		if (0 <= SL) {
			result.Fluxes = FL;
		};

		std::vector<double> DR(5);
		DR[0] = 0;
		DR[1] = f.FaceNormal.x;
		DR[2] = f.FaceNormal.y;
		DR[3] = f.FaceNormal.z;
		DR[4] = uR;
		std::vector<double> FR = uR*UR + pR*DR;
		if ( 0 >= SR ) {
			result.Fluxes = FR;
		};

		//Toro formulation 10.41 for HLLC flux 
		std::vector<double> DStar(5);
		DStar[0] = 0;
		DStar[1] = f.FaceNormal.x;
		DStar[2] = f.FaceNormal.y;
		DStar[3] = f.FaceNormal.z;
		DStar[4] = SStar;		
		if (( SL < 0 ) && ( 0 <= SStar )) {
			result.Fluxes = SStar * (SL * UL - FL) + SL * (pL + roL*(SL - uL)*(SStar - uL)) * DStar;
			result.Fluxes /= SL - SStar;
		};
		if (( SStar < 0 ) && ( 0 < SR )) {
			result.Fluxes = SStar * (SR * UR - FR) + SR * (pR + roR*(SR - uR)*(SStar - uR)) * DStar;
			result.Fluxes /= SR - SStar;
		};

		//Estimate for maximal speed
		result.MaxEigenvalue = max(abs(SL), abs(SR));

		//Estimate for velocity
		result.Velocity = velocity;

		//DEBUG
		/*Roe3DSolverPerfectGas rSolver;
		rSolver.SetGamma(1.4);
		rSolver.SetOperatingPressure(0.0);
		rSolver.SetHartenEps(0.0);
		std::vector<double> roeFlux = rSolver.ComputeFlux(UL, UR, f);*/
		return result;
	};
};

#endif
#ifndef TURBO_RIEMANNSOLVERS_GENERAL_EOS_HLLC_3D
#define TURBO_RIEMANNSOLVERS_GENERAL_EOS_HLLC_3D

#include "datatypes.h"
#include "basetypes.h"
#include "grid.h"
#include "GasModel.h"
#include "RiemannSolver.h"

class HLLCSolverGeneralEOS : public RiemannSolver {
	//Private helpers
	inline double takeRoeAverage(double& roLRoot, double roRRoot, double fL, double fR) {
		double roeAvg = (roLRoot * fL + roRRoot * fR) / (roLRoot + roRRoot);
		return roeAvg;
	};
public:	
	using RiemannSolver::RiemannSolver; //Inherit constructor

	//Check gas models
	virtual bool BindGasModels(std::vector<GasModel*>& gasModels) {
		bool isCorrect = true;
		for (GasModel* gasModel : gasModels) {
			//Check if number of conservative variables is equal to 5
			if (gasModel->nConservativeVariables != 5) {
				isCorrect = false;
				break;
			};
		};

		//If everything is ok bind gas models
		if (isCorrect) {
			_gasModels = gasModels;
		};
		return isCorrect;
	};

	//Load settings from configuration object
	bool loadConfiguration(Logger* logger, RiemannSolverConfiguration configuration) {
		return true;
	};

	//Solve riemann problem	
	RiemannProblemSolutionResult Solve(int nmatL, const GasModel::ConservativeVariables& UL, int nmatR, const GasModel::ConservativeVariables& UR, const Face& f, double ALEindicator) {
		RiemannProblemSolutionResult result;
		result.FluxesLeft.resize(_gasModels[nmatL]->nConservativeVariables, 0);
		result.FluxesRight.resize(_gasModels[nmatR]->nConservativeVariables, 0);

		//The Harten, Lax, and van Leer with contact restoration (HLLC) Riemann solver
		//Left and right states
		double roL = UL.ro;
		Vector vL = Vector(UL.rou / roL, UL.rov / roL, UL.row / roL);
		double uL = vL * f.FaceNormal;
		double pL = 0;
		double cL = 0;
		double GrL = 0;	
		assert(roL > 0);
		_gasModels[nmatL]->GetPressureAndSoundSpeed(UL, pL, cL, GrL);				
		double phiL = cL*cL - GrL * pL / roL;

		double roR = UR.ro;
		Vector vR = Vector(UR.rou / roR, UR.rov / roR, UR.row / roR);
		double uR = vR * f.FaceNormal;
		double pR = 0;
		double cR = 0;
		double GrR = 0;
		double GammaR = 0;				
		assert(roR > 0);
		_gasModels[nmatR]->GetPressureAndSoundSpeed(UR, pR, cR, GrR);			
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

		//Sample state and pressure
		std::vector<double> D(5);
		D[0] = 0;
		D[1] = f.FaceNormal.x;
		D[2] = f.FaceNormal.y;
		D[3] = f.FaceNormal.z;
		D[4] = 0;
		double u = 0;
		GasModel::ConservativeVariables U;
		D[4] = uL;
		std::vector<double> FL = uL*UL + pL*D;
		if (0 <= SL) {
			U = UL;
			result.Pressure = pL; //Pressure estimate
			u = uL; 
		};
	
		D[4] = uR;
		std::vector<double> FR = uR*UR + pR*D;
		if ( 0 >= SR ) {
			U = UR;
			result.Pressure = pR; //Pressure estimate
			u = uR;
		};

		//Toro formulation 10.41 for HLLC flux 
		if (( SL < 0 ) && ( 0 <= SStar )) {			
			result.Pressure = pL + SL * roL*(SL - uL)*(SStar - uL); //Pressure estimate
			U = (SStar * (SL * UL - FL)) / (SL - SStar);
			u = SStar;
			/*result.FluxesLeft = SStar * (SL * UL - FL);
			result.FluxesLeft = SL * (pL + roL*(SL - uL)*(SStar - uL)) * DStar;
			result.FluxesLeft /= SL - SStar;
			result.FluxesRight = result.FluxesLeft;*/
		};
		if (( SStar < 0 ) && ( 0 < SR )) {
		/*	result.FluxesLeft = SStar * (SR * UR - FR) + SR * (pR + roR*(SR - uR)*(SStar - uR)) * DStar;
			result.FluxesLeft /= SR - SStar;
			result.FluxesRight = result.FluxesLeft;	*/                                  		
			result.Pressure = pR + SR * roR * (SR - uR) * (SStar - uR); //Pressure estimate
			U = (SStar * (SR * UR - FR)) / (SR - SStar);
			u = SStar;
		};

		//And choose estimate for interface velocity
		Vector velocity = ALEindicator * u * f.FaceNormal; 
		double uRelative = (u * (1.0 - ALEindicator));
		assert(velocity.x == velocity.x);
	
		//Now compute flux
		D[4] = uRelative;
		result.Fluxes = (uRelative)*U + result.Pressure*D;

		//Estimate for maximal speed
		result.MaxEigenvalue = max(abs(SL), abs(SR));

		//Estimate for velocity
		result.Velocity = velocity;		

		return result;
	};
};

#endif
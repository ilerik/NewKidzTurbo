#ifndef TURBO_RIEMANNSOLVERS_GENERAL_EOS_HLLC_3D
#define TURBO_RIEMANNSOLVERS_GENERAL_EOS_HLLC_3D

#include "datatypes.h"
#include "basetypes.h"
#include "utilityfunctions.h"
#include "grid.h"
#include "GasModel.h"
#include "RiemannSolver.h"

class HLLCSolverGeneralEOS : public RiemannSolver {	
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
	RiemannProblemSolutionResult Solve(int nmatL, const GasModel::ConservativeVariables& UL, int nmatR, const GasModel::ConservativeVariables& UR, const Face& f, Vector faceVelocity) {
		RiemannProblemSolutionResult result;
		result.FluxesLeft.resize(_gasModels[nmatL]->nConservativeVariables, 0);
		result.FluxesRight.resize(_gasModels[nmatR]->nConservativeVariables, 0);

		//The Harten, Lax, and van Leer with contact restoration (HLLC) Riemann solver
		//Left and right states
		double roL = UL.ro;
		Vector vL = Vector(UL.rou / roL, UL.rov / roL, UL.row / roL);
		double uL = (vL - faceVelocity) * f.FaceNormal;
		double pL = 0;
		double cL = 0;
		double GrL = 0;	
		assert(roL > 0);
		_gasModels[nmatL]->GetPressureAndSoundSpeed(UL, pL, cL, GrL);				
		double phiL = cL*cL - GrL * pL / roL;

		double roR = UR.ro;
		Vector vR = Vector(UR.rou / roR, UR.rov / roR, UR.row / roR);
		double uR = (vR - faceVelocity) * f.FaceNormal;
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
		double uRoe = takeRoeAverage(roLRoot, roRRoot, vL * f.FaceNormal, vR * f.FaceNormal) - faceVelocity * f.FaceNormal; //Roe averaged velocity (Hu et al (17))
		double pOverRoRoe = takeRoeAverage(roLRoot, roRRoot, pL/roL, pR/roR) + 0.5*pow((uR - uL) / (roRRoot + roLRoot), 2.0); //(Hu et al (18))		
		double phiRoe = takeRoeAverage(roLRoot, roRRoot, phiL, phiR); //Roe averaged phi (Hu et al (25))
		double GrRoe = takeRoeAverage(roLRoot, roRRoot, GrL, GrR); //Roe averaged Gruneisen coefficient (Hu et al (25))
		double cRoe = phiRoe + GrRoe * pOverRoRoe; //Roe averaged sound speed (Hu et al (16))
		//cRoe = sqrt(cRoe);
		cRoe = takeRoeAverage(roLRoot, roRRoot, cL, cR); //TO DO check
		cRoe = sqrt(cRoe);
		assert(cRoe > 0);

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

		if (0 <= SL) {
			std::vector<double> DL(5);
			DL[0] = 0;
			DL[1] = f.FaceNormal.x;
			DL[2] = f.FaceNormal.y;
			DL[3] = f.FaceNormal.z;
			DL[4] = vL*f.FaceNormal;
			std::vector<double> FL = uL*UL + pL*DL;
			result.Fluxes = FL;
		};

		//
		if ( 0 >= SR ) {
			std::vector<double> DR(5);
			DR[0] = 0;
			DR[1] = f.FaceNormal.x;
			DR[2] = f.FaceNormal.y;
			DR[3] = f.FaceNormal.z;
			DR[4] = vR * f.FaceNormal;
			std::vector<double> FR = uR*UR + pR*DR;
			result.Fluxes = FR;
		};

		//Toro formulation 10.41 for HLLC flux 
		std::vector<double> DStar(5);
		DStar[0] = 0;
		DStar[1] = f.FaceNormal.x;
		DStar[2] = f.FaceNormal.y;
		DStar[3] = f.FaceNormal.z;
		DStar[4] = SStar + faceVelocity * f.FaceNormal;
		
 		if (( SL < 0 ) && ( 0 <= SStar )) {
			double pLStar = pL + roL*(SL - uL)*(SStar - uL);
			std::vector<double> ULStar(5);
			ULStar = UL * (SL - uL);
			ULStar[1] += (pLStar - pL)*f.FaceNormal.x;
			ULStar[2] += (pLStar - pL)*f.FaceNormal.y;
			ULStar[3] += (pLStar - pL)*f.FaceNormal.z;
			ULStar[4] += -pL*uL + pLStar*SStar;
			ULStar /= (SL - SStar);
			
			result.Fluxes = SStar * ULStar + pLStar * DStar;
			result.Pressure = pLStar;
		};
		if (( SStar < 0 ) && ( 0 < SR )) {
			double pRStar = pR + roR*(SR - uR)*(SStar - uR);
			std::vector<double> URStar(5);
			URStar = UR * (SR - uR);
			URStar[1] += (pRStar - pR)*f.FaceNormal.x;
			URStar[2] += (pRStar - pR)*f.FaceNormal.y;
			URStar[3] += (pRStar - pR)*f.FaceNormal.z;
			URStar[4] += -pR*uR + pRStar*SStar;
			URStar /= (SR - SStar);
			
			result.Fluxes = SStar * URStar + pRStar * DStar;
			result.Pressure = pRStar;
		};

		//Estimate for maximal speed
		result.MaxEigenvalue = max(abs(SL), abs(SR));

		return result;
	};
};

#endif
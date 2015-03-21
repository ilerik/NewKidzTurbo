#ifndef TURBO_RIEMANNSOLVERS_PERFECT_GAS_GODUNOV_3D
#define TURBO_RIEMANNSOLVERS_PERFECT_GAS_GODUNOV_3D

#include "datatypes.h"
#include "basetypes.h"
#include "utilityfunctions.h"
#include "grid.h"
#include "GasModel.h"
#include "RiemannSolver.h"

class Godunov3DSolverPerfectGas : public RiemannSolver {
	//Required info
	double gamma;
public:
	using RiemannSolver::RiemannSolver; //Inherit constructor

	//Check gas models
	virtual bool BindGasModels(std::vector<GasModel*>& gasModels) {
		bool isCorrect = true;

		//Pointer to perfect gas model
		PerfectGasModel* perfectGasModel;
		std::vector<double> gammas;

		//Only several perfect gas models with common value of gamma allowed
		for (GasModel* gasModel : gasModels) {
			//Check if it's perfect gas 
			if (perfectGasModel = dynamic_cast<PerfectGasModel*>(gasModel)) {
				//and save gamma value
				gammas.push_back(perfectGasModel->Gamma);
			} else {
				//It's not perfect gas model
				return false;
			};
		};

		//Check if all gammas are the same
		for (int i = 0; i<gammas.size() - 1; i++) {
			if (gammas[i] != gammas[i-1]) {
				isCorrect = false;
				break;
			};
		};

		//If everything is ok bind gas models
		if (isCorrect) {
			_gasModels = gasModels;
			//And set gamma
			gamma = gammas[0];
		};
		return isCorrect;
	};

	//Load settings from configuration object
	bool loadConfiguration(Logger* logger, RiemannSolverConfiguration configuration) {
		//Load properties
		std::pair<double, bool> res;
		return true;
	};

	//empty constructor
	Godunov3DSolverPerfectGas(){
	};

	//case without vacuum
	double GetPressure(const GasModel::ConservativeVariables& celldata) {
		double ro = celldata.ro;
		double vx = celldata.rou/celldata.ro;
		double vy = celldata.rov/celldata.ro;
		double vz = celldata.row/celldata.ro;
		double E = celldata.roE/celldata.ro;
		double P = (gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
		return P;
	};

	//Numerical flux
	std::vector<double> F(GasModel::ConservativeVariables U, Vector n)
	{		
		std::vector<double> res(5,0);		
		double ro = U.ro;
		double vx = U.rou/ro;
		double vy = U.rov/ro;
		double vz = U.row/ro;
		double roe = U.roE;	//ro*e
		double p = (gamma-1.0)*(roe-ro*(vx*vx+vy*vy+vz*vz)/2.0);		
		res[0] = n.x * (ro*vx) + n.y*(ro*vy) + n.z*(ro*vz);
		double vn = vx*n.x + vy*n.y + vz*n.z;
		res[1] = n.x * (ro*vx*vx+p) + n.y*(ro*vx*vy) + n.z*(ro*vx*vz);
		res[2] = n.x * (ro*vy*vx) + n.y*(ro*vy*vy+p) + n.z*(ro*vy*vz);
		res[3] = n.x * (ro*vz*vx) + n.y*(ro*vz*vy) + n.z*(ro*vz*vz+p);
		res[4] = (n.x * vx + n.y * vy + n.z * vz)*(roe+p);		
		return res;
	};

	struct fParameters {
		double roL;
		double pL;
		double uL;
		double roR;
		double pR;
		double uR;
		double deltaU;		
	} params;

	//fL part of algebraic equation for pressure in exact RP solver		(proposition 4.2.1 from Toro)
	double f1(double roL, double uL, double pL, double pStar) 
	{
		double res = 0;
		
		if (pStar > pL) {
			//Left shock
			double AL = 2 / ((gamma + 1.0) * roL);
			double BL = pL * (gamma - 1.0) / (gamma + 1.0);
			res = (pStar - pL) * sqrt(AL / (pStar + BL));
		} else {
			//Left rarefaction
			double aL = sqrt(gamma * pL / roL);
			res = pow(pStar/pL, (gamma - 1.0)/(2*gamma)) - 1.0;
			res *= 2*aL/(gamma - 1.0);
		};

		return res;
	};

	//fR part
	double f2(double roR, double uR, double pR, double pStar) 
	{
		double res = 0;

		if (pStar > pR) {
			//Right shock
			double AR = 2 / ((gamma + 1.0) * roR);
			double BR = pR * (gamma - 1.0) / (gamma + 1.0);
			res = (pStar - pR) * sqrt(AR / (pStar + BR));
		} else {
			//Right rarefaction
			double aR = sqrt(gamma * pR / roR);
			res = pow(pStar/pR, (gamma - 1.0)/(2*gamma)) - 1.0;
			res *= 2*aR/(gamma - 1.0);
		};

		return res;
	};

	//Target function for pressure equation in exact RP solver
	static void fFunction(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *object)
	{
		Godunov3DSolverPerfectGas* solver = (Godunov3DSolverPerfectGas*)object;
		//
		// this callback calculates		
		// f(pStar) = f1(pStar, WL) + f2(pStar, WR) + deltaU
		//
		double pStar = x[0];
		fParameters par = solver->params;
		double res = solver->f1(par.roL, par.uL, par.pL, pStar) + solver->f2(par.roR, par.uR, par.pR, pStar) + par.deltaU;

		fi[0] = res;		
	}

	//Types of nonlinear waves
	enum WaveType {
		Shock,
		Rarefaction
	};

	//Find star region quantites
	struct StarVariables {
		double pStar;
		double uStar;
		double roStarL;
		double roStarR;
		WaveType leftWave;
		WaveType rightWave;

		//Propagation speeds
		double SL;
		double SHL;
		double STL;

		double SR;
		double SHR;
		double STR;

		double MaxSpeed;
	};	

	StarVariables ComputeStarVariables(double roL, double uL, double pL, double roR, double uR, double pR) {				
		StarVariables result;

		//Left state
		params.roL = roL;
		params.uL = uL;
		params.pL = pL;

		//Right state
		params.roR = roR;
		params.uR = uR;
		params.pR = pR;		

		//Normal velocity difference
		params.deltaU = params.uR - params.uL;

		//Compute star region pressure
		alglib::real_1d_array x;
		alglib::real_1d_array bl;
		alglib::real_1d_array bu;
		x.setlength(1);		
		bl.setlength(1);
		bu.setlength(1);

		//Set boundary constraints (non negative pressure)
		bl[0] = 0.0;
		bu[0] = std::numeric_limits<double>::max();

		//Initial guess (possible othe choises, the simplest for now)
		double pStar = 0.5*(params.pL + params.pR);
		x[0] = pStar;

		//Iterative scheme parameters
		double epsg = 1e-10;
		double epsf = 0;
		double epsx = 0;
		double diffstep = 1e-6;
		alglib::ae_int_t maxits = 0;
		alglib::minlmstate state;
		alglib::minlmreport rep;

		alglib::minlmcreatev(1, x, diffstep, state);
		alglib::minlmsetbc(state, bl, bu);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);				
		alglib::minlmoptimize(state, fFunction, NULL, (void*)this);		
		alglib::minlmresults(state, x, rep);

		//Check if solution converged (TO DO)		
		result.pStar = x[0];		

		//Compute star region velocity
		result.uStar = 0.5*(uL + uR) + 0.5*(f2(roR, uR, pR, result.pStar) - f1(roL, uL, pL, result.pStar));

		//Determine nonlinear wave type and properties
		double C = (gamma - 1.0) / (gamma + 1.0);
		result.MaxSpeed = 0;

		//Left side of contact
		double pRatioL = result.pStar / pL;
		double aL = sqrt(gamma * pL / roL);
		if (result.pStar > pL) {
			//Left shock
			result.leftWave = Shock;

			//Determine density in star region			
			result.roStarL = pRatioL + C;
			result.roStarL /= C * pRatioL + 1.0;
			result.roStarL *= roL;

			//Determine shock propagation speed			
			result.SL = sqrt((gamma + 1) * pRatioL / (2*gamma) + (gamma - 1) / (2*gamma));
			result.SL = uL - aL * result.SL;

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SL));
		} else {
			//Left rarefaction
			result.leftWave = Rarefaction;

			//Determine density in star region
			result.roStarL = roL * pow(pRatioL , 1.0 / gamma);

			//Determine rarefaction head propagation speed		
			result.SHL = uL - aL;

			//Determine rarefaction tail propagation speed		
			double aStarL = aL * pow(pRatioL, (gamma - 1) / (2*gamma));
			result.STL = result.uStar - aStarL;		

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SHL));
		};

		//Right side of contact
		double pRatioR = result.pStar / pR;
		double aR = sqrt(gamma * pR / roR);
		if (result.pStar > pR) {
			//Right shock
			result.rightWave = Shock;

			//Determine density in star region			
			result.roStarR = pRatioR + C;
			result.roStarR /= C * pRatioR + 1.0;
			result.roStarR *= roR;

			//Determine shock propagation speed			
			result.SR = sqrt((gamma + 1) * pRatioR / (2*gamma) + (gamma - 1) / (2*gamma));
			result.SR = uL + aR * result.SR;

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SR));
		} else {
			//Left rarefaction
			result.rightWave = Rarefaction;

			//Determine density in star region
			result.roStarR = roR * pow(pRatioR , 1.0 / gamma);

			//Determine rarefaction head propagation speed		
			result.SHR = uL + aR;

			//Determine rarefaction tail propagation speed		
			double aStarR = aR * pow(pRatioL, (gamma - 1) / (2*gamma));
			result.STR = result.uStar + aStarR;			

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SHR));
		};
				
		return result;
	};
	  	
	//Solve riemann problem
	std::vector<double> ComputeFlux(const GasModel::ConservativeVariables& UL, const GasModel::ConservativeVariables& UR, const Face& f) {
		std::vector<double> res(5,0);		
		Vector velocityL = Vector(UL.rou, UL.rov, UL.row) / UL.ro;
		Vector velocityR = Vector(UR.rou, UR.rov, UR.row) / UR.ro;

		//Compute normal velocity
		double uL = velocityL * f.FaceNormal;
		double uR = velocityR * f.FaceNormal;

		//Left state
		double roL = UL.ro;
		double pL = GetPressure(UL);

		//Right state
		double roR = UR.ro;
		double pR = GetPressure(UR);			

		//Solve riemann problem
		StarVariables starValues = ComputeStarVariables(roL, uL, pL, roR, uR, pR);
		
		//Compute maximum propagation speed		
		MaxEigenvalue = starValues.MaxSpeed;

		//Compute flux (Toro p. 219) 
		//Sample exact solution at S = x/t = 0
		double S = 0;
		double ro;
		double u;
		double p;
		
		if (starValues.uStar >= S) {
			//Left side of contact

			//Shock wave
			if (starValues.leftWave == Godunov3DSolverPerfectGas::Shock) {
				if (starValues.SL >= S) {
					//Case a1
					//Supersonic shock
					ro = roL;
					u = uL;
					p = pL;
				} else {
					//Case a2
					//Subsonic shock
					ro = starValues.roStarL;
					u = starValues.uStar;
					p = starValues.pStar;
				};
			};

			//Rarefaction wave
			if (starValues.leftWave == Godunov3DSolverPerfectGas::Rarefaction) {
				if (starValues.SHL > S) {
					//Left region
					ro = roL;
					u = uL;
					p = pL;
				} else if ( S > starValues.STL) {
					//Star region
					ro = starValues.roStarL;
					u = starValues.uStar;
					p = starValues.pStar;
				} else {
					//Rarefaction fan region
					double aL = sqrt(gamma * pL / roL);
					double C = 2 / (gamma + 1) + (gamma-1) * (uL - S) / ((gamma + 1) * aL);
					C = pow(C, 2 / (gamma - 1));

					//Density
					ro = C * roL;

					//Velocity
					u = aL + (gamma - 1) * uL/ 2.0 + S;
					u *= 2 / (gamma + 1);

					//Pressure					
					p = C * pL;
				};
			};

		} else {
			//Right side of contact

			//Shock wave
			if (starValues.rightWave == Godunov3DSolverPerfectGas::Shock) {
				if (starValues.SR <= S) {
					//Case a1
					//Supersonic shock
					ro = roR;
					u = uR;
					p = pR;
				} else {
					//Case a2
					//Subsonic shock
					ro = starValues.roStarR;
					u = starValues.uStar;
					p = starValues.pStar;
				};
			};

			//Rarefaction wave
			if (starValues.rightWave == Godunov3DSolverPerfectGas::Rarefaction) {
				if (starValues.SHR < S) {
					//Right region
					ro = roR;
					u = uR;
					p = pR;
				} else if ( S < starValues.STR) {
					//Star region
					ro = starValues.roStarR;
					u = starValues.uStar;
					p = starValues.pStar;
				} else {
					//Rarefaction fan region
					double aR = sqrt(gamma * pR / roR);
					double C = 2 / (gamma + 1) - (gamma-1) * (uR - S) / ((gamma + 1) * aR);
					C = pow(C, 2 / (gamma - 1));

					//Density
					ro = C * roR;

					//Velocity
					u = -aR + (gamma - 1) * uR/ 2.0 + S;
					u *= 2 / (gamma + 1);

					//Pressure					
					p = C * pR;
				};
			};
		};

		//Compute flux given variables values at point S = 0
		GasModel::ConservativeVariables U;
		U.ro = ro;
		Vector V; 
		if (starValues.uStar > S) {
			V = velocityL - uL * f.FaceNormal + u * f.FaceNormal;			
		} else if (starValues.uStar < S) {
			V = velocityR - uR * f.FaceNormal + u * f.FaceNormal;
		} else {
			V = 0.5 * ((velocityR - uR * f.FaceNormal) + (velocityL - uL * f.FaceNormal)) + u * f.FaceNormal;
		};
		U.rou = ro * V.x;
		U.rov = ro * V.y;
		U.row = ro * V.z;
		double e = p / (ro * (gamma - 1.0));
		U.roE = (e + V*V / 2.0) * ro;
		
		return F(U, f.FaceNormal);
	};

	//Public properties
	double MaxEigenvalue;

	//Solve riemann problem	
	RiemannProblemSolutionResult Solve(int nmatL, const GasModel::ConservativeVariables& UL, int nmatR, const GasModel::ConservativeVariables& UR, const Face& f,  Vector faceVelocity) {
		RiemannProblemSolutionResult result;		

		result.Fluxes = ComputeFlux(UL, UR, f);
		result.MaxEigenvalue = MaxEigenvalue;
		double ro = result.Fluxes[0];
		double u = result.Fluxes[1]/ro;
		double v = result.Fluxes[2]/ro;
		double w = result.Fluxes[3]/ro;
		double roe = result.Fluxes[4] - ro * (u*u + v*v + w*w) / 2.0;
		result.Velocity = Vector(u,v,w);
		result.Pressure = (gamma - 1) * roe;

		return result;
	};
};

#endif
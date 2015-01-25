#ifndef TURBO_RIEMANNSOLVERS_PERFECT_GAS_ROE_3D
#define TURBO_RIEMANNSOLVERS_PERFECT_GAS_ROE_3D

#include "datatypes.h"
#include "basetypes.h"
#include "grid.h"
#include "RiemannSolver.h"
#include "GasModel.h"
#include "PerfectGasModel.h"

class Roe3DSolverPerfectGas : public RiemannSolver {
	//Internal storage //TO DO remove
	Vector velocity; //Velocity estimate
	double MaxEigenvalue; //

	//Required info
	double gamma; //Specific heat ratio (inherited from gas model)
	double eps;	//Harten entropy correction coefficient (optional, 0 by default)
	double operatingPressure; //Operating pressure (optional, 0 by default)
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

		//Load operating pressure
		res = configuration.GetPropertyValue("OperatingPressure");
		if (res.second) {
			operatingPressure = res.first;
		} else {
			//Default
			operatingPressure = 0;
			logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Roe3DSolverPerfectGas Riemann solver operating pressure is not specified, default value of 0 assumed.");
		};

		//Load Harten's epsilon
		res = configuration.GetPropertyValue("HartenEpsilon");
		if (res.second) {
			eps = res.first;
		} else {
			//Default
			eps = 0;
			logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Roe3DSolverPerfectGas Riemann solver Harten's correction epsilon is not specified, default value of 0 assumed.");
		};

		return true;
	};

	//Numerical flux
	std::vector<double> F(GasModel::ConservativeVariables U, Vector n)
	{		
		std::vector<double> res(5,0);		
		double ro = U.ro;
		double vx = U.rou/ro;
		double vy = U.rov/ro;
		double vz = U.row/ro;
		double roE = U.roE;	//ro*e
		double p = (gamma-1.0)*(roE-ro*(vx*vx+vy*vy+vz*vz)/2.0) - operatingPressure;		
		double vn = vx*n.x + vy*n.y + vz*n.z;

		/*res[0] = n.x * (ro*vx) + n.y*(ro*vy) + n.z*(ro*vz);		
		res[1] = n.x * (ro*vx*vx+p) + n.y*(ro*vx*vy) + n.z*(ro*vx*vz);
		res[2] = n.x * (ro*vy*vx) + n.y*(ro*vy*vy+p) + n.z*(ro*vy*vz);
		res[3] = n.x * (ro*vz*vx) + n.y*(ro*vz*vy) + n.z*(ro*vz*vz+p);
		res[4] = (n.x * vx + n.y * vy + n.z * vz)*(roe+p);		*/

		res[0] = ro*vn;
		res[1] = ro*vn*vx + n.x*p;
		res[2] = ro*vn*vy + n.y*p;
		res[3] = ro*vn*vz + n.z*p;
		res[4] = vn*(roE + p);
		return res;
	};

	//Solve riemann problem
	//std::vector<double> ComputeFlux(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f) {
	std::vector<double> ComputeFlux(const GasModel::ConservativeVariables& UL, const GasModel::ConservativeVariables& UR, const Face& f) {
		std::vector<double> res(5,0);						

		//Calculate symmetric flux part	
		std::vector<double> symFluxPart = 1.0*(F(UL, (1.0)*f.FaceNormal) + F(UR, (1.0)*f.FaceNormal)); //is 1.0 or 0.5 TODO	
		res += symFluxPart;

		//for (int i = 0; i<nv; i++) printf("%lg\n", res[i]);
		//Calculates stabilization term which is a part of numerical
		//flux vector i.e. |A|(Q{R}-Q{L})
		// Roe type averaging procedure first
		double ro_l = UL.ro;
		double ro_r = UR.ro;
		Vector velocity_l, velocity_r;
		velocity_l.x = UL.rou/ro_l;
		velocity_l.y = UL.rov/ro_l;
		velocity_l.z = UL.row/ro_l;
		velocity_r.x = UR.rou/ro_r;
		velocity_r.y = UR.rov/ro_r;
		velocity_r.z = UR.row/ro_r;
		double e_l = UL.roE/ro_l;
		double e_r = UR.roE/ro_r;
		double k;
		k = 0.5*(velocity_l*velocity_l);		//kinetik energy
		double h_l = (e_l-k)*gamma + k;	//enthalpy
		k = 0.5*(velocity_r*velocity_r);		//kinetik energy
		double h_r = (e_r-k)*gamma + k;	//enthalpy
		double ro = sqrt(ro_l*ro_r);          // (Roe averaged) density		
		double ql  = sqrt(ro_l)/(sqrt(ro_l)+sqrt(ro_r));
		double qr  = sqrt(ro_r)/(sqrt(ro_l)+sqrt(ro_r));
		velocity = ql*velocity_l + qr*velocity_r;	// (Roe averaged) velocity	
		double h  = ql*h_l + qr*h_r;  // (Roe averaged) total enthalpy
		//Proceed to solution
		double phi2 = 0.5*(gamma - 1)*(velocity*velocity);
		double dn = f.FaceNormal.mod();
		double c = sqrt((gamma - 1)*h - phi2);	//acoustic velocity
		//Debug	
		double uw = velocity * f.FaceNormal;
		double eig_max = fabs(uw)+c*dn;	
		double AA1 = Harten(uw, eps*eig_max);       // AA1, AA3, AA1 -
		double AA3 = Harten(uw+c*dn, eps*eig_max);  // eigenvalues of a flux vector
		double AA4 = Harten(uw-c*dn, eps*eig_max);  // Jacobian matrix
		double Eig1= AA1;
		double Eiga=(AA3 - AA4)*0.5/(dn*c);
		double Eigb=(AA3 + AA4)*0.5 - AA1;
		//parametrs vectors Qa and Qb (i guess)
		std::vector<double> Qa(5, 0);
		std::vector<double> Qb(5, 0);
		Qa[0]=0;
		Qa[1]=f.FaceNormal.x;
		Qa[2]=f.FaceNormal.y;
		Qa[3]=f.FaceNormal.z;
		Qa[4]=uw;    
		Qb[0]=1;
		Qb[1]=velocity.x;
		Qb[2]=velocity.y;
		Qb[3]=velocity.z;
		Qb[4]=h;
		//Calculate solution
		//Some quotients
		double R1 =phi2*ro_r-(gamma-1)*ro_r*(velocity*velocity_r-e_r);	//PR =R1
		double D1 =ro_r*f.FaceNormal*(velocity_r - velocity);				//DR =D1

		double R2 =phi2*ro_l-(gamma-1)*ro_l*(velocity*velocity_l-e_l);	//PR =R1
		double D2 =ro_l*f.FaceNormal*(velocity_l - velocity);		

		double C2 = 1.0/(c*c);
		double DN2= 1.0/(dn*dn);

		std::vector<double> ul(5,0);
		ul[0] = UL.ro;
		ul[1] = UL.rou;
		ul[2] = UL.rov;
		ul[3] = UL.row;
		ul[4] = UL.roE;
		std::vector<double> ur(5,0);
		ur[0] = UR.ro;
		ur[1] = UR.rou;
		ur[2] = UR.rov;
		ur[3] = UR.row;
		ur[4] = UR.roE;
		for(int i=0; i<5; i++){
				  res[i]+=(Eig1*ul[i] + Eiga*(R2*Qa[i]     +D2*Qb[i])
									+ Eigb*(R2*Qb[i]*C2  +D2*Qa[i]*DN2)
						-Eig1*ur[i] - Eiga*(R1*Qa[i]     +D1*Qb[i])
									- Eigb*(R1*Qb[i]*C2  +D1*Qa[i]*DN2));
				  res[i]*= 0.5;
		};	

		//Preconditioning
	/*	const double Cp = 1006.43;
		const double Cv = Cp / gamma;
		const double R = Cp - Cv;
		double T = h / Cp;
		double P = ro * R * T;
		double dro_dT =  - ro / T;
		double Ur = uw;		
		double teta = (1.0/(Ur*Ur) - dro_dT/(ro*Cp));*/		

		MaxEigenvalue = eig_max;	
		return res;
	};

	//fabs() and Harten's entropy correction procedure
	double Harten(double z, double eps) 
	{
		z = fabs(z);
		if (z<eps) z = ((z*z)/eps + eps)*0.5;
		return z;
	};

	//Solve riemann problem	
	RiemannProblemSolutionResult Solve(int nmatL, const GasModel::ConservativeVariables& UL, int nmatR, const GasModel::ConservativeVariables& UR, const Face& f, double ALEindicator) {
		RiemannProblemSolutionResult result;
		result.FluxesLeft.resize(_gasModels[nmatL]->nConservativeVariables, 0);
		result.FluxesRight.resize(_gasModels[nmatR]->nConservativeVariables, 0);

		result.Fluxes = ComputeFlux(UL, UR, f);
		result.MaxEigenvalue = MaxEigenvalue;
		result.Velocity = velocity;
		result.Pressure = 0;

		return result;
	};
};

#endif
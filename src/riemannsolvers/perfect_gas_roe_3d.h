#ifndef TURBO_RIEMANNSOLVERS_PERFECT_GAS_ROE_3D
#define TURBO_RIEMANNSOLVERS_PERFECT_GAS_ROE_3D

#include "datatypes.h"
#include "basetypes.h"
#include "grid.h"

class Roe3DSolverPerfectGas {
	//Required info
	double gamma;
	double eps;	
public:	
	Roe3DSolverPerfectGas(){
		eps = 0.05;
	};

	void SetGamma(double _g) {
		gamma = _g;
	};

	void SetHartenEps(double _eps) {
		eps = _eps;
	};

	//Numerical flux
	std::vector<double> F(ConservativeVariables U, Vector n)
	{		
		std::vector<double> res(5,0);		
		double ro = U.ro;
		double vx = U.rou/ro;
		double vy = U.rov/ro;
		double vz = U.row/ro;
		double roE = U.roE;	//ro*e
		double p = (gamma - 1.0)*(roE - ro*(vx*vx + vy*vy + vz*vz)/2.0);
		double vn = vx*n.x + vy*n.y + vz*n.z;

		res[0] = ro*vn;
		res[1] = ro*vn*vx + n.x*p;
		res[2] = ro*vn*vy + n.y*p;
		res[3] = ro*vn*vz + n.z*p;
		res[4] = vn*(roE + p);
		return res;
	};

	//Solve riemann problem
	std::vector<double> ComputeFlux(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f) {
		std::vector<double> res(5,0);
		
		//Calculate symmetric flux part		
		res += 1.0*(F(UL, (1.0)*f.FaceNormal) + F(UR, (1.0)*f.FaceNormal)); //is 1.0 or 0.5 TODO		
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
		Vector velocity = ql*velocity_l + qr*velocity_r;	// (Roe averaged) velocity	
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

		MaxEigenvalue = eig_max;	
		return res;
	};
	std::vector<double> ComputeFlux2D(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f) {
		
		//face flux
		std::vector<double> res(5,0);
		//Tangent vector
		Vector n = f.FaceNormal;
		Vector m(0, 0, 0);
		m.x = -n.y;
		m.y = n.x;
		
		//Primitive and other variables
		//Left state
		double roL = UL.ro;
		double uL = UL.rou/UL.ro;
		double vL = UL.rov/UL.ro;
		double unL = uL*n.x + vL*n.y;
		double umL = uL*m.x + vL*m.y;
		double EL = UL.roE/UL.ro;
		double pL = (gamma - 1.0)*(roL*EL - 0.5*roL*(uL*uL + vL*vL));
		double aL = sqrt(gamma*pL/roL);
		double HL = EL + pL/roL;

		//Right state
		double roR = UR.ro;
		double uR = UR.rou/UR.ro;
		double vR = UR.rov/UR.ro;
		double unR = uR*n.x + vR*n.y;
		double umR = uR*m.x + vR*m.y;
		double ER = UR.roE/UR.ro;
		double pR = (gamma - 1.0)*(roR*ER - 0.5*roR*(uR*uR + vR*vR));
		double aR = sqrt(gamma*pR/roR);
		double HR = ER + pR/roR;

		//compute the Roe Averages
		double RT = sqrt(roR/roL);
		double ro = RT*roL;
		double u = (uL + RT*uR)/(1.0 + RT);
		double v = (vL + RT*vR)/(1.0 + RT);
		double H = (HL + RT*HR)/(1.0 + RT);
		double a = sqrt((gamma - 1.0)*(H - 0.5*(u*u + v*v)));
		double un = u*n.x + v*n.y;
		double um = u*m.x + v*m.y;

		//Wave Strengths
		double  dro = roR - roL; 
		double dp = pR - pL;
		double dun = unR - unL;
		double dum = umR - umL;

		//Rieman's variables
		std::vector<double> LdU(4);
		LdU[0] = (dp - ro*a*dun )/(2.0*a*a);
		LdU[1] = ro*dum;
		LdU[2] =  dro - dp/(a*a);
		LdU[3] = (dp + ro*a*dun )/(2.0*a*a);

		//Wave Speed
		std::vector<double> ws(4);
		ws[0] = abs(un-a);
		ws[1] = abs(un);
		ws[2] = abs(un);
		ws[3] = abs(un+a);

		//Harten's Entropy Fix JCP(1983), 49, pp357-393
		if(ws[0]<eps) ws[0] = 0.5*(ws[0]*ws[0]/eps + eps);
		if(ws[3]<eps) ws[3] = 0.5*(ws[3]*ws[3]/eps + eps);

		//Right Eigenvectors
		std::vector<std::vector<double>> Rv(4);
		for(int i=0; i<4; i++) Rv[i].resize(4);
		Rv[0][0] = 1.0;    
		Rv[1][0] = u - a*n.x;
		Rv[2][0] = v - a*n.y;
		Rv[3][0] = H - un*a;

		Rv[0][1] = 0;
		Rv[1][1] = m.x;
		Rv[2][1] = m.y;
		Rv[3][1] = um;

		Rv[0][2] = 1.0;
		Rv[1][2] = u;
		Rv[2][2] = v;
		Rv[3][2] = 0.5*(u*u+v*v);

		Rv[0][3] = 1.0;
		Rv[1][3] = u + a*n.x;
		Rv[2][3] = v + a*n.y;
		Rv[3][3] = H + un*a;

		//Dissipation Term
		std::vector<double> diss(4, 0);
		for(int i=0; i<4; i++)
		   for(int j=0; j<4; j++)
			diss[i] += ws[j]*LdU[j]*Rv[i][j];

		//Compute right and left flux
		std::vector<double> fL(4);
		fL[0] = roL*unL;
		fL[1] = roL*unL * uL + pL*n.x;
		fL[2] = roL*unL * vL + pL*n.y;
		fL[3] = roL*unL * HL;

		std::vector<double> fR(4);
		fR[0] = roR*unR;
		fR[1] = roR*unR * uR + pR*n.x;
		fR[2] = roR*unR * vR + pR*n.y;
		fR[3] = roR*unR * HR;

		//compute fluxes
		std::vector<double> flux = 0.5*(fL + fR - diss);
		res[0] = flux[0];
		res[1] = flux[1];
		res[2] = flux[2];
		res[4] = flux[3];
		MaxEigenvalue =(abs(un) + a);  //Normal max wave speed

		return res;
	};
	std::vector<double> ComputeFluxDirect(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f, const Vector& c1, const Vector& c2) {
		std::vector<double> res(5,0);
		Vector F_center = f.FaceCenter;
		double r_l = (c1 - F_center).mod();
		double r_r = (c2 - F_center).mod();
		double r = r_l + r_r;
		r_l /= r;
		r_r /= r;

		ConservativeVariables U_av;
		U_av.ro = 0.5*(UL.ro + UR.ro);
		U_av.roE = 0.5*(UL.roE + UR.roE);
		U_av.rou = 0.5*(UL.rou + UR.rou);
		U_av.rov = 0.5*(UL.rov + UR.rov);
		U_av.row = 0.5*(UL.row + UR.row);
		res += F(U_av, f.FaceNormal);

		double F_ro = r_r*UL.ro + r_l*UR.ro;
		double F_vx =r_r*(UL.rou/UL.ro) + r_l*(UR.rou/UL.ro);
		double F_vy =r_r*(UL.rov/UL.ro) + r_l*(UR.rov/UL.ro);
		double F_vz =r_r*(UL.row/UL.ro) + r_l*(UR.row/UL.ro);
		double e_in_l = UL.roE - 0.5*(UL.rou*UL.rou + UL.rov*UL.rov + UL.row*UL.row)/UL.ro;
		e_in_l /= UL.ro;
		double e_in_r = UR.roE - 0.5*(UR.rou*UR.rou + UR.rov*UR.rov + UR.row*UR.row)/UR.ro;
		e_in_r /= UR.ro;
		double F_e = r_r*e_in_l + r_l*e_in_r;		//specific internal energy
		double F_E = F_e + 0.5*(F_vx*F_vx + F_vy*F_vy + F_vz*F_vz);
		double F_p = (gamma-1.0)*F_ro*F_e;

		Vector n = f.FaceNormal;

		res[0] = n.x * (F_ro*F_vx) + n.y*(F_ro*F_vy) + n.z*(F_ro*F_vz);	//WID vn is not used
		res[1] = n.x * (F_ro*F_vx*F_vx + F_p) + n.y*(F_ro*F_vx*F_vy) + n.z*(F_ro*F_vx*F_vz);
		res[2] = n.x * (F_ro*F_vy*F_vx) + n.y*(F_ro*F_vy*F_vy + F_p) + n.z*(F_ro*F_vy*F_vz);
		res[3] = n.x * (F_ro*F_vz*F_vx) + n.y*(F_ro*F_vz*F_vy) + n.z*(F_ro*F_vz*F_vz + F_p);
		res[4] = (n.x * F_vx + n.y * F_vy + n.z * F_vz)*(F_ro*F_e + F_p);	//WID	
		
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
		Vector velocity = ql*velocity_l + qr*velocity_r;	// (Roe averaged) velocity	
		double h  = ql*h_l + qr*h_r;  // (Roe averaged) total enthalpy
		//Proceed to solution
		double phi2 = 0.5*(gamma - 1)*(velocity*velocity);
		double dn = f.FaceNormal.mod();
		double c = sqrt((gamma - 1)*h - phi2);	//acoustic velocity
		//Debug	
		double uw = velocity * f.FaceNormal;
		double eig_max = fabs(uw)+c*dn;	

		MaxEigenvalue = eig_max;	
		return res;
	};
	std::vector<double> ComputeFluxSemiDirect(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f) {
		std::vector<double> res(5,0);
		
		//Calculate symmetric flux part	
		ConservativeVariables U_av = UL;
		U_av.ro = 0.5*(UL.ro + UR.ro);
		U_av.roE = 0.5*(UL.roE + UR.roE);
		U_av.rou = 0.5*(UL.rou + UR.rou);
		U_av.rov = 0.5*(UL.rov + UR.rov);
		U_av.row = 0.5*(UL.row + UR.row);
		res += F(U_av, f.FaceNormal);

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
		Vector velocity = ql*velocity_l + qr*velocity_r;	// (Roe averaged) velocity	
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
		};	

		MaxEigenvalue = eig_max;	
		return res;
	};
	std::vector<double> ComputeFluxTrivial(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f) {
		std::vector<double> res(5,0);
		
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
		Vector velocity = ql*velocity_l + qr*velocity_r;	// (Roe averaged) velocity	
		double h  = ql*h_l + qr*h_r;  // (Roe averaged) total enthalpy
		//Proceed to solution
		double phi2 = 0.5*(gamma - 1)*(velocity*velocity);
		double dn = f.FaceNormal.mod();
		double c = sqrt((gamma - 1)*h - phi2);	//acoustic velocity
		//Debug	
		double uw = velocity * f.FaceNormal;
		double eig_max = fabs(uw)+c*dn;	
		
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


	//Public properties
	double MaxEigenvalue;
};

#endif
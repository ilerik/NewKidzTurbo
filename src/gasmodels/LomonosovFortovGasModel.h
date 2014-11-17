#ifndef TURBO_GasModels_LomonosovFortovGasModel
#define TURBO_GasModels_LomonosovFortovGasModel

#include "GasModel.h"
#include "utilityfunctions.h"

class LomonosovFortovGasModel : public GasModel {
private:	
	//    List of materials:
    //    1 - Be   2 - Mg     3 - Al     4 - e-Fe     5 - Cr
    //    6 - Co   7 - Ni     8 - Cu     9 - Pb      10 - Mo
    //    11 - Cd  12 - Ta    13 - W     14 - UO2     15 - Ti
	//    16 - Au  17 - LiF   18 - corund 19 - quartz  20 - lix
	int _nmat; // Index of material
public:
	//Constructor specifies material to use
	LomonosovFortovGasModel(int nmat) {
		//Model information
		GasModelName = "LomonosovFortovGasModel";
		GasModelType = ModelType_t::ModelTypeUserDefined;

		//Material index
		_nmat = nmat;

		//Fill material dependent constants
		_prepareData();
	};

	//Load configuration
	void loadConfiguration(Configuration _conf) {
		//Not implemented yet
		nConservativeVariables = 5;
	};

	//Get pressure
	double GetPressure(GasModel::ConservativeVariables U) {
		double ro = U.ro;
		double u = U.rou / ro;
		double v = U.rov / ro;
		double w = U.row / ro;
		double E = U.roE / ro;
		double e = E - (u*u + v*v + w*w) / 2.0;
		double P;
		double c;
		double Gr;
		bool NF;
		e /= 1e6; //J/kg -> kJ/g
		ro /= 1000; //kg/cc -> g/cc
		EOSE5(ro, e, P, c, Gr, NF);
		P *= 1e9;  //GPa -> Pa
		c *= 1000; //km\s -> m\s
		return P;
	};

	//Get pressure and sound speed
	void GetPressureAndSoundSpeed(GasModel::ConservativeVariables U, double& pressure, double& soundspeed, double& gruneisen) {
		double ro = U.ro;
		double u = U.rou / ro;
		double v = U.rov / ro;
		double w = U.row / ro;
		double E = U.roE / ro;
		double e = E - (u*u + v*v + w*w) / 2.0;		
		bool NF; //If solution is non-physical
		e /= 1e6; //J/kg -> kJ/g
		ro /= 1000; //kg/cc -> g/cc
		EOSE5(ro, e, pressure, soundspeed, gruneisen, NF);
		pressure *= 1e9;  //GPa -> Pa
		soundspeed *= 1000; //km\s -> m\s	
	};

	//Get temperature

	//Get internal energy	
	struct EOSIdentityParams {
		double ro;
		double p;
		LomonosovFortovGasModel* gasModel;
	};

	static void EOSIndentity(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *object) {
		EOSIdentityParams* params = (EOSIdentityParams*)object;		
		double e = x[0];
		double P;
		double c;
		double Gr;
		bool NF;
		params->gasModel->EOSE5(params->ro, e, P, c, Gr, NF);
		fi[0] = P - params->p;
		if (NF) fi[0] += 10.0; //Non physical penalty
	};

	//Given pressure and density numerically find specific internal energy
	double FindInternalEnergy(double ro, double p) {
		//Convert quatities to new system of units
		ro /= 1e3; //to g/cc
		p /= 1e9; //to GPa

		//Init vectors
		alglib::real_1d_array x;
		alglib::real_1d_array bl;
		alglib::real_1d_array bu;
		x.setlength(1);		
		bl.setlength(1);
		bu.setlength(1);

		//Set boundary constraints (non negative specific energy)
		bl[0] = 0.0;
		bu[0] = 100.0;

		//Initial guess (possible othe choises, the simplest for now)
		x[0] = 1.0;
		
		EOSIdentityParams params;
		params.ro = ro;
		params.p = p;
		params.gasModel = this;

		//Iterative scheme parameters
		double epsg = 1e-14;
		double epsf = 0;
		double epsx = 0;
		double diffstep = 1e-10;
		alglib::ae_int_t maxits = 0;
		alglib::minlmstate state;
		alglib::minlmreport rep;
		
		alglib::minlmcreatev(1, x, diffstep, state);
		alglib::minlmsetbc(state, bl, bu);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);				
		alglib::minlmoptimize(state, EOSIndentity, NULL, (void*)this);		
		alglib::minlmresults(state, x, rep);

		double e = x[0];

		//convert enery units to SI
		e *= 1e6; //kJ/g to J/kg
		return e;
	};

	//    Equation of state
	//    ---------- P=P(E,X), X=LN(V0/V) -----------
	//    Input parameters:
	//    NMAT - material number, real
	//    RO [g/cc] density (initial density R0=1./V0, see array V0
	//    EE [kJ/g] internal energy
	//    Output paramters:
	//    PP [GPa] pressure
	//    C  [km/s] isentropic sound speed
	//    NF [-] physical state, 0=unphysical, 1=ok
	//	  g [-] Gruniesen coefficient ?
	//    List of materials:
	void EOSE5(double ro, double e, double& P, double& C, double& g, bool& nonPhysical) {
		int nmet = _nmat;
		double x = log(ro * V0[nmet]);
		e = e + E0[nmet];
		double vr = exp(x);
		double v = DX[nmet] * vr;
		double v0 = V0[nmet]; //specific material volume at T = 0
		double px = 0;
		double gx = 0;
		double c2x = 0;
		double et = 0;
		double g1x = 0;

		if ((1.0 - v) > 0) {
			//Cold lattice contribution for liquid
			double cm = 0;
			if (x > -100.0/GM[nmet]) cm = CMN[nmet]*pow(v, GM[nmet]);
			double cn = 0;
			if (x > -100.0/GN[nmet]) cn = CMN[nmet]*pow(v, GN[nmet]);
			px = cm - cn;
			double ex = cm / GM[nmet] - cn / GN[nmet] + ES[nmet];
			c2x = px + cm * GM[nmet] - cn * GN[nmet];
			px = px * v / (v0 * DX[nmet]);
			double sq = GC[nmet] * pow(v, QS[nmet]);
			double sr = sq*pow(v/SM[nmet], RS[nmet]);
			double g = sq/QS[nmet];
			double gr = sr/(QS[nmet] + RS[nmet]);
			gx = GGI + g - gr;
			et = e - ex;
			double g1x = sq - sr;
		} else {
			//Cold lattice contribution for solid
			double v13 = pow(v,R13);
			double v23 = v13*v13;
			double v33 = v;
			double v43 = v23*v23;
			double v53 = v23*v33;
			double v1 = A1[nmet]*v13;
			double v2 = A2[nmet]*v23;
			double v3 = A3[nmet]*v33;
			double v4 = A4[nmet]*v43;
			double v5 = A5[nmet]*v53;
			double ex = 3.0 * (v1-A1[nmet] 
			+ 0.5*(v2-A2[nmet]) 
			+ (v3 - A3[nmet])/3.0
			+ 0.25*(v4-A4[nmet]) 
			+ 0.2*(v5-A5[nmet])
			);
			double r0 = v1 + v2 + v3 + v4 + v5;
			double r1 = 4./3.*v1+5./3.*v2+2.*v3+7./3.*v4+8./3.*v5;
			double r2 = 4./9.*v1+10./9.*v2+2.*v3+28./9.*v4+40./9.*v5;
			double r3 = -8./27.*v1-10./27.*v2+28./27.*v4+80./27.*v5;
			ex = ex*V0[nmet]*DX[nmet];
			px = r0*v;
			double c2x = r1*V0[nmet]*DX[nmet];
			r3 = r3 + TA[nmet] * ((r1*(TA[nmet] + 1.0) - r2) * 3.
			- r0 * (TA[nmet] + 1.) * (TA[nmet] + 2.0));
			r2 = r2 + TA[nmet] * (r0*(TA[nmet] + 1.0) - 2.*r1);
			r1 = r1 - r0*TA[nmet];
			double gr = r2/r1;
			gx = 0.5 * (TA[nmet] + gr + GGI);
			et = e - ex;
			g1x = 0.5 * (gr*(1.-gr)+r3/r1);
		};

		// Full pressure, sound speed and g calculations ---
		double eav = EA[nmet] * pow(v, GGI);
		v = V0[nmet]/vr;
		g = GI[nmet] + (gx - GI[nmet]) * GF[nmet]/(1.0+et/eav);
		double pt = g*et/v;
		P = pt + px;
		e = e - E0[nmet];
		double c2s = 0;

		// Check pressure position  ---
		if (et < 0) {
			//Non physical state
			nonPhysical = true;
			P=1.E-8; //dummy value for pressure
			C=1.E-4; //dummy value for sound speed
			return;
		} else {
			//Correct physical state
			nonPhysical = false;
			if ((P < 1.E-4) && (vr < 0.7)) { //TO DO why hardcode?
				c2s = CR[nmet]*CR[nmet]*GI[nmet]*(et+pt*v)/C1R[nmet];
				c2s = sign(c2s) * sqrt(abs(c2s));
			} else {
				double e1x = GGI*EA[nmet];
				double ga = g - GF[nmet] * (gx - GI[nmet])/(1.+et/eav)/(1.+et/eav)*et/eav;
				c2s = c2x 
					+ ga*P*v 
					+ g*(et-px*v) 
					+ et*GF[nmet]/(1.+et/eav)*(g1x+(gx-GI[nmet])/(et+eav)*(px*v+et*e1x/eav));            
				c2x = sign(c2x) * sqrt(abs(c2x)); 
				c2s = sign(c2s) * sqrt(abs(c2s));
			};
		};
		C = c2s;
		return;
	};
private:
	//Model constants for all materials
	const static double R13;
	const static double GGI;
		
	//vectors of usefull data
	std::vector<double> V0;
	std::vector<double> DX;
	std::vector<double> CMN;
	std::vector<double> GM;
	std::vector<double> GN;
	std::vector<double> ES;
	std::vector<double> E0;
	std::vector<double> A1;
	std::vector<double> A2;
	std::vector<double> A3;
	std::vector<double> A4;
	std::vector<double> A5;
	std::vector<double> TA;
	std::vector<double> GC;
	std::vector<double> GF;
	std::vector<double> GI;
	std::vector<double> QS;
	std::vector<double> RS;
	std::vector<double> SM;
	std::vector<double> EA;
	std::vector<double> CR;
	std::vector<double> C1R;
	
	// fill corresponding vectors of coeffiecients
	void _prepareData()
	{
		//0 - stainless steel
		V0.push_back(0.12700000E+00);
		DX.push_back(0.98664743E+00);
		CMN.push_back(0.95555325E+01);
		GM.push_back(0.30000000E+01);
		GN.push_back(0.90187567E+00);
		ES.push_back(0.74099998E+01);
		E0.push_back(0.13392337E+00);
		A1.push_back(-.10226000E+04);
		A2.push_back(0.37132200E+04);
		A3.push_back(-.50760498E+04);
		A4.push_back(0.26228401E+04);
		A5.push_back(-.23741000E+03);
		TA.push_back(0.62708032E+00);
		GC.push_back(0.21484993E+02);
		GF.push_back(0.10031745E+01);
		GI.push_back(0.50000000E+00);
		QS.push_back(0.22000000E+01);
		RS.push_back(0.50000000E+00);
		SM.push_back(0.88413364E+00);
		EA.push_back(0.42000000E+02);
		CR.push_back(0.34523895E+01);
		C1R.push_back(0.32979987E+01);

		//1 - Pb
		V0.push_back(0.88200003E-01);
		DX.push_back(0.97725165E+00);
		CMN.push_back(-0.14261333E+02);
		GM.push_back(0.20000000E+01);
		GN.push_back(0.23021955E+01);
		ES.push_back(0.93599999E+00);
		E0.push_back(0.36263671E-01);
		A1.push_back(0.16000900E+03);
		A2.push_back(-.28083600E+03);
		A3.push_back(-.83797997E+02);
		A4.push_back(0.22006700E+03);
		A5.push_back(-.15442000E+02);
		TA.push_back(-.33601224E+00);
		GC.push_back(0.39365627E+02);
		GF.push_back(0.10038391E+01);
		GI.push_back(0.46000001E+00);
		QS.push_back(0.20000000E+01);
		RS.push_back(0.50000000E+00);
		SM.push_back(0.78492051E+00);
		EA.push_back(0.93000002E+01);
		CR.push_back(0.16488581E+01);
		C1R.push_back(0.81879050E+00);
		/*V0.push_back(0.54200000E+00);V0.push_back(0.57700002E+00);V0.push_back(0.36899999E+00);V0.push_back(0.12000000E+00);
		V0.push_back(0.13900000E+00);V0.push_back(0.11300000E+00);V0.push_back(0.11200000E+00);V0.push_back(0.11200000E+00);
		V0.push_back(0.88200003E-01);V0.push_back(0.97800002E-01);V0.push_back(0.11400000E+00);V0.push_back(0.59700001E-01);
		V0.push_back(0.51899999E-01);V0.push_back(0.91200002E-01);V0.push_back(0.22200000E+00);V0.push_back(0.51800001E-01);
		V0.push_back(0.37900001E+00);V0.push_back(0.25099999E+00);V0.push_back(0.23400000E+00);V0.push_back(0.11200000E+01);
    
		DX.push_back(0.98145735E+00);DX.push_back(0.97675198E+00);DX.push_back(0.97923738E+00);DX.push_back(0.98654467E+00);
		DX.push_back(0.98867989E+00);DX.push_back(0.98776722E+00);DX.push_back(0.98788631E+00);DX.push_back(0.98512793E+00);
		DX.push_back(0.97725165E+00);DX.push_back(0.99516839E+00);DX.push_back(0.97772300E+00);DX.push_back(0.99392861E+00);
		DX.push_back(0.99559391E+00);DX.push_back(0.97596389E+00);DX.push_back(0.99185765E+00);DX.push_back(0.98891979E+00);
		DX.push_back(0.98592889E+00);DX.push_back(0.99903792E+00);DX.push_back(0.99866945E+00);DX.push_back(0.95384234E+00);
    
		CMN.push_back(0.42994564E+02);CMN.push_back(0.11484624E+02);CMN.push_back(0.13236353E+02);CMN.push_back(0.93908024E+01);
		CMN.push_back(0.17328915E+02);CMN.push_back(0.13348364E+02);CMN.push_back(0.40190891E+02);CMN.push_back(0.32090862E+02);
		CMN.push_back(-0.14261333E+02);CMN.push_back(0.15142752E+02);CMN.push_back(0.16921934E+02);CMN.push_back(0.55232458E+01);
		CMN.push_back(0.88132391E+01);CMN.push_back(0.94326097E+00);CMN.push_back(0.10463772E+02);CMN.push_back(0.80316305E+01);
		CMN.push_back(0.95852518E+01);CMN.push_back(0.27481421E+02);CMN.push_back(-0.65462183E+03);CMN.push_back(0.10416588E+02);
    
		GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);
		GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.20000000E+01);GM.push_back(0.20000000E+01);
		GM.push_back(0.20000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);
		GM.push_back(0.30000000E+01);GM.push_back(0.90000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);
		GM.push_back(0.30000000E+01);GM.push_back(0.30000000E+01);GM.push_back(0.10000000E+01);GM.push_back(0.30000000E+01);
    
		GN.push_back(0.14658144E+01);GN.push_back(0.11745121E+01);GN.push_back(0.79678905E+00);GN.push_back(0.88841671E+00);
		GN.push_back(0.14138776E+01);GN.push_back(0.12774221E+01);GN.push_back(0.14686730E+01);GN.push_back(0.15083531E+01);
		GN.push_back(0.23021955E+01);GN.push_back(0.12716897E+01);GN.push_back(0.25520797E+01);GN.push_back(0.89647335E+00);
		GN.push_back(0.11707672E+01);GN.push_back(0.39387769E+00);GN.push_back(0.79043901E+00);GN.push_back(0.18141516E+01);
		GN.push_back(0.11521416E+00);GN.push_back(0.30822426E+00);GN.push_back(0.11392220E+01);GN.push_back(0.38477069E+00);
    
		ES.push_back(0.15000000E+02);ES.push_back(0.59499998E+01);ES.push_back(0.12200000E+02);ES.push_back(0.74400001E+01);
		ES.push_back(0.64800000E+01);ES.push_back(0.60000000E+01);ES.push_back(0.72700000E+01);ES.push_back(0.52300000E+01);
		ES.push_back(0.93599999E+00);ES.push_back(0.68600001E+01);ES.push_back(0.99000001E+00);ES.push_back(0.43200002E+01);
		ES.push_back(0.45900002E+01);ES.push_back(0.22900000E+01);ES.push_back(0.97500000E+01);ES.push_back(0.17500000E+01);
		ES.push_back(0.80000000E+02);ES.push_back(0.80000000E+02);ES.push_back(0.80000000E+02);ES.push_back(0.23600000E+02);
    
		E0.push_back(0.82284731E+00);E0.push_back(0.30740389E+00);E0.push_back(0.27874166E+00);E0.push_back(0.13305190E+00);
		E0.push_back(0.14209883E+00);E0.push_back(0.12565036E+00);E0.push_back(0.12638120E+00);E0.push_back(0.11601043E+00);
		E0.push_back(0.36263671E-01);E0.push_back(0.76777481E-01);E0.push_back(0.66890948E-01);E0.push_back(0.40647842E-01);
		E0.push_back(0.40003933E-01);E0.push_back(0.90129100E-01);E0.push_back(0.15371029E+00);E0.push_back(0.37791826E-01);
		E0.push_back(0.29866642E+00);E0.push_back(0.71235187E-01);E0.push_back(0.12139247E+00);E0.push_back(0.84046292E+00);
    
		A1.push_back(-0.13212900E+04);A1.push_back(0.35834999E+02);A1.push_back(0.39950000E+03);A1.push_back(-0.10416100E+04);
		A1.push_back(0.69557001E+03);A1.push_back(0.70087701E+03);A1.push_back(0.18928400E+04);A1.push_back(-0.65371002E+03);
		A1.push_back(0.16000900E+03);A1.push_back(0.45810999E+03);A1.push_back(0.84567999E+03);A1.push_back(-0.25452299E+02);
		A1.push_back(-0.40144601E+03);A1.push_back(0.74053003E+03);A1.push_back(-0.29326501E+03);A1.push_back(0.18746400E+04);
		A1.push_back(0.45814001E+03);A1.push_back(-0.81706201E+03);A1.push_back(-0.23182000E+04);A1.push_back(0.13791600E+03);
    
		A2.push_back(0.31578201E+04);A2.push_back(-0.19723599E+03);A2.push_back(-0.11932000E+04);A2.push_back(0.42018701E+04);
		A2.push_back(-0.20732300E+04);A2.push_back(-0.20388700E+04);A2.push_back(-0.45849902E+04);A2.push_back(0.28425300E+04);
		A2.push_back(-0.28083600E+03);A2.push_back(-0.19048800E+04);A2.push_back(-0.19349200E+04);A2.push_back(-0.31007999E+03);
		A2.push_back(0.34389801E+03);A2.push_back(-0.16957900E+04);A2.push_back(0.34880301E+03);A2.push_back(-0.42837700E+04);
		A2.push_back(-0.12964900E+04);A2.push_back(0.67171399E+03);A2.push_back(0.37505901E+04);A2.push_back(-0.41737799E+03);
    
		A3.push_back(-0.27733201E+04);A3.push_back(0.17582001E+03);A3.push_back(0.95741998E+03);A3.push_back(-0.60544502E+04);
		A3.push_back(0.14612500E+04);A3.push_back(0.13551400E+04);A3.push_back(0.29190200E+04);A3.push_back(-0.44766201E+04);
		A3.push_back(-0.83797997E+02);A3.push_back(0.16349200E+04);A3.push_back(0.11284700E+04);A3.push_back(0.97382004E+02);
		A3.push_back(-0.44822000E+03);A3.push_back(0.89040997E+03);A3.push_back(-0.13732100E+03);A3.push_back(0.23809600E+04);
		A3.push_back(0.10074200E+04);A3.push_back(0.23222800E+03);A3.push_back(-0.17300400E+04);A3.push_back(0.35004999E+03);
    
		A4.push_back(0.98634003E+03);A4.push_back(-0.14870000E+02);A4.push_back(-0.17536000E+03);A4.push_back(0.31672300E+04);
		A4.push_back(-0.85099998E+02);A4.push_back(-0.15187000E+02);A4.push_back(-0.23344000E+03);A4.push_back(0.26115000E+04);
		A4.push_back(0.22006700E+03);A4.push_back(-0.19434000E+03);A4.push_back(-0.38910000E+02);A4.push_back(0.24988499E+03);
		A4.push_back(0.53452899E+03);A4.push_back(0.70820000E+02);A4.push_back(0.86291000E+02);A4.push_back(0.33049999E+02);
		A4.push_back(-0.17992999E+03);A4.push_back(-0.96350998E+02);A4.push_back(0.31112000E+03);A4.push_back(-0.76129997E+02);
    
		A5.push_back(-0.49549999E+02);A5.push_back(0.45100001E+00);A5.push_back(0.11640000E+02);A5.push_back(-0.27304001E+03);
		A5.push_back(0.15100000E+01);A5.push_back(-0.19600000E+01);A5.push_back(0.65700002E+01);A5.push_back(-0.32370001E+03);
		A5.push_back(-0.15442000E+02);A5.push_back(0.61900001E+01);A5.push_back(-0.31999999E+00);A5.push_back(-0.11734700E+02);
		A5.push_back(-0.28761000E+02);A5.push_back(-0.59699998E+01);A5.push_back(-0.45079999E+01);A5.push_back(-0.48800001E+01);
		A5.push_back(0.10860000E+02);A5.push_back(0.94709997E+01);A5.push_back(-0.13470000E+02);A5.push_back(0.55419998E+01);
    
		TA.push_back(-0.41694292E+00);TA.push_back(0.40492105E+00);TA.push_back(-0.14242107E+00);TA.push_back(0.11158217E+01);
		TA.push_back(-0.25900945E+00);TA.push_back(-0.29996026E+00);TA.push_back(0.12563297E+01);TA.push_back(0.88399196E+00);
		TA.push_back(-0.33601224E+00);TA.push_back(0.29947156E+00);TA.push_back(0.12490201E+01);TA.push_back(0.10827350E+00);
		TA.push_back(-0.15773018E+00);TA.push_back(0.13764737E+01);TA.push_back(0.44122058E+00);TA.push_back(0.17943431E+00);
		TA.push_back(0.17006929E+01);TA.push_back(0.66973019E+00);TA.push_back(0.17068748E+00);TA.push_back(0.11692852E+01);
    
		GC.push_back(0.11733573E+02);GC.push_back(0.13373021E+02);GC.push_back(0.27289089E+02);GC.push_back(0.21815765E+02);
		GC.push_back(0.28064201E+02);GC.push_back(0.11755576E+02);GC.push_back(0.21993205E+03);GC.push_back(0.20048616E+02);
		GC.push_back(0.39365627E+02);GC.push_back(0.29287760E+02);GC.push_back(0.37389256E+02);GC.push_back(0.31843267E+02);
		GC.push_back(0.18399178E+02);GC.push_back(0.98581100E+00);GC.push_back(0.76616817E+01);GC.push_back(0.45784122E+02);
		GC.push_back(0.98853817E+01);GC.push_back(0.47292509E+01);GC.push_back(0.61693435E+01);GC.push_back(0.11861769E+02);
    
		GF.push_back(0.10456556E+01);GF.push_back(0.10218971E+01);GF.push_back(0.10197377E+01);GF.push_back(0.10031537E+01);
		GF.push_back(0.10070708E+01);GF.push_back(0.10065769E+01);GF.push_back(0.10041945E+01);GF.push_back(0.10064120E+01);
		GF.push_back(0.10038391E+01);GF.push_back(0.10030688E+01);GF.push_back(0.10066030E+01);GF.push_back(0.10016240E+01);
		GF.push_back(0.10016652E+01);GF.push_back(0.10074447E+01);GF.push_back(0.10076891E+01);GF.push_back(0.10014995E+01);
		GF.push_back(0.10049790E+01);GF.push_back(0.10007124E+01);GF.push_back(0.10020235E+01);GF.push_back(0.10279195E+01);
    
		GI.push_back(0.56000000E+00);GI.push_back(0.50000000E+00);GI.push_back(0.56000000E+00);GI.push_back(0.47999999E+00);
		GI.push_back(0.50000000E+00);GI.push_back(0.44999999E+00);GI.push_back(0.60000002E+00);GI.push_back(0.50000000E+00);
		GI.push_back(0.46000001E+00);GI.push_back(0.55000001E+00);GI.push_back(0.50000000E+00);GI.push_back(0.55000001E+00);
		GI.push_back(0.38000000E+00);GI.push_back(0.50000000E+00);GI.push_back(0.34999999E+00);GI.push_back(0.55000001E+00);
		GI.push_back(0.50000000E+00);GI.push_back(0.44999999E+00);GI.push_back(0.50000000E+00);GI.push_back(0.50000000E+00);
    
		QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);QS.push_back(0.22000000E+01);
		QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);QS.push_back(0.30000000E+01);QS.push_back(0.20000000E+01);
		QS.push_back(0.20000000E+01);QS.push_back(0.30000000E+01);QS.push_back(0.20000000E+01);QS.push_back(0.30000000E+01);
		QS.push_back(0.20000000E+01);QS.push_back(0.50000000E+00);QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);
		QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);QS.push_back(0.20000000E+01);
    
		RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);
		RS.push_back(0.50000000E+00);RS.push_back(0.15000000E+01);RS.push_back(0.10000000E+00);RS.push_back(0.50000000E+00);
		RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);
		RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+01);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);
		RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);RS.push_back(0.50000000E+00);
    
		SM.push_back(0.84969723E+00);SM.push_back(0.84408921E+00);SM.push_back(0.79834670E+00);SM.push_back(0.88117069E+00);
		SM.push_back(0.79869562E+00);SM.push_back(0.84020197E+00);SM.push_back(0.86620516E+00);SM.push_back(0.84974420E+00);
		SM.push_back(0.78492051E+00);SM.push_back(0.90510607E+00);SM.push_back(0.77372539E+00);SM.push_back(0.90547490E+00);
		SM.push_back(0.82473016E+00);SM.push_back(0.76297575E+00);SM.push_back(0.87047338E+00);SM.push_back(0.77146745E+00);
		SM.push_back(0.83532584E+00);SM.push_back(0.86597520E+00);SM.push_back(0.80329478E+00);SM.push_back(0.84098595E+00);
    
		EA.push_back(0.18000000E+02);EA.push_back(0.14000000E+02);EA.push_back(0.14000000E+02);EA.push_back(0.42000000E+02);
		EA.push_back(0.20000000E+02);EA.push_back(0.19000000E+02);EA.push_back(0.30000000E+02);EA.push_back(0.18000000E+02);
		EA.push_back(0.93000002E+01);EA.push_back(0.25000000E+02);EA.push_back(0.10000000E+02);EA.push_back(0.25000000E+02);
		EA.push_back(0.24000000E+02);EA.push_back(0.12000000E+02);EA.push_back(0.20000000E+02);EA.push_back(0.25000000E+02);
		EA.push_back(0.60000000E+02);EA.push_back(0.10000000E+03);EA.push_back(0.60000000E+02);EA.push_back(0.30000000E+02);
    
		CR.push_back(0.45034027E+01);CR.push_back(0.29177928E+01);CR.push_back(0.41500096E+01);CR.push_back(0.34477754E+01);
		CR.push_back(0.37076824E+01);CR.push_back(0.35676460E+01);CR.push_back(0.38610494E+01);CR.push_back(0.31175354E+01);
		CR.push_back(0.16488581E+01);CR.push_back(0.35281672E+01);CR.push_back(0.17087947E+01);CR.push_back(0.26236989E+01);
		CR.push_back(0.28201730E+01);CR.push_back(0.13170183E+01);CR.push_back(0.31558700E+01);CR.push_back(0.23904498E+01);
		CR.push_back(0.41986074E+01);CR.push_back(0.61198459E+01);CR.push_back(0.77948456E+01);CR.push_back(0.39223223E+01);
    
		C1R.push_back(0.16976673E+02);C1R.push_back(0.44948950E+01);C1R.push_back(0.52029748E+01);C1R.push_back(0.31223445E+01);
		C1R.push_back(0.52075500E+01);C1R.push_back(0.39344752E+01);C1R.push_back(0.39021666E+01);C1R.push_back(0.25285969E+01);
		C1R.push_back(0.81879050E+00);C1R.push_back(0.54939213E+01);C1R.push_back(0.18226281E+01);C1R.push_back(0.20840724E+01);
		C1R.push_back(0.22750535E+01);C1R.push_back(0.13973175E+01);C1R.push_back(0.30427899E+01);C1R.push_back(0.19867337E+01);
		C1R.push_back(0.78179336E+01);C1R.push_back(0.19883865E+02);C1R.push_back(0.40253441E+02);C1R.push_back(0.52068253E+01);*/
	};
};

//Constants
const double LomonosovFortovGasModel::R13 = 1.0/3.0;
const double LomonosovFortovGasModel::GGI = 2.0/3.0;

#endif
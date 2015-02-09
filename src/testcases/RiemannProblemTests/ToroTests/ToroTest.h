#ifndef TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest
#define TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest

#include "TestCase.h"
#include "gengrid1D.h"


namespace ToroTests {

//common structure to pass input parameters into test
struct ToroTestsInit{
	int nCells;
	double Lx;					//solve problem in [0; Lx] segment
	double discontinuity_pos;		//posiyion of initial discontinuity
	double TimeMax;				//time of solution

	//Left state density, pressure and velocity
	double roL;
	double pL;
	double uL;

	//Right state density, pressure and velocity
	double roR;
	double pR;
	double uR;
};

//Base class for all test from Toro book (the paragraph 4.3.3 pp. 129 - 133)
class ToroTest : public TestCase {
protected:
	Kernel* _kernel; //Computational kernel object
	Grid _grid;					  //Grid object	
	Configuration _configuration; //Configuration object
	RiemannSolverConfiguration::RiemannSolverType _riemannSolverType;		//type of Riemann Solver Problem
	std::vector<GasModel::ConservativeVariables> U_an;			//analitical solution
	std::string solution_name;

public:
	//Test parameters
	int nCells;
	double Lx;					//solve problem in [0; Lx] segment
	double discontinuity_pos;		//posiyion of initial discontinuity
	double TimeMax;				//time of solution

	//Left state density, pressure and velocity
	double roL;
	double pL;
	double uL;

	//Right state density, pressure and velocity
	double roR;
	double pR;
	double uR;

	//constructors
	ToroTest() {
		solution_name = "solution";
	};
	ToroTest(ToroTestsInit &params) {
		this->SetParams(params);
		solution_name = "solution";
	};

	//set parameters for test 1
	void SetParams(ToroTestsInit &params) {
		//Test constant's
		nCells = params.nCells;
		Lx = params.Lx;
		discontinuity_pos = params.discontinuity_pos;
		TimeMax = params.TimeMax;

		//Left state density, pressure and velocity
		roL = params.roL;
		pL = params.pL;
		uL = params.uL;

		//Right state density, pressure and velocity
		roR = params.roR;
		pR = params.pR;
		uR = params.uR;
	};

	//Prepare computational grid
	void PrepareGrid() {		
		_grid = GenGrid1D(_kernel->getParallelHelper(), nCells, 0.0, Lx, false);
		_kernel->BindGrid(&_grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Godunov;

		//Availible gas models
		_configuration.AddGasModel("Air");			

		//Air (ideal gas)
		_configuration.GasModelsConfiguration["Air"].GasModelName = "PerfectGasModel";
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("IdealGasConstant", 8.3144621);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatPressure", 1006.43);			

		//Boundary conditions				
		_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["left"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["left"].MaterialName = "Air";

		_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["right"].MovementType = BoundaryConditionMovementType::Fixed;
		_configuration.BoundaryConditions["right"].MaterialName = "Air";

		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.1;
		_configuration.RungeKuttaOrder = 1;		

		//ALE settings
		_configuration.ALEConfiguration.ALEMotionType = "Eulerian";		

		//Run settings
		_configuration.MaxIteration = 1000000;
		_configuration.MaxTime = TimeMax;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = 0;			

		_kernel->BindConfiguration(_configuration);	

		return _configuration;
	};

	//Main interface function for running test case code
	void RunTestWithKernel(Kernel* kernel) {		
		_kernel = kernel; //Save reference to kernel
		
		//Get grid
		PrepareGrid();

		//Set parameters
		PrepareConfiguration();
		
		//Initialize calculation
		_kernel->InitCalculation();

		//Initial conditions
		_kernel->GenerateInitialConditions(new TestCaseInitialConditions(uL, roL, pL, 0, uR, roR, pR, 0, Lx, discontinuity_pos));	

		//Run computational cycle
		_kernel->RunCalculation();

		//Finilize
		_kernel->FinalizeCalculation();

		//Check results
		//Output result
		_kernel->SaveGrid(solution_name + ".cgns");
		_kernel->SaveSolution(solution_name + ".cgns", solution_name);
	};

	//Run test with program arguments
	void RunTest(int* argc, char** argv[]) {
		_kernel = new Kernel();
		_kernel->Initilize(argc, argv);

		//Call main function
		RunTestWithKernel(_kernel);

		_kernel->Finalize();
	};

	//Specify initial conditions
	//Sod's shock tube
	class TestCaseInitialConditions : public InitialConditions::InitialConditions
	{
	public:
		//Left state density, pressure, velocity and material index
		double _uLeft;
		double _roLeft;
		double _pLeft;
		int _nmatLeft;

		//Right state density, pressure, velocity and material index
		double _uRight;
		double _roRight;
		double _pRight;
		int _nmatRight;

		//length of domain and initial discontinuity position
		double _lx;
		double _x0;

		//constructor
		TestCaseInitialConditions(double uLeft, double roLeft, double pLeft, int nmatLeft, double uRight, double roRight, double pRight, int nmatRight, double lx, double x0) :
		_uLeft(uLeft),
		_roLeft(roLeft),
		_pLeft(pLeft),
		_nmatLeft(nmatLeft),
		_uRight(uRight),
		_roRight(roRight),
		_pRight(pRight),
		_nmatRight(nmatRight),
		_lx(lx),
		_x0(x0)
		{ };

		virtual int getInitialGasModelIndex(const Cell& cell) {
			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
		
			return 0; 
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
			double roE = 0;

			//TO DO insert a gamma
			if (x <= _x0) {
				ro = _roLeft;
				u = _uLeft;			
				e = _pLeft / (0.4 * ro);
			} else {
				ro = _roRight;
				u = _uRight;						
				e = _pRight / (0.4 * ro);
			};
			
			//Convert to conservative variables						
			roE = ro*(e + (u*u + v*v + w*w) / 2.0);
			initValues.resize(_gasModels[nmat]->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;

			return initValues;
		};
	}; //Initial conditions

	//Exact Solution Block//
	//star values structure
	Godunov3DSolverPerfectGas::StarVariables starValues;

	//fL part of algebraic equation for pressure in exact RP solver	(proposition 4.2.1 from Toro)
	double f1(double pStar) 
	{
		double res = 0;
		double gamma = _configuration.GasModelsConfiguration["Air"].GetPropertyValue("SpecificHeatRatio").first;	//ask EKR

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

	//fR part of algebraic equation for pressure in exact RP solver	(proposition 4.2.1 from Toro)
	double f2(double pStar) 
	{
		double res = 0;
		double gamma = _configuration.GasModelsConfiguration["Air"].GetPropertyValue("SpecificHeatRatio").first;	//ask EKR

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

	//Target function for pressure equation in exact RP solver TO DO ask Ekr about object
	static void fFunction(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *object)
	{
		ToroTest* params = (ToroTest*)object;
		//
		// this callback calculates		
		// f(pStar) = f1(pStar, WL) + f2(pStar, WR) + deltaU
		//
		double pStar = x[0];
		double deltaU = params->uR - params->uL;
		double res = params->f1(pStar) + params->f2(pStar) + deltaU;

		fi[0] = res;		
	}

	//compute variables in the star region
	Godunov3DSolverPerfectGas::StarVariables ComputeStarVariables(double eps) {				

		//Normal velocity difference
		double deltaU = uR - uL;

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
		double pStar = 0.5*(pL + pR);
		x[0] = pStar;

		//Iterative scheme parameters
		double epsg = eps;
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
		starValues.pStar = x[0];		

		//Compute star region velocity
		starValues.uStar = 0.5*(uL + uR) + 0.5*(f2(starValues.pStar) - f1(starValues.pStar));

		//Determine nonlinear wave type and properties
		double gamma = _configuration.GasModelsConfiguration["Air"].GetPropertyValue("SpecificHeatRatio").first;	//ask EKR
		double C = (gamma - 1.0) / (gamma + 1.0);
		starValues.MaxSpeed = 0;

		//Left side of contact
		double pRatioL = starValues.pStar / pL;
		double aL = sqrt(gamma * pL / roL);
		if (starValues.pStar > pL) {
			//Left shock
			starValues.leftWave = Godunov3DSolverPerfectGas::Shock;

			//Determine density in star region			
			starValues.roStarL = pRatioL + C;
			starValues.roStarL /= C * pRatioL + 1.0;
			starValues.roStarL *= roL;

			//Determine shock propagation speed			
			starValues.SL = sqrt((gamma + 1) * pRatioL / (2*gamma) + (gamma - 1) / (2*gamma));
			starValues.SL = uL - aL * starValues.SL;

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SL));
		} else {
			//Left rarefaction
			starValues.leftWave = Godunov3DSolverPerfectGas::Rarefaction;

			//Determine density in star region
			starValues.roStarL = roL * pow(pRatioL , 1.0 / gamma);

			//Determine rarefaction head propagation speed		
			starValues.SHL = uL - aL;

			//Determine rarefaction tail propagation speed		
			double aStarL = aL * pow(pRatioL, (gamma - 1) / (2*gamma));
			starValues.STL = starValues.uStar - aStarL;		

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SHL));
		};

		//Right side of contact
		double pRatioR = starValues.pStar / pR;
		double aR = sqrt(gamma * pR / roR);
		if (starValues.pStar > pR) {
			//Right shock
			starValues.rightWave = Godunov3DSolverPerfectGas::Shock;

			//Determine density in star region			
			starValues.roStarR = pRatioR + C;
			starValues.roStarR /= C * pRatioR + 1.0;
			starValues.roStarR *= roR;

			//Determine shock propagation speed			
			starValues.SR = sqrt((gamma + 1) * pRatioR / (2*gamma) + (gamma - 1) / (2*gamma));
			starValues.SR = uL + aR * starValues.SR;

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SR));
		} else {
			//Left rarefaction
			starValues.rightWave = Godunov3DSolverPerfectGas::Rarefaction;

			//Determine density in star region
			starValues.roStarR = roR * pow(pRatioR , 1.0 / gamma);

			//Determine rarefaction head propagation speed		
			starValues.SHR = uL + aR;

			//Determine rarefaction tail propagation speed		
			double aStarR = aR * pow(pRatioL, (gamma - 1) / (2*gamma));
			starValues.STR = starValues.uStar + aStarR;			

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SHR));
		};
				
		return starValues;
	};

	//compute ro u and p in x point at moment t ([0] - r0, [1] - u, [2] - p)
	std::vector<double> ComputeExactValuesInCell(double x, double t) {
		//Compute flux (Toro p. 219) 
		//Sample exact solution at S = x/t
		double S = x/t;
		double ro;
		double u;
		double p;
		double gamma = _configuration.GasModelsConfiguration["Air"].GetPropertyValue("SpecificHeatRatio").first;	//ask EKR
		
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

		std::vector<double> result;
		result.push_back(ro);
		result.push_back(u);
		result.push_back(p);
		return result;
	};
	
	void ComputeExactSolution(double eps = 1.0e-10) {
		//Add specifik heat ratio of our gas
		_configuration.AddGasModel("Air");
		_configuration.GasModelsConfiguration["Air"].SetPropertyValue("SpecificHeatRatio", 1.4);

		//compute star variables at first		
		starValues = ComputeStarVariables(eps);

		std::string exact_sol_file = "exact_solution.dat";
		std::ofstream ofs(exact_sol_file);
		for (int cellIndex = 0; cellIndex < nCells; cellIndex++) {
			double x;
			double hx = Lx/nCells;	//grid step for UNIFORM GRID
			x = hx*(0.5 + cellIndex);
			//compute ro u and p in x point at moment t
			std::vector<double> res =  ComputeExactValuesInCell(x - discontinuity_pos, TimeMax);
			ofs << x << ' ' << res[0] << ' ' << res[1] << ' ' << res[2] << '\n';
		};
		ofs.close();
	};

}; //TestCase

}; //namespace


#endif
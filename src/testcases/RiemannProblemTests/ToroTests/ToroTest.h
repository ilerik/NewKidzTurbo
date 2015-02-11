#ifndef TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest
#define TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest

#include "TestCase.h"
#include "gengrid1D.h"


namespace ToroTests {

//common structure to pass input parameters into test
struct ToroTestsInit{
	int nCells;
	double Lx;						//solve problem in [0; Lx] segment
	double discontinuityPosition;	//position of initial discontinuity
	double TimeMax;					//time of solution

	//Left state density, pressure and velocity
	double gammaL;
	double roL;
	double pL;
	double uL;

	//Right state density, pressure and velocity
	double gammaR;
	double roR;
	double pR;
	double uR;
};

//Class for storing results of testing
struct RiemanProblemTestReport {	
	double LInf; //Infinity norm of error
	double L2;	 //L2 norm of error
};

//Base class for all test from Toro book (the paragraph 4.3.3 pp. 129 - 133)
class ToroTest : public TestCase {
protected:	
	Configuration _configuration;											//Configuration object
	RiemannSolverConfiguration::RiemannSolverType _riemannSolverType;		//Type of Riemann Solver Problem	
	std::string _solutionName;												//Solution name
	double _eps;															//Iterative solver precision

	//Test settings
	ToroTestsInit _initSettings;
public:	
	//Test parameters
	int nCells;
	double Lx;					//solve problem in [0; Lx] segment
	double discontinuity_pos;	//position of initial discontinuity
	double TimeMax;				//time of solution

	//Left state density, pressure and velocity
	double gammaL;
	double roL;
	double pL;
	double uL;

	//Right state density, pressure and velocity
	double gammaR;
	double roR;
	double pR;
	double uR;

	//Constructor with test parameters
	ToroTest(ToroTestsInit &params) {
		_initSettings = params;
		this->SetParams(_initSettings);
		_eps = 1e-10;
		_solutionName = "Solution";
	};

	//set parameters for test 1
	void SetParams(ToroTestsInit &params) {
		//Test constant's
		nCells = params.nCells;
		Lx = params.Lx;
		discontinuity_pos = params.discontinuityPosition;
		TimeMax = params.TimeMax;

		//Left state density, pressure and velocity
		gammaL = params.gammaL;
		roL = params.roL;
		pL = params.pL;
		uL = params.uL;

		//Right state density, pressure and velocity
		gammaR = params.gammaR;
		roR = params.roR;
		pR = params.pR;
		uR = params.uR;
	};

	//Prepare computational grid
	void PrepareGrid() {		
		Grid& grid = GenGrid1D(_kernel->getParallelHelper(), nCells, 0.0, Lx, false);
		_kernel->BindGrid(&grid);
	};

	//Prepare configuration object and set all parameters
	Configuration& PrepareConfiguration() {
		//Hardcode configuration for now
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = _riemannSolverType;

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

	//Main interface function for running test case code (grid already defined)
	virtual void RunTestWithKernel(Kernel* kernel) override {		
		_kernel = kernel; //Save reference to kernel

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

		//Output result		
		_kernel->SaveSolution("result.cgns", _solutionName);
	};

	//Run test with program arguments
	void RunTest(int* argc, char** argv[]) {
		//Initialize kernel
		_kernel = new Kernel();
		_kernel->Initilize(argc, argv);

		//Prepare computational grid
		PrepareGrid();

		//Call main function
		RunTestWithKernel(_kernel);

		//Finilize
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

		if (pStar > pL) {
			//Left shock
			double AL = 2 / ((gammaL + 1.0) * roL);
			double BL = pL * (gammaL - 1.0) / (gammaL + 1.0);
			res = (pStar - pL) * sqrt(AL / (pStar + BL));
		} else {
			//Left rarefaction
			double aL = sqrt(gammaL * pL / roL);
			res = pow(pStar/pL, (gammaL - 1.0)/(2*gammaL)) - 1.0;
			res *= 2*aL/(gammaL - 1.0);
		};

		return res;
	};

	//fR part of algebraic equation for pressure in exact RP solver	(proposition 4.2.1 from Toro)
	double f2(double pStar) 
	{
		double res = 0;		

		if (pStar > pR) {
			//Right shock
			double AR = 2 / ((gammaR + 1.0) * roR);
			double BR = pR * (gammaR - 1.0) / (gammaR + 1.0);
			res = (pStar - pR) * sqrt(AR / (pStar + BR));
		} else {
			//Right rarefaction
			double aR = sqrt(gammaR * pR / roR);
			res = pow(pStar/pR, (gammaR - 1.0)/(2*gammaR)) - 1.0;
			res *= 2*aR/(gammaR - 1.0);
		};

		return res;
	};

	//Target function callback for pressure equation in exact RP solver
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
		starValues.MaxSpeed = 0;

		//Left side of contact
		double pRatioL = starValues.pStar / pL;
		double aL = sqrt(gammaL * pL / roL);
		double CL = (gammaL - 1.0) / (gammaL + 1.0);
		if (starValues.pStar > pL) {
			//Left shock
			starValues.leftWave = Godunov3DSolverPerfectGas::Shock;

			//Determine density in star region			
			starValues.roStarL = pRatioL + CL;
			starValues.roStarL /= CL * pRatioL + 1.0;
			starValues.roStarL *= roL;

			//Determine shock propagation speed			
			starValues.SL = sqrt((gammaL + 1) * pRatioL / (2*gammaL) + (gammaL - 1) / (2*gammaL));
			starValues.SL = uL - aL * starValues.SL;

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SL));
		} else {
			//Left rarefaction
			starValues.leftWave = Godunov3DSolverPerfectGas::Rarefaction;

			//Determine density in star region
			starValues.roStarL = roL * pow(pRatioL , 1.0 / gammaL);

			//Determine rarefaction head propagation speed		
			starValues.SHL = uL - aL;

			//Determine rarefaction tail propagation speed		
			double aStarL = aL * pow(pRatioL, (gammaL - 1) / (2*gammaL));
			starValues.STL = starValues.uStar - aStarL;		

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SHL));
		};

		//Right side of contact
		double pRatioR = starValues.pStar / pR;
		double aR = sqrt(gammaR * pR / roR);
		double CR = (gammaR - 1.0) / (gammaR + 1.0);
		if (starValues.pStar > pR) {
			//Right shock
			starValues.rightWave = Godunov3DSolverPerfectGas::Shock;

			//Determine density in star region			
			starValues.roStarR = pRatioR + CR;
			starValues.roStarR /= CR * pRatioR + 1.0;
			starValues.roStarR *= roR;

			//Determine shock propagation speed			
			starValues.SR = sqrt((gammaR + 1) * pRatioR / (2*gammaR) + (gammaR - 1) / (2*gammaR));
			starValues.SR = uL + aR * starValues.SR;

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SR));
		} else {
			//Left rarefaction
			starValues.rightWave = Godunov3DSolverPerfectGas::Rarefaction;

			//Determine density in star region
			starValues.roStarR = roR * pow(pRatioR , 1.0 / gammaR);

			//Determine rarefaction head propagation speed		
			starValues.SHR = uL + aR;

			//Determine rarefaction tail propagation speed		
			double aStarR = aR * pow(pRatioL, (gammaR - 1) / (2*gammaR));
			starValues.STR = starValues.uStar + aStarR;			

			starValues.MaxSpeed = max(starValues.MaxSpeed, std::abs(starValues.SHR));
		};
				
		return starValues;
	};

	//compute ro u and p in x point at moment t ([0] - r0, [1] - u, [2] - p) star variables must be computed in advance
	std::vector<double> ComputeExactValuesInCell(double x, double t) {
		//Compute flux (Toro p. 219) 
		//Sample exact solution at S = x/t
		double S = x/t;
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
					double aL = sqrt(gammaL * pL / roL);
					double CL = 2 / (gammaL + 1) + (gammaL-1) * (uL - S) / ((gammaL + 1) * aL);
					CL = pow(CL, 2 / (gammaL - 1));

					//Density
					ro = CL * roL;

					//Velocity
					u = aL + (gammaL - 1) * uL/ 2.0 + S;
					u *= 2 / (gammaL + 1);

					//Pressure					
					p = CL * pL;
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
					double aR = sqrt(gammaR * pR / roR);
					double CR = 2 / (gammaR + 1) - (gammaR-1) * (uR - S) / ((gammaR + 1) * aR);
					CR = pow(CR, 2 / (gammaR - 1));

					//Density
					ro = CR * roR;

					//Velocity
					u = -aR + (gammaR - 1) * uR/ 2.0 + S;
					u *= 2 / (gammaR + 1);

					//Pressure					
					p = CR * pR;
				};
			};
		};

		std::vector<double> result;
		result.push_back(ro);
		result.push_back(u);
		result.push_back(p);
		return result;
	};

	//compute conservative variables at some point, star variables must be computed in advance
	std::vector<double> ComputeExactConservativeValuesInCell(Vector r, double time) {
		double x = r.x - discontinuity_pos;
		double gamma = _configuration.GasModelsConfiguration["Air"].GetPropertyValue("SpecificHeatRatio").first;	//ask EKR
		std::vector<double> result = ComputeExactValuesInCell(x, time);
		std::vector<double> U(5);
		double ro = U[0] = result[0];
		double u = result[1];
		double P = result[2];
		U[1] = ro * u;
		U[2] = 0;
		U[3] = 0;
		U[4] = P/(gamma - 1.0) + ro * u * u / 2.0;
		return U;
	};			

	//Calculate and write exact solution to cgns
	void WriteExactSolution(std::string filename, std::string solutionName) {
		starValues = ComputeStarVariables(_eps);
		std::function<std::vector<double>(Vector)> FVariables = std::bind(&ToroTest::ComputeExactConservativeValuesInCell, this, std::placeholders::_1, TimeMax);
		std::function<int(Vector)> FMaterial = [](Vector r) { return 0; };
		_kernel->FillInitialConditions(FVariables, FMaterial);
		_kernel->SaveSolution(filename, solutionName);
	};

	//Calculate various error measures
	RiemanProblemTestReport CompareWithAnalyticalSolution() {
		//Result structure initialization
		RiemanProblemTestReport report;
		report.L2 = 0;
		report.LInf = 0;

		//Star variables
		starValues = ComputeStarVariables(_eps);

		//Compute local part of error for each cell
		for (int cellIndex = 0; cellIndex < _kernel->GetLocalCellsNumber(); cellIndex++) {
			//Obtain current solution values
			double x = _kernel->GetLocalCell(cellIndex).CellCenter.x;
			double pFact = _kernel->GetCellPressure(cellIndex);
			//double uFact = _kernel->GetCellVelocityX(cellIndex);
			//double roFact = _kernel->GetCellDensity(cellIndex);

			//Compute exact solution
			std::vector<double> exactValues = ComputeExactValuesInCell(x - _initSettings.discontinuityPosition, TimeMax);
			double pExact = exactValues[2];
			double uExact = exactValues[1];
			double roExact = exactValues[0];

			//Compute particular errors
			double dp = pFact - pExact;
			double dx = _kernel->GetLocalCell(cellIndex).CellVolume;
			report.L2 += std::pow(dp * dx, 2.0);
			report.LInf = std::max(std::abs(dp), report.LInf);
		};

		//Syncronize and aggregate
		report.L2 = _kernel->getParallelHelper()->SumDouble(report.L2);
		report.L2 = std::sqrt(report.L2);
		report.LInf = _kernel->getParallelHelper()->Max(report.LInf);

		//Return result
		return report;
	};

	//Test result function
	virtual TestCaseResultInfo GetTestCaseResultInfo() override //Get results of test run 
	{
		RiemanProblemTestReport report = CompareWithAnalyticalSolution();
		std::cout<<"L2norm = "<<report.L2<<" "<<"LInfnorm = "<<report.LInf<<"\n";

		TestCaseResultInfo result;
		result.testStatus = TestCaseResultInfo::TestStatusType::Passed;
		return result;
	};

}; //TestCase

}; //namespace


#endif
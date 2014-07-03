#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "cgnsload.h"
#include "grid.h"
#include "gridgeneration.h"
#include "model.h"
#include "riemannsolvers.h"
#include "garbarukmodel.h"
#include "SAmodel.h"
#include "biffplatemodel.h"
#include "incompress2D_Flow.h"
#include "triangleload.h"

template< typename T >
std::string int_to_hex( T i )
{
  std::stringstream stream;
  stream << "0x" 
         << std::setfill ('0') << std::setw(sizeof(T)*2) 
         << std::hex << i;
  return stream.str();
}

//WID отладочные функции
void ComputeTVisc(char* file_name)
{
	function VEL;
	VEL.ReadData1(file_name);
	double K = 0.41;
	double A = 12.0;
	double rho = 1.12;
	double vis_d = 1.7894e-03;
	double vis_k = vis_d/rho;
	double tau_w = vis_d*VEL.value[1]/VEL.x[1];
	double U_fr = sqrt(tau_w);
	double h = VEL.x[VEL.size-1];
	double U_inf = VEL.value[VEL.size-1];
	function DT = (1.0 - VEL/U_inf);
	double delta = (1.0 - VEL/U_inf).TrapeziumIntegral();

	function t_vis = VEL;
	for(int i=0; i<t_vis.size; i++)
	{
		double D = (1.0 - exp((-1.0)*U_fr*t_vis.x[i]/(vis_k*A)));
		D = D*D*D;
		double gamma = t_vis.x[i]/h;
		gamma = 1.0 + 5.5*gamma*gamma*gamma*gamma*gamma*gamma;
		gamma = 1/gamma;
		t_vis.value[i] = K*U_fr*min(t_vis.x[i]*D, delta*gamma);
	};
	t_vis.WriteData("t_vis_new.dat");
};
void TestSolver()
{
	BiffTurbulentPlateModel<Roe3DSolverPerfectGas> model(1375.0, 7562500.0, 2.0);
	model.SetViscosity(1.7894e-03);

	function vel;
	BM_FUN_t_visc.ReadData1("t_visc.dat");
	vel.ReadData1("t_velocity.dat");
	BM_FUN_t_dif.ReadData1("t_visc.dat");
	BM_FUN_lh_vel_der.ReadData1("t_LeeHarsha.dat");

	BM_FUN_t_dif_der = BM_FUN_t_dif.Derivate();
	BM_FUN_vel_der = vel.Derivate();
	BM_FUN_lh_vel_der = BM_FUN_lh_vel_der.Derivate();
	
	/*int N = 100;
	std::vector<double> subgrid(N+1);
	double h = vel.x[vel.size-1];
	for(int i=0; i<subgrid.size(); i++) subgrid[i] = i*h/N;
	BM_FUN_t_visc = BM_FUN_t_visc.BindToNewGridFO(subgrid);
	BM_FUN_vel_der = BM_FUN_vel_der.BindToNewGridFO(subgrid);
	BM_FUN_t_dif = BM_FUN_t_dif.BindToNewGridFO(subgrid);
	BM_FUN_t_dif_der = BM_FUN_t_dif_der.BindToNewGridFO(subgrid);
	BM_FUN_lh_vel_der = BM_FUN_lh_vel_der.BindToNewGridFO(subgrid);*/

	function res = model.SolverTest(vel);
	res.WriteData("tau_un_test.dat");
	std::getchar();
};

//Reimann Problem solver tests
ConservativeVariables SODTestInitDistribution(Vector r, void *par) {	
	//Gamma
	double gamma = 1.4;

	//Left state x <= 0.5
	double roL = 1.0;
	double pL = 1.0;
	Vector vL = Vector(0,0,0);		
	ConservativeVariables UL;
	UL.ro = roL;
	UL.rou = roL * vL.x;
	UL.rov = roL * vL.y;
	UL.row = roL * vL.z;
	UL.roE = pL/(gamma-1) + roL * vL.mod() * vL.mod() / 2.0;
	//Right state x > 0.5
	double roR = 0.125;
	double pR = 0.1;
	Vector vR = Vector(0,0,0);
	ConservativeVariables UR;
	UR.ro = roR;
	UR.rou = roR * vR.x;
	UR.rov = roR * vR.y;
	UR.row = roR * vR.z;
	UR.roE = pR/(gamma-1) + roR * vR.mod() * vR.mod() / 2.0;

	if (r.x <= 0.5) {
		return UL;
	} else {
		return UR;
	};
};
ConservativeVariables ToroTestInitDistribution1(Vector r, void *par) {	
	//Gamma
	double gamma = 1.4;

	//Left state x <= 0.3
	double roL = 1.0;
	double pL = 1.0;
	Vector vL = Vector(0.75,0,0);		
	ConservativeVariables UL;
	UL.ro = roL;
	UL.rou = roL * vL.x;
	UL.rov = roL * vL.y;
	UL.row = roL * vL.z;
	UL.roE = pL/(gamma-1) + roL * vL.mod() * vL.mod() / 2.0;
	//Right state x > 0.3
	double roR = 0.125;
	double pR = 0.1;
	Vector vR = Vector(0,0,0);
	ConservativeVariables UR;
	UR.ro = roR;
	UR.rou = roR * vR.x;
	UR.rov = roR * vR.y;
	UR.row = roR * vR.z;
	UR.roE = pR/(gamma-1) + roR * vR.mod() * vR.mod() / 2.0;

	if (r.x <= 0.3) {
		return UL;
	} else {
		return UR;
	};
};
//на первом порядке разваливается
ConservativeVariables ToroTestInitDistribution2(Vector r, void *par) {	
	//Gamma
	double gamma = 1.4;

	//Left state x <= 0.5
	double roL = 1.0;
	double pL = 0.4;
	Vector vL = Vector(-2.0,0,0);		
	ConservativeVariables UL;
	UL.ro = roL;
	UL.rou = roL * vL.x;
	UL.rov = roL * vL.y;
	UL.row = roL * vL.z;
	UL.roE = pL/(gamma-1) + roL * vL.mod() * vL.mod() / 2.0;
	//Right state x > 0.5
	double roR = 1.0;
	double pR = 0.4;
	Vector vR = Vector(2.0,0,0);
	ConservativeVariables UR;
	UR.ro = roR;
	UR.rou = roR * vR.x;
	UR.rov = roR * vR.y;
	UR.row = roR * vR.z;
	UR.roE = pR/(gamma-1) + roR * vR.mod() * vR.mod() / 2.0;

	if (r.x <= 0.5) {
		return UL;
	} else {
		return UR;
	};
};
ConservativeVariables ToroTestInitDistribution4(Vector r, void *par) {	
	//Gamma
	double gamma = 1.4;

	//Left state x <= 0.4
	double roL = 5.99924;
	double pL = 460.894;
	Vector vL = Vector(19.5975, 0, 0);		
	ConservativeVariables UL;
	UL.ro = roL;
	UL.rou = roL * vL.x;
	UL.rov = roL * vL.y;
	UL.row = roL * vL.z;
	UL.roE = pL/(gamma-1) + roL * vL.mod() * vL.mod() / 2.0;
	//Right state x > 0.4
	double roR = 5.99242;
	double pR = 46.095;
	Vector vR = Vector(-6.19633,0,0);
	ConservativeVariables UR;
	UR.ro = roR;
	UR.rou = roR * vR.x;
	UR.rov = roR * vR.y;
	UR.row = roR * vR.z;
	UR.roE = pR/(gamma-1) + roR * vR.mod() * vR.mod() / 2.0;

	if (r.x <= 0.4) {
		return UL;
	} else {
		return UR;
	};
};

bool RunRiemannProblemTest(ConservativeVariables(*funcInitValue)(Vector, void *), int N_cells, double time, int order) {
	Model<Roe3DSolverPerfectGas> model;
	Grid grid;

	//Shock tube problem setting
	double lBegin = 0;
	double lEnd = 1.0;
	Vector direction = Vector(1,0,0);
	grid = GenGrid1D(N_cells, lBegin, lEnd, direction);

	//grid = GenGrid2D(100, 10, 1.0, 1.0, 1.0, 1.0);
	
	//Set fluid properties	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.DisableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.005);
	model.SetSchemeOrder(order);

	//Bind computational grid
	model.BindGrid(grid);
	if(order > 1) model.EnableLimiter();
	else model.DisableLimiter();

	//Set initial conditions
	model.SetInitialConditions(funcInitValue);

	//free boundary condition
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);
	model.SetBoundaryCondition("left", NaturalBC);
	model.SetBoundaryCondition("right", NaturalBC);

	model.SaveToTechPlot("RP_Test_Init.dat");

	//Total time
		for (int i = 0; i < 200000000; i++) {
		model.Step();
		if (i % 1 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		model.SaveToTechPlot("RP_Test.dat");
		if (model.totalTime > time) break;
	};

	//Save result to techplot
	model.SaveToTechPlot("RP_Test.dat");
	return true;
};

struct PoiseuilleSettingStruct {
	double PIn;
	double POut;
	double YMin;
	double YMax;
	double XMin;
	double XMax;
	double Viscosity;
	double Temperature;
	double Gamma;
	double InitDensity;
};

Vector PoiseuilleVelocityDistribution(Vector position, void *params) {
	Vector velocity(0,0,0);
	PoiseuilleSettingStruct info = *(PoiseuilleSettingStruct *)params;

	//Velocity profile	
	double L = info.XMax - info.XMin;
	double deltaP = info.PIn - info.POut;	
	double R = (info.YMax - info.YMin) / 2.0;
	double YMid = (info.YMax + info.YMin) / 2.0;
	double r = abs(position.y - YMid);
	//velocity.x = deltaP * ( R*R - r*r)/ ( 4 * info.Viscosity * L);	//for tube
	velocity.x = deltaP * ( R*R - r*r)/ ( 2.0 * info.Viscosity * L);	//for channel

	return velocity;
};

ConservativeVariables PoiseuilleTestInitDistribution(Vector r, void *params) {	
	PoiseuilleSettingStruct info = *(PoiseuilleSettingStruct *)params;


	//Compute conservative variables
	double ro = 1.0;
	double p = 1.0;
	Vector v = PoiseuilleVelocityDistribution(r, params);

	ConservativeVariables U;
	U.ro = info.InitDensity;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = info.PIn/(info.Gamma-1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};

struct ShearFlowSettingStruct {
	double P;
	double U_inf;
	double YMin;
	double YMax;
	double XMin;
	double XMax;
	double Temperature;
	double Density;
	double e;
};

Vector ShearFlowVelocityDistribution(Vector position, void *params) {
	Vector velocity(0,0,0);
	ShearFlowSettingStruct info = *(ShearFlowSettingStruct *)params;

	double Ly = info.YMax - info.YMin;
	//Velocity profile	
	velocity.x = info.U_inf*position.y/Ly;	//for channel

	return velocity;
};

ConservativeVariables ShearFlowInitDistribution(Vector r, void *params) {	
	ShearFlowSettingStruct info = *(ShearFlowSettingStruct *)params;
	
	//Compute conservative variables
	Vector v = ShearFlowVelocityDistribution(r, params);

	ConservativeVariables U;
	U.ro = info.Density;
	U.rou = info.Density * v.x;
	U.rov = info.Density * v.y;
	U.row = info.Density * v.z;
	U.roE = info.Density*info.e + info.Density * v.mod() * v.mod() / 2.0;
	
	return U;	
};

bool RunPoiseuilleTest() {
	Model<Roe3DSolverPerfectGas> model;
	Grid grid;

	//Poiseuille problem setting
	PoiseuilleSettingStruct info;
	double lBegin = 0;
	double lEnd = 1.0;
	info.PIn = 100020;
	info.POut = 100000;
	info.YMin = 0.0;
	info.YMax = 1.0;
	info.Viscosity = 0.5;
	info.XMin = 0.0;
	info.XMax = 1.0;		
	info.Temperature = 300.0;
	info.Gamma = 1.4;		

	grid = GenGrid2D(50, 50, info.XMax, info.YMax, 1.0, 1.0);
	
	//Set fluid properties	
	model.SetGamma(info.Gamma);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.SetViscosity(info.Viscosity);
	model.SetThermalConductivity(0.0);
	model.EnableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.00);

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	ConservativeVariables inletValues = model.PrimitiveToConservativeVariables(Vector(0,0,0), info.PIn, info.Temperature, model.medium);
	info.InitDensity = inletValues.ro;
	model.SetInitialConditions(PoiseuilleTestInitDistribution, &info);

	//Boundary conditions

	//Inlet	
	Model<Roe3DSolverPerfectGas>::SubsonicInletBoundaryCondition InletBC(model);	
	InletBC.setParams(info.PIn, info.Temperature, Vector(0,0,0));
	InletBC.setVelocityDistribution( PoiseuilleVelocityDistribution, &info);
	model.SetBoundaryCondition("left", InletBC);

	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);	
	OutletBC.setParams(info.POut);
	model.SetBoundaryCondition("right", OutletBC);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
	model.SetBoundaryCondition("top", NoSlipBC);
	model.SetBoundaryCondition("bottom", NoSlipBC);	

	//model.LoadSolution("Poiseuille.txt" );
	model.ComputeGradients();
	model.SaveToTechPlot("PoiseuilleInit.dat");

	//Total time
	double maxTime = 0.2;
	for (int i = 0; i < 200000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";			
		};		
		if ((i % 100 == 0) && (i != 0)) {
			model.SaveToTechPlot("Poiseuille.dat");
			model.SaveSolution("Poiseuille.txt" );
		};
		if (model.totalTime > 0.9) break;
	};

	//Save result to techplot
	model.SaveToTechPlot("Poiseuille.dat");
	return true;
};

bool RunShearFlowTest() {
	Model<Roe3DSolverPerfectGas> model;
	Grid grid;

	//Problem setting
	ShearFlowSettingStruct info;
	double lBegin = 0;
	double lEnd = 1.0;
	info.P = 100000;
	info.YMin = 0.0;
	info.YMax = 1.0;
	info.XMin = 0.0;
	info.XMax = 1.0;		
	info.Temperature = 300.0;
	info.U_inf = 2.0;

	grid = GenGrid2D(50, 20, info.XMax, info.YMax, 1.0, 0.95);
	
	//Set fluid properties	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / model.medium.Gamma);
	model.SetMolecularWeight(28.966);
	model.SetViscosity(0.5);
	model.SetThermalConductivity(0.0);
	model.EnableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.00);

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	info.e = model.medium.Cv * info.Temperature;
	info.Density = info.P/(info.e*(model.medium.Gamma - 1.0));		
	model.SetInitialConditions(ShearFlowInitDistribution, &info);

	//Boundary conditions

	//Inlet	
	Model<Roe3DSolverPerfectGas>::SubsonicInletBoundaryCondition InletBC(model);	
	InletBC.setParams(info.P, info.Temperature, Vector(0,0,0));
	InletBC.setVelocityDistribution( ShearFlowVelocityDistribution, &info);
	model.SetBoundaryCondition("left", InletBC);

	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);	
	OutletBC.setParams(info.P);
	model.SetBoundaryCondition("right", OutletBC);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
	model.SetBoundaryCondition("bottom", NoSlipBC);	

	//top boundary condition
	Vector Velocity(info.U_inf, 0, 0);
	Model<Roe3DSolverPerfectGas>::ConstantVelocityBoundaryCondition TopBC(model, Velocity);
	model.SetBoundaryCondition("top", TopBC);
	
	//model.LoadSolution("Poiseuille.txt" );
	model.ComputeGradients();
	model.SaveToTechPlot("ShearFlowInit.dat");

	//Total time
	double maxTime = 0.2;
	for (int i = 0; i < 200000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";			
		};		
		if (i % 100 == 0) {
			model.SaveToTechPlot("ShearFlow.dat");
			model.SaveSolution("ShearFlow.txt" );
		};
		if (model.totalTime > 0.821) break;
	};

	//Save result to techplot
	model.SaveToTechPlot("ShearFlow.dat");
	return true;
};

void RunGAWCalculation() {
	//Load cgns grid	
	std::string solFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Projects\\Grids\\GAW-1\\solution.cgns"; 
	Grid grid = LoadCGNSGrid(solFile);

	// Initialize medium model and place boundary conditions
	Model<Roe3DSolverPerfectGas> model;

	//Set fluid properties
	//Air
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.SetViscosity(1.7894e-03);
	model.SetThermalConductivity(0.0242);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);

	//Bind computational grid
	model.BindGrid(grid);

	//Set wall boundaries and compute wall distances		
	model.SetWallBoundary("gaw1", true);
	model.ComputeWallDistances();
	model.DistanceSorting();

	//Enable viscous part
	model.EnableViscous();

	//Load solution
	model.ReadSolutionFromCGNS(solFile);

	//Output solution to tecplot
	model.SaveWallPressureDistribution("GAW1P.dat");
	model.SaveToTecplot25D("GAW1.dat");

	//Boundary conditions
	//Inlet boundary
	Vector velocity = 70 * Vector(0.9877,-0.1564,0); // 9 degrees
	double pressure = 101579;
	double temperature = 300.214;
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);	

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("inlet", InletBC);
	model.SetBoundaryCondition("outlet", OutletBC);
	model.SetBoundaryCondition("gaw1", NoSlipBC);	
	model.SetBoundaryCondition("sym1", SymmetryBC);
	model.SetBoundaryCondition("sym2", SymmetryBC);
	model.SetBoundaryCondition("top", NoSlipBC);
	model.SetBoundaryCondition("bottom", NoSlipBC);
	return;
};

void RunSAFlatPlate() {
	//Load cgns grid		
	std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineSA.cgns";
	Grid grid = LoadCGNSGrid(solutionFile);

	//Initialize model and bind grid
	SAModel<Roe3DSolverPerfectGas> model;
	model.BindGrid(grid);
	model.Init();

	//Set fluid properties
	//Air
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.SetViscosity(1.7894e-05);
	model.SetThermalConductivity(0.0242);

	//Set wall boundaries and compute wall distances		
	model.SetWallBoundary("plate", true);
	model.ComputeWallDistances();
	model.DistanceSorting();

	//Read initial solution from CGNS
	model.ReadSolutionFromCGNS(solutionFile);		

	//Set boundary conditions			
	Vector velocity(70,0,0);
	double pressure = 101579;
	double temperature = 300.214;
	ConservativeVariables inletValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);				

	//Boundary conditions
	//Inlet boundary
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);	

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("inlet", InletBC);
	model.SetBoundaryCondition("outlet", OutletBC);
	model.SetBoundaryCondition("plate", NoSlipBC);	
	model.SetBoundaryCondition("symmetry", SymmetryBC);
	model.SetBoundaryCondition("top_left", SymmetryBC);
	model.SetBoundaryCondition("top_right", SymmetryBC);

	//Inlet turbulence level
	model.SetBoundarySAValue("inlet", 0);
	model.SaveToTechPlot("solSA.dat");
	//model.ComputeWallVariables();

	//Run some SA steps with frozen flow field
	//model.SAStep();

	return;
};

void RunBiffFlatPlate() {
	std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineSA.cgns";
	//std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineLam.cgns";
	//std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns";	
	Grid grid = LoadCGNSGrid(solutionFile);

	// Initialize medium model and place boundary conditions
	BiffTurbulentPlateModel<Roe3DSolverPerfectGas> model(320000.0, 8337150864.0, 2.0);

	//Set fluid properties
	//Air
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);

	model.SetViscosity(1.7894e-05);	
	model.SetThermalConductivity(0.0242);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);

	////Bind computational grid
	model.BindGrid(grid);			

	////Initial conditions
	//Read solution from CGNS
	model.ReadSolutionFromCGNS(solutionFile);	
	model.LoadSolution("C:\\Users\\Erik\\Dropbox\\Science\\!Projects\\solution.txt");

	//Boundary conditions
	//Inlet boundary
	Vector velocity(70,0,0);
	double pressure = 101579;
	double temperature = 300.214;	
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("inlet", InletBC);
	model.SetBoundaryCondition("outlet", OutletBC);
	model.SetBoundaryCondition("plate", NoSlipBC);	
	model.SetBoundaryCondition("symmetry", SymmetryBC);
	model.SetBoundaryCondition("top_left", SymmetryBC);
	model.SetBoundaryCondition("top_right", SymmetryBC);

	//Set wall boundaries		
	model.SetWallBoundary("plate", true);
	model.ComputeWallDistances();
	model.DistanceSorting();
	model.EnableViscous();	

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init2.dat");

	//Output dimensionless profile for point around x = 1.0
	model.ComputeWallVariables();
	

	std::string outputSolutionFile = "solution2";	

	//Run simulation
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
			model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
			model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
		};
		if (model.totalTime > 10000) break;
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
		model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
		model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
	};
};

void RunBlasiusTest(){	

	//Load cgns grid			
	Grid grid = GenGrid2D(30, 40, 1.2, 0.5, 1.0, 1.0);

	// Initialize medium model and place boundary conditions
	Model<Roe3DSolverPerfectGas> model;

	//Set fluid properties
	//Air
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);

	//model.SetViscosity(1.7894e-05);
	model.SetViscosity(1.7894e-03);
	model.SetThermalConductivity(0.0242);
	model.SetAngle(0.000001);


	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);

	const double xPlateStart = 0.2;

	std::vector<Face*> faces = grid.faces.getLocalNodes();
	for (int i = 0; i<faces.size(); i++) {
			Face& f = *faces[i];
			if((f.BCMarker==3)&&(f.FaceCenter.x < xPlateStart)){
				f.BCMarker = 5;
			};
	};

	////Bind computational grid
	model.BindGrid(grid);	

	////Initial conditions
	ConservativeVariables initValues(0);

	Vector velocity(10.0,0,0);
	double pressure = 101579;
	double temperature = 300.214;

	initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);
	model.SetInitialConditions(initValues);	

	//Boundary conditions
	//Inlet boundary
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", InletBC);
	model.SetBoundaryCondition("right", OutletBC);
	model.SetBoundaryCondition("bottom", NoSlipBC);
	model.SetBoundaryCondition("top", SymmetryBC);
	model.SetBoundaryCondition("bottom_left", SymmetryBC);

	//Set wall boundaries		
	model.SetWallBoundary("bottom", true);
	model.ComputeWallDistances();
	model.DistanceSorting();
	model.EnableViscous();
	//model.DisableViscous();
	model.SchemeOrder = 2;
	model.DisableLimiter();

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "solution";
	//model.LoadSolution(outputSolutionFile+".txt");

	//Run simulation
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
			model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
			model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
		};
		if (model.totalTime > 10000) break;
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
		model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
		model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
	};
}

void RunIncompressibleBlasius()
{
	Flow2D model;
	int Nx = 30;
	int My = 40;
	double Lx = 1.2;
	double Ly = 0.5;
	//Grid grid = GenGrid2D(Nx, My, Lx, Ly, 1.02, 1.05);
	Grid grid = GenGrid2D(Nx, My, Lx, Ly, 1.0, 1.0);
	model.BindGrid(grid);
	model.SetGridParameters(Nx, My, Lx, Ly);

	//set medium and flow properties
	model.medium.Set_Cv(1006.43 / 1.4);
	model.medium.Set_density(1.17);
	model.medium.Set_Gamma(1.4);
	model.medium.Set_Temperature(300.0);
	model.medium.Set_Uinf(10.0);
	model.medium.Set_viscosity(1.78e-3);

	//set initial conditions
	model.ComputeKineticViscosity();
	model.SetInitialCondition(0);
	model.SaveToTechPlot("2dFlowInit.dat");
	model.CreateConvectivityList();

	//model.LoadSolution("Flow2D.txt");

	//model.ComputeVfromU();
	//model.ComputeUfromV();
	model.ComputeBlasius();
	model.SaveSolution("Flow2D.txt");
	Flow2DResWrite("Flow2Dres.dat");
	model.SaveToTechPlot("2dFlow.dat");
	std::getchar();
};

void CellGradientTest(){

	Grid grid = Load2DTriangleGrid("D:\\Projects\\NewKidzTurbo\\Grids\\Bump Triangles\\bump.grid");

	//fill in BC Markers
	std::vector<Face*> faces = grid.faces.getLocalNodes();
	for(int i=0; i<grid.faces.size(); i++){
		Face& f = *faces[i];
		if(f.isExternal!=1) continue;
		//top border
		if(f.FaceCenter.y==2.0){
			f.BCMarker = 4;
			continue;
		};
		//left border
		if(f.FaceCenter.x==-2.0){
			f.BCMarker = 1;
			continue;
		};
		//right border
		if(f.FaceCenter.x>2.999){
			f.BCMarker = 2;
			continue;
		};
		//bottom border
		f.BCMarker = 3;
	};
	//Add Patches
	grid.addPatch("left", 1);
	grid.addPatch("right", 2);
	grid.addPatch("bottom", 3);
	grid.addPatch("top", 4);
	grid.ConstructAndCheckPatches();

	//create and set model
	Model<Roe3DSolverPerfectGas> model;	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	
	model.SetViscosity(0);
	model.SetThermalConductivity(0.0242);
	model.SetAngle(0.000001);
	
	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);
	model.BindGrid(grid);

	////Initial conditions
	ConservativeVariables initValues(0);
	double density = 1.0;
	std::vector<Cell*> cells = grid.cells.getLocalNodes();
	for(int i=0; i<cells.size(); i++)
	{
		Vector Point = cells[i]->CellCenter;
		int I = cells[i]->GlobalIndex;
		
		initValues.GlobalIndex = I;
		initValues.ro = 1.0;
		initValues.rou = Point.x;
		initValues.rov = Point.y;
		initValues.roE = Point.x + 2*Point.y;

		model.U[i] = initValues;
	};

	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", NoSlipBC);
	model.SetBoundaryCondition("right", NoSlipBC);
	model.SetBoundaryCondition("bottom", NoSlipBC);
	model.SetBoundaryCondition("top", NoSlipBC);

	model.ComputeConservativeVariablesGradients();
	model.SaveToTechPlot("grad_test_fields.dat");
	model.SaveCellGradientsToTechPlot2D("gradients.dat");

	return;
};

void RunBumpFlow(){

	Grid grid = Load2DTriangleGrid("D:\\Projects\\NewKidzTurbo\\Grids\\Bump Triangles\\bump.grid");
	//Grid grid = Load2DTriangleGrid("D:\\Bump\\bump.grid");	//for ICAD comp
	//check_grid(grid);

	//fill in BC Markers
	std::vector<Face*> faces = grid.faces.getLocalNodes();
	for(int i=0; i<grid.faces.size(); i++){
		Face& f = *faces[i];
		if(f.isExternal!=1) continue;
		//top border
		if(f.FaceCenter.y==2.0){
			f.BCMarker = 4;
			continue;
		};
		//left border
		if(f.FaceCenter.x==-2.0){
			f.BCMarker = 1;
			continue;
		};
		//right border
		if(f.FaceCenter.x>2.999){
			f.BCMarker = 2;
			continue;
		};
		//bottom border
		f.BCMarker = 3;
	};
	//Add Patches
	grid.addPatch("left", 1);
	grid.addPatch("right", 2);
	grid.addPatch("bottom", 3);
	grid.addPatch("top", 4);
	grid.ConstructAndCheckPatches();

	//create and set model
	Model<Roe3DSolverPerfectGas> model;	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	
	model.SetViscosity(0);
	model.SetThermalConductivity(0.0242);
	model.SetAngle(0.000001);
	
	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);
	model.BindGrid(grid);

	////Initial conditions
	ConservativeVariables initValues(0);

	//compute flow prarameters and dimensional values
	double Mach_inf = 0.3;
	double temperature = 300.214;
	double pressure = 101579;
	double sound_speed = sqrt(model.medium.Gamma*(model.medium.Gamma - 1.0)*model.medium.Cv*temperature);
	double u_inf = Mach_inf*sound_speed;
	double ro_inf = model.medium.Gamma*pressure/(sound_speed*sound_speed);
	double pressure_gamma = pressure*model.medium.Gamma;
	Vector velocity(u_inf,0,0);
	ConservativeVariables DimVal(0);
	DimVal.ro = ro_inf;
	DimVal.rou = ro_inf*sound_speed;
	DimVal.rov = ro_inf*sound_speed;
	DimVal.row = ro_inf*sound_speed;
	DimVal.roE = ro_inf*sound_speed*sound_speed;

	initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);
	model.SetInitialConditions(initValues);	

	//Boundary conditions
	//Inlet boundary
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	/*model.SetBoundaryCondition("left", InletBC);
	model.SetBoundaryCondition("right", OutletBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);*/
	model.SetBoundaryCondition("left", InletBC);
	model.SetBoundaryCondition("right", InletBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", InletBC);

	//Set wall boundaries		
	model.SetWallBoundary("bottom", true);
	model.ComputeWallDistances();
	model.DistanceSorting();
	model.DisableViscous();
	model.SchemeOrder = 2;

	//Load solution
	model.LoadSolution("solution.txt");
	std::string outputSolutionFile = "solution";

	//Run simulation
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlotUndim(outputSolutionFile+".dat", DimVal);
		};
		if (model.totalTime > 10000) break;
	};
	std::getchar();

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlotUndim(outputSolutionFile+".dat", DimVal);
	};

	return;
};

void RunSteadyShock(){

	double xl = -1.0;
	double xr = 1.0;
	int N_cells = 81;	//in Katate code N_nodes is 81
	double h = (xr - xl)/(N_cells - 1);

	//real size
	xl -= 0.5*h;
	xr += 0.5*h;

	//generate the same grid (3D)
	Grid grid = GenGrid1D(N_cells, xl, xr, Vector(1,0,0));

	//read solution file
	double read;	//variable for reading
	std::vector<double> den, vel, press, temp;
	std::ifstream ifs("D:\\Projects\\NewKidzTurbo\\Grids\\SteadyShock1D\\ns_shock_structure.dat");
	ifs >> read;
	while(!ifs.eof())
	{
		ifs >> read;
		ifs >> read;
		den.push_back(read);
		ifs >> read;
		vel.push_back(read);
		ifs >> read;
		press.push_back(read);
		ifs >> read;
		temp.push_back(read);
		ifs >> read >> read;
	};
	//create and set model
	Model<Roe3DSolverPerfectGas> model;	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	
	model.SetViscosity(1.0);
	model.SetThermalConductivity(0.0242);
	model.SetAngle(0.000001);
	
	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);
	model.BindGrid(grid);

	////Initial conditions
	double rho_inf = 1.17;		//my value
	double R = (model.medium.Gamma - 1.0)*model.medium.Cv; // R/M for air
	double T_inf = 400.0;
	double c_inf = sqrt(model.medium.Gamma*R*T_inf); //speed of sound for air while 400K
	double p_inf = rho_inf*(model.medium.Gamma - 1.0)*model.medium.Cv*T_inf;

	//set initial condition according to Katate solution
	std::vector<Cell*> cells = grid.cells.getLocalNodes();
	for(int i=0; i<cells.size(); i++)
	{
		Cell c = *cells[i];
		ConservativeVariables VAR;
		VAR.GlobalIndex = c.GlobalIndex;
		VAR.ro = den[i]*rho_inf;
		VAR.rou = VAR.ro*c_inf*vel[i];
		VAR.rov = 0;
		VAR.row = 0;
		double pressure = press[i]*p_inf;
		double roe = pressure/(model.medium.Gamma - 1.0);
		VAR.roE = roe + 0.5*VAR.rou*VAR.rou/VAR.ro;
		model.U[c.GlobalIndex] = VAR;
	};

	//save our left and right values
	double PRESSURE_INLET_l = press[0]*p_inf;
	double TEMPERATURE_INLET_l = PRESSURE_INLET_l/((model.medium.Gamma - 1.0)*den[0]*rho_inf*model.medium.Cv);
	double VELOCITY_INLET_l = c_inf*vel[0];

	double PRESSURE_INLET_r = press[N_cells - 1]*p_inf;
	double TEMPERATURE_INLET_r = PRESSURE_INLET_r/((model.medium.Gamma - 1.0)*den[N_cells - 1]*rho_inf*model.medium.Cv);
	double VELOCITY_INLET_r = c_inf*vel[N_cells - 1];

	//Boundary conditions
	//Inlet boundary
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(PRESSURE_INLET_l, TEMPERATURE_INLET_l, Vector(VELOCITY_INLET_l, 0, 0));

	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition OutletBC(model);
	OutletBC.setParams(PRESSURE_INLET_r, TEMPERATURE_INLET_r, Vector(VELOCITY_INLET_r, 0, 0));

	//Set boundary conditions
	model.SetBoundaryCondition("left", InletBC);
	model.SetBoundaryCondition("right", OutletBC);

	//Set wall boundaries		
	model.EnableViscous();
	model.SchemeOrder = 1;

	//Load solution
	//model.LoadSolution("solution.txt");
	model.SaveToTechPlot("init.dat");
	std::string outputSolutionFile = "solution";

	//Run simulation
	bool isSave = true;	
	for (int i = 0; i < 1901; i++) {
		model.Step();	
		if (i % 50 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 1 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			//model.SaveToTechPlotUndim(outputSolutionFile+".dat", DimVal);
			model.SaveToTechPlot(outputSolutionFile+".dat");
		};
		if (model.totalTime > 10000) break;
	};
	std::getchar();

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		//model.SaveToTechPlotUndim(outputSolutionFile+".dat", DimVal);
	};

	return;
};

void RunVoronka()
{
	double alpha = 4*PI/180;
	std::string solutionFile = "D:\\Projects\\NewKidzTurbo\\Solutions\\Voronka.cgns";
	Grid grid = LoadCGNSGrid(solutionFile);
	grid.addPatch("side_wall", 6);
	grid.ComputeGeometryToAxisymmetric(alpha, 6);

	// Initialize medium model and place boundary conditions
	Model<Roe3DSolverPerfectGas> model;

	//Set fluid properties
	//Air
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);

	//model.SetViscosity(1.7894e-05);
	model.SetViscosity(1.7894e-03);
	model.SetThermalConductivity(0.0242);
	model.SetAngle(0.000001);
	model.SetAxiSymmetric();
	
	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);

	////Bind computational grid
	model.BindGrid(grid);

	////Initial conditions
	ConservativeVariables initValues(0);

	Vector velocity(-10.0,0,0);
	double pressure = 101579;
	double temperature = 300;

	//Set initial conditions
	initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);
	model.SetInitialConditions(initValues);

	//Read solution from CGNS
	model.ReadSolutionFromCGNS(solutionFile);

	//Boundary conditions
	//Inlet boundary
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Periodic boundary condition
	Model<Roe3DSolverPerfectGas>::PeriodicAxisymmetricBoundaryCondition PeriodicBC(model);
	PeriodicBC.setParams(alpha);

	//Set boundary conditions
	model.SetBoundaryCondition("inlet", InletBC);
	model.SetBoundaryCondition("outlet", OutletBC);
	model.SetBoundaryCondition("wall", NoSlipBC);
	model.SetBoundaryCondition("front_wall", NoSlipBC);
	model.SetBoundaryCondition("symmetry", SymmetryBC);
	model.SetBoundaryCondition("side_wall", PeriodicBC);

	//Set wall boundaries		
	//model.SetWallBoundary("wall", true);
	//model.SetWallBoundary("front_wall", true);
	//model.ComputeWallDistances();
	//model.DistanceSorting();
	model.DisableViscous();
	model.SchemeOrder = 2;	

	//Save initial solution
	//model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "solution";
	//model.LoadSolution("solution.txt");
	//return 0;

	//Run simulation
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 10 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
		};
		if (model.totalTime > 10000) break;
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");;
	};
};

void ConvertGrid(char* file_name)
{
	std::string gridFile = "D:\\Projects\\NewKidzTurbo\\Grids\\SimpleCircle.cgns";
	Grid grid = LoadCGNSGrid(gridFile);

	//write full info about grid
	std::ofstream ofs(file_name);
	std::vector<Node*> nodes = grid.nodes.getLocalNodes();
	std::vector<Cell*> cells = grid.cells.getLocalNodes();
	std::vector<Face*> faces = grid.faces.getLocalNodes();

	//write nodes info
	ofs << nodes.size() << '\n';
	for(int i=0; i<nodes.size(); i++)
	{
		ofs << nodes[i]->GlobalIndex << '\n';
		ofs << nodes[i]->P.x << ' ';
		ofs << nodes[i]->P.y << ' ';
		ofs << nodes[i]->P.z << '\n';
	};

	//write cells info
	ofs << cells.size() << '\n';
	for(int i=0; i<cells.size(); i++)
	{
		ofs << cells[i]->GlobalIndex << '\n';
		ofs << cells[i]->CellCenter.x << ' ';
		ofs << cells[i]->CellCenter.y << ' ';
		ofs << cells[i]->CellCenter.z << '\n';
		ofs << cells[i]->CellVolume << '\n';
		int faces_size = cells[i]->Faces.size();
		ofs << faces_size;
		for(int j=0; j<faces_size; j++) ofs << ' ' << cells[i]->Faces[j];
		ofs << '\n';
		int nodes_size = cells[i]->Nodes.size();
		ofs << nodes_size;
		for(int j=0; j<nodes_size; j++) ofs << ' ' << cells[i]->Nodes[j];
		ofs << '\n';
	};

	//write faces info
	ofs << faces.size() << '\n';
	for(int i=0; i<faces.size(); i++)
	{
		ofs << faces[i]->GlobalIndex << '\n';
		ofs << faces[i]->FaceCenter.x << ' ';
		ofs << faces[i]->FaceCenter.y << ' ';
		ofs << faces[i]->FaceCenter.z << '\n';
		ofs << faces[i]->isExternal << '\n';
		ofs << faces[i]->BCMarker << '\n';
		ofs << faces[i]->FaceSquare << '\n';
		ofs << faces[i]->FaceNormal.x << ' ';
		ofs << faces[i]->FaceNormal.y << ' ';
		ofs << faces[i]->FaceNormal.z << '\n';
		ofs << faces[i]->FaceCell_1 << ' ';
		ofs << faces[i]->FaceCell_2 << '\n';
		int nodes_size = faces[i]->FaceNodes.size();
		ofs << nodes_size;
		for(int j=0; j<nodes_size; j++) ofs << ' ' << faces[i]->FaceNodes[j];
		ofs << '\n';
	};
	ofs.close();

	return;
};

//Main program ))
int main(int argc, char *argv[]) {
	//RunSAFlatPlate();
	//RunGAWCalculation();
	//RunPoiseuilleTest();
	//RunRiemannProblemTest(ToroTestInitDistribution1, 100, 0.15, 2);
	//RunRiemannProblemTest(SODTestInitDistribution, 400, 0.2, 2);
	//RunShearFlowTest();
	//RunBlasiusTest();		//test has bad grid (not from ANSYS)
	//RunIncompressibleBlasius();
	//RunBumpFlow();
	//CellGradientTest();
	//RunSteadyShock();
	//RunVoronka();
	//ConvertGrid("SimpleCircle.dat");
	//return 0;

	//Load cgns grid			
	//std::string solutionFile = "D:\\Projects\\NewKidzTurbo\\Solutions\\Laminar_70ms_Air.cgns";
	//std::string solutionFile = "D:\\BlasLimiters\\solution.cgns";
	std::string solutionFile = "D:\\Projects\\NewKidzTurbo\\Solutions\\solution.cgns";
	Grid grid = LoadCGNSGrid(solutionFile);
	//check_grid(grid);

	// Initialize medium model and place boundary conditions
	Model<Roe3DSolverPerfectGas> model;

	//Set fluid properties
	//Air
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);

	//model.SetViscosity(1.7894e-05);
	model.SetViscosity(1.7894e-03);
	model.SetThermalConductivity(0.0242);
	model.SetAngle(0.000001);
	
	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.000);

	////Bind computational grid
	model.BindGrid(grid);	

	////Initial conditions
	ConservativeVariables initValues(0);

	Vector velocity(10.0,0,0);
	//Vector velocity(70,0,0);
	double pressure = 101579;
	double temperature = 300.214;

	initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);
	model.SetInitialConditions(initValues);	
	//Determine plate start coordinate
	const double xPlateStart = 0.2;
	//model.SetInitialConditionsBlasius(xPlateStart, velocity.x ,initValues.ro, initValues.roE);
	//model.SetInitialConditionsLinear(initValues);
	//model.SetInitialParabola();
	//model.ComputeParabolaGradients();	

	//Read solution from CGNS
	model.ReadSolutionFromCGNS(solutionFile);	
	//model.SaveSliceToTechPlot("uVisc.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
	//model.ReadSolutionFromCGNS("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\solutionNoVisc.cgns");
	//model.SaveSliceToTechPlot("uNoVisc.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
	//model.ReadSolutionFromCGNS("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\solutionInit.cgns");
	//model.Init();	
	//return 0;

	//Boundary conditions
	//Inlet boundary
	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
	InletBC.setParams(pressure, temperature, velocity);

	//Outlet boundary
	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
	OutletBC.setParams(pressure);

	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("inlet", InletBC);
	model.SetBoundaryCondition("outlet", OutletBC);
	model.SetBoundaryCondition("plate", NoSlipBC);
	//model.SetBoundaryCondition("plate", SymmetryBC);
	model.SetBoundaryCondition("symmetry", SymmetryBC);
	model.SetBoundaryCondition("top_left", SymmetryBC);
	model.SetBoundaryCondition("top_right", SymmetryBC);

	//Set wall boundaries		
	model.SetWallBoundary("plate", true);
	model.ComputeWallDistances();
	model.DistanceSorting();
	model.EnableViscous();
	model.SchemeOrder = 2;
	model.EnableLimiter();
	//model.DisableLimiter();
	//model.DisableViscous();
	
	//Init model
	//model.Init();		

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//model.SaveSliceToTechPlot("u0_8.dat", 0.2, 72.0, 0.76, 0.8, 0, 0.06);
	//model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
	//model.ComputeBoundaryLayerHeight(0.0001);
	//model.SaveBoundaryLayerHeightToTechPlot("BL_height.dat");

	//Load solution
	std::string outputSolutionFile = "solution";
	//model.LoadSolution("solution.txt");
	//return 0;

	//Run simulation
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
			model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
			model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
		};
		if (model.totalTime > 10000) break;
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
		model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
		model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
	};
}
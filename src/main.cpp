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

bool RunSODTest() {
	Model<Roe3DSolverPerfectGas> model;
	Grid grid;

	//Shock tube problem setting
	double lBegin = 0;
	double lEnd = 1.0;
	Vector direction = Vector(1,0,0);
	grid = GenGrid1D(1000, lBegin, lEnd, direction);
	//grid = GenGrid2D(100, 10, 1.0, 1.0, 1.0, 1.0);
	
	//Set fluid properties	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.DisableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.005);

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(SODTestInitDistribution);

	//Boundary conditions
	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
	model.SetBoundaryCondition("left", NoSlipBC);
	model.SetBoundaryCondition("right", NoSlipBC);

	/*Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	model.SetBoundaryCondition("top", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);*/

	model.SaveToTechPlot("SODInit.dat");

	//Total time
	double maxTime = 0.2;
	for (int i = 0; i < 200000000; i++) {
		model.Step();
		if (i % 1 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		model.SaveToTechPlot("SOD.dat");
		if (model.totalTime > 0.2) break;
	};

	//Save result to techplot
	model.SaveToTechPlot("SOD.dat");
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

void RunBlaiusTest(){	

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
		if (i % 1 == 0) {
			std::cout<<"Interation = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 1 == 0) && (isSave)) {
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

void RunBumpFlow(){

	Grid grid = Load2DTriangleGrid("D:\\Projects\\NewKidzTurbo\\Grids\\Bump Triangles\\bump.grid");
	//Grid grid = Load2DTriangleGrid("D:\\Projects\\NewKidzTurbo\\Grids\\Bump Triangles\\Test.grid");
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
	model.SetBoundaryCondition("right", OutletBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Set wall boundaries		
	model.SetWallBoundary("bottom", true);
	model.ComputeWallDistances();
	model.DistanceSorting();
	model.DisableViscous();

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

//Main program ))
int main(int argc, char *argv[]) {
	//RunSAFlatPlate();
	//RunGAWCalculation();
	//RunPoiseuilleTest();
	//RunSODTest();
	//RunShearFlowTest();
	//RunBlaiusTest();
	//RunIncompressibleBlasius();
	RunBumpFlow();
	return 0;

	//Load cgns grid			
	//std::string solutionFile = "D:\\Projects\\NewKidzTurbo\\Solutions\\Laminar_70ms_Air.cgns";
	std::string solutionFile = "D:\\Projects\\NewKidzTurbo\\Solutions\\solutionUG.cgns";	
	Grid grid = LoadCGNSGrid(solutionFile);
	check_grid(grid);

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
	//model.DisableViscous();
	
	//Init model
	//model.Init();		

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");
	std::vector<Cell*> cells = grid.cells.getLocalNodes();
	int Lb = 0;
	for(int i=0; i<cells.size(); i++)
	{
		Cell c_b = *cells[i];
		for(int j=0; j<4; j++)
		{
			if((c_b.Faces[j]==31)&&(c_b.GlobalIndex!=77))  Lb = c_b.GlobalIndex;
		};
	};
	Cell c = grid.cells[77];
	Face f_up = grid.faces[c.Faces[0]];
	Face f_left = grid.faces[c.Faces[1]];
	Face f_down = grid.faces[c.Faces[2]];
	Face f_right = grid.faces[c.Faces[3]];

	Cell cb = grid.cells[10];
	Face fb_right = grid.faces[cb.Faces[0]];
	Face fb_up = grid.faces[cb.Faces[1]];
	Face fb_left = grid.faces[cb.Faces[2]];
	Face fb_down = grid.faces[cb.Faces[3]];

	//model.SaveSliceToTechPlot("u0_8.dat", 0.2, 72.0, 0.76, 0.8, 0, 0.06);
	//model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
	//model.ComputeBoundaryLayerHeight(0.0001);
	//model.SaveBoundaryLayerHeightToTechPlot("BL_height.dat");

	//Load solution
	std::string outputSolutionFile = "solution";
	//model.LoadSolution("solution.txt");
	//model.LoadSolution(solutionFile+".txt");
	//model.LoadSolution("sol.txt");
	//model.LoadSolution("solInviscid.txt");
	//model.LoadSolution("solViscousBase.txt");
	//model.LoadSolution("solViscousTest2.txt");
	//model.SaveSliceToTechPlot("u.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
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
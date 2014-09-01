#include "stdafx.h"

#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include "GridLoading\cgnsload.h"
#include "grid.h"
#include "gridgeneration.h"
#include "model.h"
#include "riemannsolvers.h"
#include "garbarukmodel.h"
#include "SAmodel.h"
#include "biffplatemodel.h"
#include "bicgstab.h"
#include "gridloading\triangleload.h"
#include "mpi.h"
#include "kernel.h"


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
//void ComputeTVisc(char* file_name)
//{
//	function VEL;
//	VEL.ReadData1(file_name);
//	double K = 0.41;
//	double A = 12.0;
//	double rho = 1.12;
//	double vis_d = 1.7894e-03;
//	double vis_k = vis_d/rho;
//	double tau_w = vis_d*VEL.value[1]/VEL.x[1];
//	double U_fr = sqrt(tau_w);
//	double h = VEL.x[VEL.size-1];
//	double U_inf = VEL.value[VEL.size-1];
//	function DT = (1.0 - VEL/U_inf);
//	double delta = (1.0 - VEL/U_inf).TrapeziumIntegral();
//
//	function t_vis = VEL;
//	for(int i=0; i<t_vis.size; i++)
//	{
//		double D = (1.0 - exp((-1.0)*U_fr*t_vis.x[i]/(vis_k*A)));
//		D = D*D*D;
//		double gamma = t_vis.x[i]/h;
//		gamma = 1.0 + 5.5*gamma*gamma*gamma*gamma*gamma*gamma;
//		gamma = 1/gamma;
//		t_vis.value[i] = K*U_fr*min(t_vis.x[i]*D, delta*gamma);
//	};
//	t_vis.WriteData("t_vis_new.dat");
//};
//void TestSolver()
//{
//	BiffTurbulentPlateModel<Roe3DSolverPerfectGas> model(1375.0, 7562500.0, 2.0);
//	model.SetViscosity(1.7894e-03);
//
//	function vel;
//	BM_FUN_t_visc.ReadData1("t_visc.dat");
//	vel.ReadData1("t_velocity.dat");
//	BM_FUN_t_dif.ReadData1("t_visc.dat");
//	BM_FUN_lh_vel_der.ReadData1("t_LeeHarsha.dat");
//
//	BM_FUN_t_dif_der = BM_FUN_t_dif.Derivate();
//	BM_FUN_vel_der = vel.Derivate();
//	BM_FUN_lh_vel_der = BM_FUN_lh_vel_der.Derivate();
//	
//	/*int N = 100;
//	std::vector<double> subgrid(N+1);
//	double h = vel.x[vel.size-1];
//	for(int i=0; i<subgrid.size(); i++) subgrid[i] = i*h/N;
//	BM_FUN_t_visc = BM_FUN_t_visc.BindToNewGridFO(subgrid);
//	BM_FUN_vel_der = BM_FUN_vel_der.BindToNewGridFO(subgrid);
//	BM_FUN_t_dif = BM_FUN_t_dif.BindToNewGridFO(subgrid);
//	BM_FUN_t_dif_der = BM_FUN_t_dif_der.BindToNewGridFO(subgrid);
//	BM_FUN_lh_vel_der = BM_FUN_lh_vel_der.BindToNewGridFO(subgrid);*/
//
//	function res = model.SolverTest(vel);
//	res.WriteData("tau_un_test.dat");
//	std::getchar();
//};
//
//
//ConservativeVariables SODTestInitDistribution(Vector r, void *par) {	
//	//Gamma
//	double gamma = 1.4;
//
//	//Left state x <= 0.5
//	double roL = 1.0;
//	double pL = 1.0;
//	Vector vL = Vector(0,0,0);
//	ConservativeVariables UL;
//	UL.ro = roL;
//	UL.rou = roL * vL.x;
//	UL.rov = roL * vL.y;
//	UL.row = roL * vL.z;
//	UL.roE = pL/(gamma-1) + roL * vL.mod() * vL.mod() / 2.0;
//	//Right state x > 0.5
//	double roR = 0.125;
//	double pR = 0.1;
//	Vector vR = Vector(0,0,0);
//	ConservativeVariables UR;
//	UR.ro = roR;
//	UR.rou = roR * vR.x;
//	UR.rov = roR * vR.y;
//	UR.row = roR * vR.z;
//	UR.roE = pR/(gamma-1) + roR * vR.mod() * vR.mod() / 2.0;
//
//	if (r.x <= 0.5) {
//		return UL;
//	} else {
//		return UR;
//	};
//};
//
//bool RunSODTest() {
//	//Model<Roe3DSolverPerfectGas> model;
//	Model<Godunov3DSolverPerfectGas> model;
//	Grid grid;
//
//	//Shock tube problem setting
//	double lBegin = 0;
//	double lEnd = 1.0;
//	Vector direction = Vector(1,0,0);
//	grid = GenGrid1D(1000, lBegin, lEnd, direction);
//	//grid = GenGrid2D(100, 10, 1.0, 1.0, 1.0, 1.0);
//	
//	//Set fluid properties	
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	model.DisableViscous();	
//
//	//Set computational settings
//	model.SetCFLNumber(0.5);
//	model.SetHartenEps(0.00);
//	model.DisableSecondOrder();
//
//	//Bind computational grid
//	model.BindGrid(grid);	
//
//	//Set initial conditions
//	model.SetInitialConditions(SODTestInitDistribution);
//
//	//Boundary conditions
//
//	//No slip boundary
//	//Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//	Model<Godunov3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//	model.SetBoundaryCondition("left", NoSlipBC);
//	model.SetBoundaryCondition("right", NoSlipBC);	
//
//	/*Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//	model.SetBoundaryCondition("top", SymmetryBC);
//	model.SetBoundaryCondition("bottom", SymmetryBC);*/
//
//	model.SaveToTechPlot("SODInit.dat");
//
//	//Total time
//	double maxTime = 0.2;
//	for (int i = 0; i < 2000000000; i++) {
//		model.ExplicitTimeStep();	
//		//model.Step();
//		if (i % 1 == 0) {
//			std::cout<<"Interation = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//		};
//		model.SaveToTechPlot("SOD.dat");
//		if (model.totalTime > 0.2) break;
//	};
//
//	//Save result to techplot
//	model.SaveToTechPlot("SOD.dat");
//	return true;
//};

//struct PoiseuilleSettingStruct {
//	double PIn;
//	double POut;
//	double YMin;
//	double YMax;
//	double XMin;
//	double XMax;
//	double Viscosity;
//	double Temperature;
//	double Gamma;
//	double InitDensity;
//};
//
//Vector PoiseuilleVelocityDistribution(Vector position, void *params) {
//	Vector velocity(0,0,0);
//	PoiseuilleSettingStruct info = *(PoiseuilleSettingStruct *)params;
//
//	//Velocity profile	
//	double L = info.XMax - info.XMin;
//	double deltaP = info.PIn - info.POut;	
//	double R = (info.YMax - info.YMin) / 2.0;
//	double YMid = (info.YMax + info.YMin) / 2.0;
//	double r = abs(position.y - YMid);
//	velocity.x = deltaP * ( R*R - r*r)/ ( 4 * info.Viscosity * L);
//
//	return velocity;
//};
//
//ConservativeVariables PoiseuilleTestInitDistribution(Vector r, void *params) {	
//	PoiseuilleSettingStruct info = *(PoiseuilleSettingStruct *)params;
//
//
//	//Compute conservative variables
//	double ro = 1.0;
//	double p = 1.0;
//	Vector v = PoiseuilleVelocityDistribution(r, params);
//
//	ConservativeVariables U;
//	U.ro = info.InitDensity;
//	U.rou = ro * v.x;
//	U.rov = ro * v.y;
//	U.row = ro * v.z;
//	U.roE = info.PIn/(info.Gamma-1.0) + ro * v.mod() * v.mod() / 2.0;
//	
//	return U;	
//};
//
//bool RunPoiseuilleTest() {		
//	Model<Roe3DSolverPerfectGas> model;
//	Grid grid;
//
//	//Poiseuille problem setting
//	PoiseuilleSettingStruct info;
//	double lBegin = 0;
//	double lEnd = 1.0;
//	info.PIn = 100020;
//	info.POut = 100000;
//	info.YMin = 0.0;
//	info.YMax = 1.0;
//	info.Viscosity = 1.0;
//	info.XMin = 0.0;
//	info.XMax = 1.0;		
//	info.Temperature = 300.0;
//	info.Gamma = 1.4;		
//	
//	grid = GenGrid2D(pHelper, 50, 50, info.XMax, info.YMax, 1.0, 1.0);
//	
//	//Set fluid properties	
//	model.SetGamma(info.Gamma);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	model.SetViscosity(info.Viscosity);
//	model.SetThermalConductivity(0.0);
//	model.EnableViscous();
//
//	//Set computational settings
//	model.SetCFLNumber(0.35);
//	model.SetHartenEps(0.00);
//
//	//Bind computational grid
//	model.BindGrid(grid);	
//
//	//Set initial conditions
//	ConservativeVariables inletValues = model.PrimitiveToConservativeVariables(Vector(0,0,0), info.PIn, info.Temperature, model.medium);
//	info.InitDensity = inletValues.ro;
//	model.SetInitialConditions(PoiseuilleTestInitDistribution, &info);
//
//	//Boundary conditions
//
//	//Inlet	
//	Model<Roe3DSolverPerfectGas>::SubsonicInletBoundaryCondition InletBC(model);	
//	InletBC.setParams(info.PIn, info.Temperature, Vector(0,0,0));
//	InletBC.setVelocityDistribution( PoiseuilleVelocityDistribution, &info);
//	model.SetBoundaryCondition("left", InletBC);
//
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);	
//	OutletBC.setParams(info.POut);
//	model.SetBoundaryCondition("right", OutletBC);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//	model.SetBoundaryCondition("top", NoSlipBC);
//	model.SetBoundaryCondition("bottom", NoSlipBC);	
//
//	//Total time
//	double maxTime = 0.2;
//	for (int i = 0; i < 200000000; i++) {
//		model.Step();	
//		if (i % 1 == 0) {
//			std::cout<<"Interation = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";			
//		};		
//		if ((i % 10 == 0) && (i != 0)) {
//			model.SaveToTechPlot("Poiseuille.dat");
//			model.SaveSolution("Poiseuille.txt" );
//		};
//		if (model.totalTime > 0.2) break;
//	};
//
//	//Save result to techplot
//	model.SaveToTechPlot("Poiseuille.dat");
//	
//	return true;
//};
//
//void RunGAWCalculation() {
//	//Load cgns grid	
//	std::string solFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Projects\\Grids\\GAW-1\\solution.cgns"; 
//	Grid grid = LoadCGNSGrid(solFile);
//
//	// Initialize medium model and place boundary conditions
//	Model<Roe3DSolverPerfectGas> model;
//
//	//Set fluid properties
//	//Air
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	model.SetViscosity(1.7894e-03);
//	model.SetThermalConductivity(0.0242);
//
//	//Set computational settings
//	model.SetCFLNumber(0.35);
//	model.SetHartenEps(0.000);
//
//	//Bind computational grid
//	model.BindGrid(grid);
//
//	//Set wall boundaries and compute wall distances		
//	model.SetWallBoundary("gaw1", true);
//	model.ComputeWallDistances();
//	model.DistanceSorting();
//
//	//Enable viscous part
//	model.EnableViscous();
//
//	//Load solution
//	model.ReadSolutionFromCGNS(solFile);
//
//	//Output solution to tecplot
//	model.SaveWallPressureDistribution("GAW1P.dat");
//	model.SaveToTecplot25D("GAW1.dat");
//
//	//Boundary conditions
//	//Inlet boundary
//	Vector velocity = 70 * Vector(0.9877,-0.1564,0); // 9 degrees
//	double pressure = 101579;
//	double temperature = 300.214;
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);	
//
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
//	OutletBC.setParams(pressure);
//
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	//Set boundary conditions
//	model.SetBoundaryCondition("inlet", InletBC);
//	model.SetBoundaryCondition("outlet", OutletBC);
//	model.SetBoundaryCondition("gaw1", NoSlipBC);	
//	model.SetBoundaryCondition("sym1", SymmetryBC);
//	model.SetBoundaryCondition("sym2", SymmetryBC);
//	model.SetBoundaryCondition("top", NoSlipBC);
//	model.SetBoundaryCondition("bottom", NoSlipBC);
//	return;
//};
//
//void RunSAFlatPlate() {
//	//Load cgns grid		
//	std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineSA.cgns";
//	Grid grid = LoadCGNSGrid(solutionFile);
//
//	//Initialize model and bind grid
//	SAModel<Roe3DSolverPerfectGas> model;
//	model.BindGrid(grid);
//	model.Init();
//
//	//Set fluid properties
//	//Air
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	model.SetViscosity(1.7894e-05);
//	model.SetThermalConductivity(0.0242);
//
//	//Set wall boundaries and compute wall distances		
//	model.SetWallBoundary("plate", true);
//	model.ComputeWallDistances();
//	model.DistanceSorting();
//
//	//Read initial solution from CGNS
//	model.ReadSolutionFromCGNS(solutionFile);		
//
//	//Set boundary conditions			
//	Vector velocity(70,0,0);
//	double pressure = 101579;
//	double temperature = 300.214;
//	ConservativeVariables inletValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);				
//
//	//Boundary conditions
//	//Inlet boundary
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);	
//
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
//	OutletBC.setParams(pressure);
//
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	//Set boundary conditions
//	model.SetBoundaryCondition("inlet", InletBC);
//	model.SetBoundaryCondition("outlet", OutletBC);
//	model.SetBoundaryCondition("plate", NoSlipBC);	
//	model.SetBoundaryCondition("symmetry", SymmetryBC);
//	model.SetBoundaryCondition("top_left", SymmetryBC);
//	model.SetBoundaryCondition("top_right", SymmetryBC);
//
//	//Inlet turbulence level
//	model.SetBoundarySAValue("inlet", 0);
//	model.SaveToTechPlot("solSA.dat");
//	//model.ComputeWallVariables();
//
//	//Run some SA steps with frozen flow field
//	//model.SAStep();
//
//	return;
//};
//
//void RunBiffFlatPlate() {
//	std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineSA.cgns";
//	//std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineLam.cgns";
//	//std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns";	
//	Grid grid = LoadCGNSGrid(solutionFile);
//
//	// Initialize medium model and place boundary conditions
//	BiffTurbulentPlateModel<Roe3DSolverPerfectGas> model(320000.0, 8337150864.0, 2.0);
//
//	//Set fluid properties
//	//Air
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//
//	model.SetViscosity(1.7894e-05);	
//	model.SetThermalConductivity(0.0242);
//
//	//Set computational settings
//	model.SetCFLNumber(0.35);
//	model.SetHartenEps(0.000);
//
//	////Bind computational grid
//	model.BindGrid(grid);			
//
//	////Initial conditions
//	//Read solution from CGNS
//	model.ReadSolutionFromCGNS(solutionFile);	
//	model.LoadSolution("C:\\Users\\Erik\\Dropbox\\Science\\!Projects\\solution.txt");
//
//	//Boundary conditions
//	//Inlet boundary
//	Vector velocity(70,0,0);
//	double pressure = 101579;
//	double temperature = 300.214;	
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);
//
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
//	OutletBC.setParams(pressure);
//
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	//Set boundary conditions
//	model.SetBoundaryCondition("inlet", InletBC);
//	model.SetBoundaryCondition("outlet", OutletBC);
//	model.SetBoundaryCondition("plate", NoSlipBC);	
//	model.SetBoundaryCondition("symmetry", SymmetryBC);
//	model.SetBoundaryCondition("top_left", SymmetryBC);
//	model.SetBoundaryCondition("top_right", SymmetryBC);
//
//	//Set wall boundaries		
//	model.SetWallBoundary("plate", true);
//	model.ComputeWallDistances();
//	model.DistanceSorting();
//	model.EnableViscous();	
//
//	//Save initial solution
//	model.ComputeGradients();
//	model.SaveToTechPlot("init2.dat");
//
//	//Output dimensionless profile for point around x = 1.0
//	model.ComputeWallVariables();
//	
//
//	std::string outputSolutionFile = "solution2";		
//
//	//Run simulation
//	bool isSave = true;	
//	for (int i = 0; i < 1000000; i++) {
//		model.Step();	
//		if (i % 10 == 0) {
//			std::cout<<"Interation = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//		};		
//		if ((i % 100 == 0) && (isSave)) {
//			model.SaveSolution(outputSolutionFile+".txt");
//			model.SaveToTechPlot(outputSolutionFile+".dat");
//			model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
//			model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
//		};
//		if (model.totalTime > 10000) break;
//	};
//
//	//Save result to techplot
//	if (isSave) {
//		model.SaveSolution(outputSolutionFile+".txt");
//		model.SaveToTechPlot(outputSolutionFile+".dat");
//		model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
//		model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);
//	};
//};
//
//void GodunovTests() {
//	Godunov3DSolverPerfectGas solver;
//	Face f;
//	f.FaceNormal = Vector(1,0,0);
//	ConservativeVariables UL;
//	ConservativeVariables UR;
//
//	////Numerical tests for exact solution of riemann proble
//	// See Toro p. 129
//	Godunov3DSolverPerfectGas::StarVariables starValues;
//	//Test 1 (SOD shock tube problem)
//	std::cout<<"Test 1 run (SOD shock tube)\n";
//	solver.SetGamma(1.4);
//	//Left state
//	double roL = 1.0;
//	double pL = 1.0;
//	double uL = 0;
//
//	//Right state
//	double roR = 0.125;
//	double pR = 0.1;
//	double uR = 0;
//	starValues = solver.ComputeStarVariables(roL, uL, pL, roR, uR, pR);
//
//	//Output results	
//	std::cout<<"pStar = "<<starValues.pStar<<"; pStarRight = "<<0.30313<<"\n";
//	std::cout<<"uStar = "<<starValues.uStar<<"; uStarRight = "<<0.92745<<"\n";
//	std::cout<<"roStarL = "<<starValues.roStarL<<"; roStarLRight = "<<0.42632<<"\n";
//	std::cout<<"roStarR = "<<starValues.roStarR<<"; roStarRRight = "<<0.26557<<"\n";
//	if (starValues.leftWave == Godunov3DSolverPerfectGas::Shock) {
//		std::cout<<"Left wave is shock.\n";
//	} else {
//		std::cout<<"Left wave is rarefaction.\n";
//	};
//	if (starValues.rightWave == Godunov3DSolverPerfectGas::Shock) {
//		std::cout<<"Right wave is shock.\n";
//	} else {
//		std::cout<<"Right wave is rarefaction.\n";
//	};
//	std::cout<<std::endl;
//	 
//	//Test 2 (123 problem)
//	std::cout<<"Test 2 run (123 problem)\n";
//	solver.SetGamma(1.4);
//	//Left state
//	roL = 1.0;
//	pL = 0.4;
//	uL = -2.0;
//
//	//Right state
//	roR = 1.0;
//	pR = 0.4;
//	uR = 2.0;
//	starValues = solver.ComputeStarVariables(roL, uL, pL, roR, uR, pR);
//
//	//Output results	
//	std::cout<<"pStar = "<<starValues.pStar<<"; pStarRight = "<<0.00189<<"\n";
//	std::cout<<"uStar = "<<starValues.uStar<<"; uStarRight = "<<0.00000<<"\n";
//	std::cout<<"roStarL = "<<starValues.roStarL<<"; roStarLRight = "<<0.02185<<"\n";
//	std::cout<<"roStarR = "<<starValues.roStarR<<"; roStarRRight = "<<0.02185<<"\n";	 
//	if (starValues.leftWave == Godunov3DSolverPerfectGas::Shock) {
//		std::cout<<"Left wave is shock.\n";
//	} else {
//		std::cout<<"Left wave is rarefaction.\n";
//	};
//	if (starValues.rightWave == Godunov3DSolverPerfectGas::Shock) {
//		std::cout<<"Right wave is shock.\n";
//	} else {
//		std::cout<<"Right wave is rarefaction.\n";
//	};
//	std::cout<<std::endl;
//
//	//Test 3
//	std::cout<<"Test 3 run\n";
//	solver.SetGamma(1.4);
//	//Left state
//	roL = 1.0;
//	pL = 1000.0;
//	uL = 0.0;
//
//	//Right state
//	roR = 1.0;
//	pR = 0.01;
//	uR = 0.0;
//	starValues = solver.ComputeStarVariables(roL, uL, pL, roR, uR, pR);
//
//	//Output results	
//	std::cout<<"pStar = "<<starValues.pStar<<"; pStarRight = "<<460.894 <<"\n";
//	std::cout<<"uStar = "<<starValues.uStar<<"; uStarRight = "<<19.5975 <<"\n";
//	std::cout<<"roStarL = "<<starValues.roStarL<<"; roStarLRight = "<<0.57506<<"\n";
//	std::cout<<"roStarR = "<<starValues.roStarR<<"; roStarRRight = "<<5.99924<<"\n";	 
//	if (starValues.leftWave == Godunov3DSolverPerfectGas::Shock) {
//		std::cout<<"Left wave is shock.\n";
//	} else {
//		std::cout<<"Left wave is rarefaction.\n";
//	};
//	if (starValues.rightWave == Godunov3DSolverPerfectGas::Shock) {
//		std::cout<<"Right wave is shock.\n";
//	} else {
//		std::cout<<"Right wave is rarefaction.\n";
//	};
//	std::cout<<std::endl;	     
//
//	return;
//};
//
//void LinearSolverTests() {
//	const int N = 2;
//	double* x = new double[N];
//	double* b = new double[N];
//	x[0] = 0.0;
//	x[1] = -0.1;
//	b[0] = 1;
//	b[1] = 4;	
//	SparseRowMatrix A(N);
//	A.setValue(0,0, 1.0);
//	A.setValue(0,1, 2.0);
//	A.setValue(1,0, 3.0);
//	A.setValue(1,1, 4.0);
//	bicgstab<SparseRowMatrix>(N, A, b, x, 1e-10, true);	
//	std::cout<<x[0]<<" "<<x[1]<<std::endl;
//};
//
//void ImplicitSolverTest() {
//	Model<Godunov3DSolverPerfectGas> model;
//	Grid grid;
//
//	//Shock tube problem setting
//	double lBegin = 0;
//	double lEnd = 1.0;
//	Vector direction = Vector(1,0,0);
//	grid = GenGrid1D(1000, lBegin, lEnd, direction);
//	//grid = GenGrid1D(2, lBegin, lEnd, direction);
//	
//	//Set fluid properties	
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	model.DisableViscous();	
//
//	//Bind computational grid
//	model.BindGrid(grid);	
//
//	//Set initial conditions
//	model.SetInitialConditions(SODTestInitDistribution);
//	/*Vector velocity(0,0,0);
//	double pressure = 101579;
//	double temperature = 300.214;
//
//	ConservativeVariables initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);		
//	model.SetInitialConditions(initValues);	*/
//
//	//Boundary conditions
//	//No slip boundary
//	Model<Godunov3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//	model.SetBoundaryCondition("left", NoSlipBC);
//	model.SetBoundaryCondition("right", NoSlipBC);
//
//	model.SaveToTechPlot("SODInit.dat");
//
//	//Total time
//	//double maxTime = 0.2;
//	//for (int i = 0; i < 2000000000; i++) {
//	//	model.ExplicitTimeStep();	
//	//	//model.Step();
//	//	if (i % 1 == 0) {
//	//		std::cout<<"Interation = "<<i<<"\n";
//	//		std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//	//		for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//	//		std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//	//	};
//	//	model.SaveToTechPlot("SOD.dat");
//	//	if (model.totalTime > 0.2) break;
//	//};
//
//	//Set computational settings
//	model.SetCFLNumber(1.0);
//	model.SetHartenEps(0.005);
//
//	//Run computation
//	model.ImplicitSteadyState(0.2, 2, 1e-15, 0.0);
//
//	//Save result to techplot
//	model.SaveToTechPlot("SOD.dat");
//	return ;
//};
//
//int BumpFlowTestExplicit() {
//	Grid grid = Load2DTriangleGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\BumpFlow\\bump.grid");
//
//	//fill in BC Markers
//	std::vector<Face*> faces = grid.faces.getLocalNodes();
//	for(int i=0; i<grid.faces.size(); i++){
//		Face& f = *faces[i];
//		if(f.isExternal!=1) continue;
//		//top border
//		if(f.FaceCenter.y==2.0){
//			f.BCMarker = 4;
//			continue;
//		};
//		//left border
//		if(f.FaceCenter.x==-2.0){
//			f.BCMarker = 1;
//			continue;
//		};
//		//right border
//		if(f.FaceCenter.x>2.999){
//			f.BCMarker = 2;
//			continue;
//		};
//		//bottom border
//		f.BCMarker = 3;
//	};
//	//Add Patches
//	grid.addPatch("left", 1);
//	grid.addPatch("right", 2);
//	grid.addPatch("bottom", 3);
//	grid.addPatch("top", 4);
//	grid.ConstructAndCheckPatches();
//
//	//create and set model
//	Model<Roe3DSolverPerfectGas> model;	
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	
//	model.SetViscosity(0);
//	model.SetThermalConductivity(0.0242);
//	model.SetAngle(0.000001);
//		
//	model.BindGrid(grid);
//
//	////Initial conditions
//	ConservativeVariables initValues(0);
//
//	double Mach_inf = 0.3;
//	double temperature = 300.214;
//	double pressure = 101579;
//	double sound_speed = sqrt(model.medium.Gamma*(model.medium.Gamma - 1.0)*model.medium.Cv*temperature);
//	double u_inf = Mach_inf*sound_speed;
//	double ro_inf = model.medium.Gamma*pressure/(sound_speed*sound_speed);
//	Vector velocity(u_inf,0,0);
//
//	//Set computational settings
//	model.SetCFLNumber(0.99);
//	model.SetHartenEps(0.000);
//	model.SetOperatingPressure(pressure);
//
//	initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);
//	model.SetInitialConditions(initValues);	
//
//	//Boundary conditions
//	//Inlet boundary
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);
//
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
//	OutletBC.setParams(pressure);
//
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	//Set boundary conditions
//	model.SetBoundaryCondition("left", InletBC);
//	model.SetBoundaryCondition("right", InletBC);
//	model.SetBoundaryCondition("bottom", SymmetryBC);
//	model.SetBoundaryCondition("top", InletBC);
//
//	//Set wall boundaries		
//	model.SetWallBoundary("bottom", true);
//	model.ComputeWallDistances();
//	model.DistanceSorting();
//	model.DisableViscous();
//	model.EnableSecondOrder();
//
//	//Save initial solution
//	model.SaveToTechPlot("init.dat");
//
//	//Load solution
//	std::string outputSolutionFile = "bumpRS";
//	model.LoadSolution("bump.txt");
//
//	//Compute flow prarameters and dimensional values	
//	ConservativeVariables DimVal(0);
//	DimVal.ro = ro_inf;
//	DimVal.rou = ro_inf*sound_speed;
//	DimVal.rov = ro_inf*sound_speed;
//	DimVal.row = ro_inf*sound_speed;
//	DimVal.roE = ro_inf*sound_speed*sound_speed;
//	//model.SaveToTechPlotUndim("bumpUD.dat", DimVal);
//	//return 0;
//
//	//Convergence history
//	std::ofstream fout("convRS.dat", std::ios::app | std::ios::out);
//
//	//Run simulation
//	bool isSave = true;	
//	for (int i = 0; i < 100000000; i++) {
//		model.ExplicitTimeStep();	
//		fout<<i<<" "<<model.stepInfo.Residual[4]<<"\n";
//		if (i % 10 == 0) {
//			std::cout<<"Iteration = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//		};
//		if ((i % 100 == 0) && (isSave)) {
//			model.SaveSolution(outputSolutionFile+".txt");
//			model.SaveToTechPlot(outputSolutionFile+".dat");
//			model.SaveToTechPlotUndim("bumpUD.dat", DimVal);
//		};
//		if (i > 15000) break;
//		if (model.totalTime > 13322) break;
//	};
//	//model.LoadSolution(outputSolutionFile+".txt");
//	//model.SetCFLNumber(1e10);
//	//model.ImplicitSteadyState(1e100, 1, 1e-15, 1.0);
//	//model.ImplicitSteadyState(1, 10, 1e-15, 1.0);
//
//	//Save result to techplot
//	if (isSave) {
//		model.SaveSolution(outputSolutionFile+".txt");
//		model.SaveToTechPlot(outputSolutionFile+".dat");
//	};
//
//	fout.close();
//
//	return 0;
//};
//
//int SimpleChannelTest() {
//	Grid grid = GenGrid2D(100, 100, 5.0, 2.0, 1.0, 1.0);
//
//	//create and set model
//	Model<Roe3DSolverPerfectGas> model;	
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	
//	model.SetViscosity(0);
//	model.SetThermalConductivity(0.0242);
//	model.SetAngle(0.000001);
//	model.DisableViscous();
//		
//	model.BindGrid(grid);
//		
//	////Initial conditions
//	ConservativeVariables initValues(0);
//
//	double Mach_inf = 0.3;
//	double temperature = 300.214;
//	double pressure = 101579;
//	double sound_speed = sqrt(model.medium.Gamma*(model.medium.Gamma - 1.0)*model.medium.Cv*temperature);
//	double u_inf = Mach_inf*sound_speed;
//	Vector velocity(u_inf,0,0);
//
//	//Set computational settings
//	model.SetCFLNumber(0.99);
//	model.SetHartenEps(0.000);
//	model.SetOperatingPressure(pressure);
//
//	//initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);
//	initValues = model.PrimitiveToConservativeVariables(Vector(0,0,0), pressure, temperature, model.medium);
//	model.SetInitialConditions(initValues);	
//
//	//Boundary conditions
//	//Inlet boundary
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);
//
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
//	OutletBC.setParams(pressure);
//
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	//Set boundary conditions
//	model.SetBoundaryCondition("left", InletBC);
//	model.SetBoundaryCondition("right", OutletBC);
//	model.SetBoundaryCondition("bottom", SymmetryBC);
//	model.SetBoundaryCondition("top", SymmetryBC);
//
//	//Set wall boundaries		
//	model.SetWallBoundary("bottom", true);
//	model.ComputeWallDistances();
//	model.DistanceSorting();
//	model.DisableViscous();
//
//	//Save initial solution
//	model.SaveToTechPlot("init.dat");
//
//	//Load solution
//	std::string outputSolutionFile = "simple";
//
//	//Run simulation
//	bool isSave = true;	
//	for (int i = 0; i < 100000000; i++) {
//		model.ExplicitTimeStep();	
//		if (i % 10 == 0) {
//			std::cout<<"Iteration = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//		};
//		if ((i % 100 == 0) && (isSave)) {
//			model.SaveSolution(outputSolutionFile+".txt");
//			model.SaveToTechPlot(outputSolutionFile+".dat");
//		};
//		if (i > 7500) break;
//		if (model.totalTime > 13322) break;
//	};
//
//	//Save result to techplot
//	if (isSave) {
//		model.SaveSolution(outputSolutionFile+".txt");
//		model.SaveToTechPlot(outputSolutionFile+".dat");
//	};
//
//	return 0;
//};
//
//int ImplicitSolverTestOneCell() {
//	Grid grid = GenGrid2D(1, 1, 1.0, 1.0, 1.0, 1.0);
//
//	// Initialize medium model and place boundary conditions
//	Model<Roe3DSolverPerfectGas> model;
//
//	//Set fluid properties
//	//Air
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);
//	model.DisableViscous();
//
//	model.SetViscosity(1.7894e-05);
//	//model.SetViscosity(1.7894e-03);
//	model.SetThermalConductivity(0.0242);
//	model.SetAngle(0.000001);
//
//	//Bind computational grid
//	model.BindGrid(grid);	
//
//	//Initial conditions
//	ConservativeVariables initValues(0);	
//	Vector velocity(10,0,0);
//	double pressure = 101579;
//	double temperature = 300.214;
//	model.SetOperatingPressure(pressure);
//
//	initValues = model.PrimitiveToConservativeVariables(velocity, pressure + 100, temperature, model.medium);		
//	model.SetInitialConditions(initValues);	
//
//	//Boundary conditions
//	//Inlet
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);	
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);	
//	OutletBC.setParams(pressure);
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	model.SetBoundaryCondition("left", InletBC);
//	model.SetBoundaryCondition("right", OutletBC);
//	model.SetBoundaryCondition("bottom", SymmetryBC);	
//	model.SetBoundaryCondition("top", SymmetryBC);
//
//	//Set computational settings
//	model.SetCFLNumber(1.0e100);
//	model.SetHartenEps(0.000);
//
//	//Run computation
//	bool isSave = true;
//	model.ImplicitSteadyState(1, 1, 1e-15, 0);
//	/*for (int i = 0; i < 1000000; i++) {
//		bool isConverged = true;
//		model.ExplicitTimeStep();	
//		if (i % 1 == 0)  {
//			std::cout<<"Interation = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) {
//				std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//				if (abs(model.stepInfo.Residual[k]) > 1e-10) isConverged = false;
//			};
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//		};
//		if (model.totalTime > 10000) break;
//		if (isConverged) break;
//	};*/
//
//	//Save result to techplot
//	model.SaveToTechPlot("solution.dat");
//	return 0;
//};
//
//int BlasiusTest() {
//	//Load cgns grid			
//	//std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\MeshUniform\\solution.cgns";
//	std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns";
//	//std::string solutionFile = "C:\\Users\\Erik\\Dropbox\\Science\\!Grants\\TSAGI\\FlatPlate\\FlatPlateFluent\\solutionSuperFineSA.cgns";	
//	Grid grid = LoadCGNSGrid(solutionFile);
//
//	// Initialize medium model and place boundary conditions
//	Model<Roe3DSolverPerfectGas> model;
//
//	//Set fluid properties
//	//Air
//	model.SetGamma(1.4);
//	model.SetCv(1006.43 / 1.4);
//	model.SetMolecularWeight(28.966);	
//
//	//model.SetViscosity(1.7894e-05);
//	model.SetViscosity(1.7894e-03);
//	model.SetThermalConductivity(0.0242);
//	model.SetAngle(0.000001);	
//
//	////Bind computational grid
//	model.BindGrid(grid);	
//
//	////Initial conditions
//	ConservativeVariables initValues(0);
//
//	//Vector velocity(0.1,0,0);
//	Vector velocity(10,0,0);
//	double pressure = 101579;
//	double temperature = 300.214;
//
//	//Set computational settings
//	model.SetCFLNumber(0.9);
//	model.SetHartenEps(0.00);
//	model.SetOperatingPressure(pressure);
//
//	/*initValues = model.PrimitiveToConservativeVariables(velocity, pressure, temperature, model.medium);		
//	model.SetInitialConditions(initValues);	*/
//	//Determine plate start coordinate
//	const double xPlateStart = 0.2;	
//
//	//Read solution from CGNS
//	model.ReadSolutionFromCGNS(solutionFile);		
//
//	//Boundary conditions
//	//Inlet boundary
//	Model<Roe3DSolverPerfectGas>::InletBoundaryCondition InletBC(model);
//	InletBC.setParams(pressure, temperature, velocity);
//
//	//Outlet boundary
//	Model<Roe3DSolverPerfectGas>::SubsonicOutletBoundaryCondition OutletBC(model);
//	//OutletBC.setParams(101604);
//	OutletBC.setParams(pressure);
//
//	//Symmetry boundary
//	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
//
//	//No slip boundary
//	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);
//
//	//Fixed velocity boundary condition
//	Model<Roe3DSolverPerfectGas>::ConstantVelocityBoundaryCondition FixedVelocityBC(model, velocity);
//
//	//Set boundary conditions
//	model.SetBoundaryCondition("inlet", InletBC);
//	model.SetBoundaryCondition("outlet", OutletBC);
//	model.SetBoundaryCondition("plate", NoSlipBC);	
//	model.SetBoundaryCondition("symmetry", SymmetryBC);
//	model.SetBoundaryCondition("top_left", OutletBC);
//	model.SetBoundaryCondition("top_right", OutletBC);
//	
//
//	//Set wall boundaries		
//	model.SetWallBoundary("plate", true);
//	model.ComputeWallDistances();
//	model.DistanceSorting();
//	model.EnableViscous();
//	//model.DisableViscous();
//	 
//	//Init model
//	//model.Init();		
//
//	//Save initial solution
//	model.ComputeGradients();
//	model.SaveToTechPlot("init.dat");
//	model.SaveSolution("init.sol");	
//
//	//Load solution
//	std::string outputSolutionFile = "solution";
//
//	////Run simulation
//	bool isSave = true;		
//	for (int i = 0; i < 1000000; i++) {
//		model.ExplicitTimeStep();			
//		if (i % 10 == 0)  {
//			std::cout<<"Interation = "<<i<<"\n";
//			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
//			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
//			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
//		};
//		if ((i % 100 == 0) && (isSave)) {
//			model.SaveSolution(outputSolutionFile+".sol");
//			model.SaveToTechPlot(outputSolutionFile+".dat");
//			/*model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
//			model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);*/
//			model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.992, 1.0, 0, 0.06);			
//		};
//		if (model.totalTime > 10000) break;
//	};
//
//	//model.LoadSolution(outputSolutionFile + ".sol");
//	//model.ImplicitSteadyState(10, 1, 1e-10, 1.0);
//
//	//Save result to techplot
//	if (isSave) {
//		model.SaveSolution(outputSolutionFile+".sol");
//		model.SaveToTechPlot(outputSolutionFile+".dat");
//		/*model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.96, 1.01, 0, 0.06);
//		model.SaveSliceToTechPlot("u0_8.dat", 0.2, 10.5, 0.76, 0.8, 0, 0.06);*/
//		model.SaveSliceToTechPlot("u1_0.dat", 0.2, 10.5, 0.992, 1.0, 0, 0.06);		
//	};
//	return 0;
//};

//SOD Test new version
class SodTestInitialConditions : public InitialConditions::InitialConditions
{
public:
	virtual std::vector<double> getInitialValues(const Cell& cell) {
		std::vector<double> initValues;
			//http://en.wikipedia.org/wiki/Sod_shock_tube
			//Grid sizes
			double Lx = 0.02; // SI 2 cm				

			//Other velocities
			double v = 0;
			double w = 0;

			//Left state
			double roL = 1.0;
			double PL = 1.0;
			double uL = 0;
			
			//Right state
			double roR = 0.125;
			double PR = 0.1;
			double uR = 0;

			//Cell center
			double x = cell.CellCenter.x;
			double y = cell.CellCenter.y;
			double z = cell.CellCenter.z;
				
			//Values
			double u = 0;
			double ro = 0;
			double roE = 0;
			if (x <= 0.5) {
				ro = roL;
				u = uL;
				roE = PL/(_gasModel->Gamma - 1.0);
			} else {
				ro = roR;
				u = uR;
				roE = PR/(_gasModel->Gamma - 1.0);
			};
			
			//Convert to conservative variables
			initValues.resize(_gasModel->nConservativeVariables);
			initValues[0] = ro;
			initValues[1] = ro * u;
			initValues[2] = ro * v;
			initValues[3] = ro * w;
			initValues[4] = roE;// + (u*u + v*v + w*w) / 2.0);

			return initValues;		
	};
};

void runSodTest(int argc, char *argv[]) {
	Kernel _kernel;	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 1000, 1, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0);
	_kernel.BindGrid(_grid);
	//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\Mixer\\Mixer.cgns");	
	//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns");	
	//_kernel.ReadInitialConditions("FlowSolution.E:1");
	_kernel.ReadConfiguration("");	
	_kernel.InitCalculation();

	//Initial conditions
	SodTestInitialConditions ic;
	_kernel.GenerateInitialConditions(ic);	

	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();	
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	//_kernel.SaveSolution();
	//_kernel.SaveGrid("grid.cgns");
	_kernel.Finalize();	
};

//Main program ))
int main(int argc, char *argv[]) {		
	//SimpleChannelTest(); return 0;
	//ImplicitSolverTestOneCell(); return 0;
	//BlasiusTest(); return 0;
	//BumpFlowTestExplicit(); return 0;
	//ImplicitSolverTest(); return 0;
	//LinearSolverTests(); return 0;
	//GodunovTests();
	//RunSAFlatPlate();
	//RunGAWCalculation();
	//RunPoiseuilleTest();
	//RunSODTest();
	//RunBiffFlatPlate();
    
	//RunSODTest();
	//return 0;

	runSodTest(argc, argv);
	return 0;
	
	Kernel _kernel;	
	_kernel.Initilize(&argc, &argv);
	Grid _grid = GenGrid2D(_kernel.getParallelHelper(), 30, 90, -0.25, 0.25, -0.75, 0.75, 1.0, 1.0);
	_kernel.BindGrid(_grid);
	//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\Mixer\\Mixer.cgns");	
	//_kernel.LoadGrid("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns");	
	//_kernel.ReadInitialConditions("FlowSolution.E:1");
	_kernel.ReadConfiguration("");
	_kernel.InitCalculation();
	_kernel.RunCalculation();
	_kernel.FinalizeCalculation();	
	_kernel.SaveGrid("result.cgns");
	_kernel.SaveSolution("result.cgns", "Solution");
	//_kernel.SaveSolution();
	//_kernel.SaveGrid("grid.cgns");
	_kernel.Finalize();	
	//_kernel.LoadGridTopologyAndInfo("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\SODShockTube\\SodShockTube.cgns");
	//std::cout<<argc<<"\n"<<argv[1]<<"\n";
	
	//_kernel.LoadGridTopology("C:\\Users\\Erik\\Dropbox\\Science\\ValidationCFD\\LaminarFlatPlate\\Mesh80\\solution.cgns");
	//_kernel.PartitionGrid();
	//_kernel.GenerateGridGeometry();
	

	return 0;	
};
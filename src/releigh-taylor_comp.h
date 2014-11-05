#ifndef TURBO_RELEIGH_TAYLOR_COMPUTATIONS
#define TURBO_RELEIGH_TAYLOR_COMPUTATIONS

#include "tests.h"
#include "RTI_solution_processor.h"

////Releigh-Taylor Instability RUN methods////

//2D initializing of init conditions
ConservativeVariables ReleighTaylorInit2D(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0.25*info.A*(1.0 + cos(2.0*PI*CellPos.x/(info.x_max - info.x_min)))*(1.0 + cos(2.0*PI*CellPos.y/(info.y_max - info.y_min)));

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//2D discrete perturbation
ConservativeVariables ReleighTaylorInit2D_discrete_perturbations(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<0.2*info.x_max)&&(abs(CellPos.y)<0.05*info.y_max)) v.y = info.A*info.x_max;
	
	//perturbation in one central cell for grid 31х91
	//if((abs(CellPos.x)<info.x_max/31.0)&&(abs(CellPos.y)<info.y_max/91.0)) v.y = info.A*info.x_max;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//2D discrete perturbation by velocity(постановка Фортовой, variant 1)
ConservativeVariables ReleighTaylorInit2D_var1(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<4.0e-4)&&(abs(CellPos.y)<4.0e-4)) v.y = info.A;
	
	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//2D discrete perturbation by velocity(постановка Фортовой, variant 2)
ConservativeVariables ReleighTaylorInit2D_var2(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<4.0e-4)&&(abs(CellPos.y)<4.0e-4)&&(CellPos.y>0)) ro = 0.1*info.ro_top;
	
	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};


//3D discrete perturbation
ConservativeVariables ReleighTaylorInit3D_discrete_perturbations(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<0.2*info.x_max)&&(abs(CellPos.z)<0.2*info.z_max)&&(abs(CellPos.y)<0.1*info.y_max)) v.y = info.A*info.x_max;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//3D smooth perturbation
ConservativeVariables ReleighTaylorInit3D(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0.125*info.A*(1.0 + cos(2.0*PI*CellPos.x/(info.x_max - info.x_min)))*(1.0 + cos(2.0*PI*CellPos.y/(info.y_max - info.y_min)))*(1.0 + cos(2.0*PI*CellPos.z/(info.z_max - info.z_min)));

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//3D discrete perturbation in z direction
ConservativeVariables ReleighTaylorInit3D_discrete_perturbations_z(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.z<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.z;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<0.2*info.x_max)&&(abs(CellPos.y)<0.2*info.y_max)&&(abs(CellPos.z)<0.1*info.z_max)) v.z = info.A*info.x_max;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//3D discrete perturbation in z direction (also x and y velocity components are not zero)
ConservativeVariables ReleighTaylorInit3D_discrete_perturbations_z_xy(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.z<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.z;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<0.15*info.x_max)&&(abs(CellPos.y)<0.15*info.y_max)&&(abs(CellPos.z)<0.1*info.z_max)) {
		v.z = info.A*info.x_max;
		v.x = v.z;
		v.y = v.z;
	};

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//3D discrete perturbation by velocity(постановка Фортовой, variant 1)
ConservativeVariables ReleighTaylorInit3D_var1(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<4.0e-4)&&(abs(CellPos.z)<4.0e-4)&&(abs(CellPos.y)<4.0e-4)) v.y = info.A;
	
	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
//3D discrete perturbation by velocity(постановка Фортовой, variant 2)
ConservativeVariables ReleighTaylorInit3D_var2(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if((abs(CellPos.x)<4.0e-4)&&(abs(CellPos.z)<4.0e-4)&&(abs(CellPos.y)<4.0e-4)&&(CellPos.y>0)) ro = 0.1*info.ro_top;
	
	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};



//2D classic test
void RunReleighTaylor2D(int N_x, int N_y)
{
	//set task parameters
	Rectangle_Area_RT info;

	//classic settings
	/*info.A = 0.01;	//perturbation amplitude
	info.g = 0.25;	//acceleration
	info.P0 = 2.5;	//zero altitude pressure
	info.ro_bot = 1.0;	//bottom density
	info.ro_top = 2.0;	//top pressure
	info.gamma = 1.4;	//constant of adiabat
	info.x_min = -0.25;
	info.x_max = (-1)*info.x_min;
	info.y_min = -0.75;
	info.y_max = (-1)*info.y_min;*/

	//Fe Pl sinusoidal perturbations
	/*info.A = 5.0e2;
	info.g = 5.0e8;
	info.gamma = 1.4;
	info.P0 = 2.0e12;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.03;
	info.y_max = 0.03;*/

	//Fe Pe discrete perturbation
	info.A = 5.0e5;
	info.g = -5.0e8;
	info.gamma = 1.4;
	info.P0 = 2.0e12;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.03;
	info.y_max = 0.03;

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.0);
	model.SchemeOrder = 2;
	model.EnableLimiter();
	//model.DisableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(ReleighTaylorInit2D_discrete_perturbations, &info);
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NoSlipBoundaryCondition NoSlipBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "Releigh-Taylor";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 1e-6;
	double TimeToSave = delta_t;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
		};
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
		};
		if (model.totalTime > 10000) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << TimeToSave;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

//2D Fortova settings var 1
void RunReleighTaylor2Dvar1(int N_x, int N_y)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation by velocity var1
	info.A = 1.0e2;
	info.g = 5.0e8;
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;
	RTI_proc proc;
	proc.task2D();
	proc.SetParameters(info);
	proc.SetInitConditionFunction(ReleighTaylorInit2D_var1);

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);
	proc.BindGrid(grid, N_x, N_y, 0);
	proc.BindModel(model);

	//write settigs in file
	proc.WriteTaskSettings();
	proc.IncrementRequired();

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(proc.GetPointerToFuncInitValues(), proc.GetInfoVar());
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	
	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "Releigh-Taylor";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 1e-6;
	double TimeToSave = delta_t;
	double Comp_Max_Time = 4.0e-5;
	int n_delta = 1;

	proc.SetTimersForSavings(model.totalTime, 0.01*delta_t);
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
			n_delta = (int)(model.totalTime/delta_t) + 1;
		};
		//increment part
		proc.IncrementWrite(model.totalTime);
		/*if((proc.isIncrementRequired == true)&&(model.totalTime > TimeToSaveInc)) {
			TimeToSaveInc += 0.01*delta_t;
			ofs.open("increment.dat", std::ofstream::out | std::ofstream::app);
			ofs << model.totalTime << '\t' << proc.GetTargetRhoPos() << '\n';
			ofs.close();
		};*/

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
		};
		if (model.totalTime > Comp_Max_Time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};
//2D Fortova settings var 2
void RunReleighTaylor2Dvar2(int N_x, int N_y)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation by velocity var1
	info.A = 1.0e2;		//not use in this test
	info.g = 5.0e8;
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;
	//set pointers for RTIprocessor 
	RTI_proc proc;
	proc.task2D();
	proc.SetParameters(info);
	proc.SetInitConditionFunction(ReleighTaylorInit2D_var2);

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);
	proc.BindGrid(grid, N_x, N_y, 0);
	proc.BindModel(model);

	//write settigs in file
	proc.WriteTaskSettings();
	proc.IncrementRequired();

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(proc.GetPointerToFuncInitValues(), proc.GetInfoVar());
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	
	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "Releigh-Taylor";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 1e-6;
	double TimeToSave = delta_t;
	double Comp_Max_Time = 4.0e-5;
	int n_delta = 1;

	proc.SetTimersForSavings(model.totalTime, 0);//0.01*delta_t);
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
			n_delta = (int)(model.totalTime/delta_t) + 1;
		};
		//increment part
		proc.IncrementWrite(model.totalTime);

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
		};
		if (model.totalTime > Comp_Max_Time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};


//3D test for Releigh-Taylor Instability with constant accelaration along y axis
void RunReleighTaylor3D(int N_x, int N_y, int N_z)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation
	/*info.A = 5.0e4;
	info.g = 5.0e8;
	info.gamma = 1.4;
	info.P0 = 2.0e12;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.03;
	info.y_max = 0.03;
	info.z_min = info.x_min;
	info.z_max = info.x_max;*/

	//Fe Pe riverse direction for g
	info.A = 5.0e4;
	info.g = -5.0e8;
	info.gamma = 1.4;
	info.P0 = 2.0e12;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.z_max = 0.03;
	info.z_min = -0.03;
	info.y_max = info.x_max;
	info.y_min = info.x_min;


	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid3D(N_x, N_y, N_z, info.x_min, info.x_max, info.y_min, info.y_max, info.z_min, info.z_max);

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, 0, -info.g);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(ReleighTaylorInit3D_discrete_perturbations_z_xy, &info);
	
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);
	model.SetBoundaryCondition("back", SymmetryBC);
	model.SetBoundaryCondition("front", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "Releigh-Taylor";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 1e-6;
	double TimeToSave = delta_t;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
		};
		model.Step();	
		if (i % 1 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
			//std::stringstream ss;
			//ss << i;
			//std::string iter = ss.str();
			//model.SaveToTechPlot(outputSolutionFile+iter+".dat");
		};
		if (model.totalTime > 10000) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << TimeToSave;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

//3D Fortova settings var 1
void RunReleighTaylor3Dvar1(int N_x, int N_y, int N_z)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation by velocity var1
	info.A = 1.0e2;
	info.g = 5.0e8;
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;
	info.z_max = 0.01;
	info.z_min = -0.01;
	RTI_proc proc;
	proc.task3D();
	proc.SetParameters(info);
	proc.SetInitConditionFunction(ReleighTaylorInit3D_var1);

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid3D(N_x, N_y, N_z, info.x_min, info.x_max, info.y_min, info.y_max, info.z_min, info.z_max);
	proc.BindGrid(grid, N_x, N_y, N_z);
	proc.BindModel(model);

	//write settigs in file
	proc.WriteTaskSettings();
	proc.IncrementRequired();

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(proc.GetPointerToFuncInitValues(), proc.GetInfoVar());
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	
	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);
	model.SetBoundaryCondition("front", SymmetryBC);
	model.SetBoundaryCondition("back", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "Releigh-Taylor";
	model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 5e-7;
	double TimeToSave = delta_t;
	double Comp_Max_Time = 4.0e-5;
	int n_delta = 1;

	proc.SetTimersForSavings(model.totalTime, 0.01*delta_t);
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
			n_delta = (int)(model.totalTime/delta_t) + 1;
		};
		//increment part
		proc.IncrementWrite(model.totalTime);

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
		};
		if (model.totalTime > Comp_Max_Time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};
//3D Fortova settings var 2
void RunReleighTaylor3Dvar2(int N_x, int N_y, int N_z)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation by velocity var1
	info.A = 1.0e2;
	info.g = 5.0e8;
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;
	info.z_max = 0.01;
	info.z_min = -0.01;
	RTI_proc proc;
	proc.task3D();
	proc.SetParameters(info);
	proc.SetInitConditionFunction(ReleighTaylorInit3D_var2);

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid3D(N_x, N_y, N_z, info.x_min, info.x_max, info.y_min, info.y_max, info.z_min, info.z_max);
	proc.BindGrid(grid, N_x, N_y, N_z);
	proc.BindModel(model);

	//write settigs in file
	proc.WriteTaskSettings();
	proc.IncrementRequired();

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(proc.GetPointerToFuncInitValues(), proc.GetInfoVar());
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	
	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);
	model.SetBoundaryCondition("back", SymmetryBC);
	model.SetBoundaryCondition("front", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "Releigh-Taylor";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 5e-7;
	double TimeToSave = delta_t;
	double Comp_Max_Time = 2.0e-5;
	int n_delta = 1;

	//increment savings
	proc.SetTimersForSavings(model.totalTime, 0.01*delta_t);
	bool isSave = true;

	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
			n_delta = (int)(model.totalTime/delta_t) + 1;
		};
		//increment part
		proc.IncrementWrite(model.totalTime);

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
		};
		if (model.totalTime > Comp_Max_Time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};





////Fe Pl plates collision part of tests////
ConservativeVariables FePuCollisionInit1D(Vector CellPos, void *params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;		// just one y axis
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0;

	//velocity
	v.x = 0;
	v.z = 0;
	if(CellPos.y <= 0) v.y = 0;
	else v.y = info.A;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};
bool RunFePuCollision1DTest(int N_cells, double time, double delta_time_to_save) {
	Model<Roe3DSolverPerfectGas> model;
	Grid grid;
	Rectangle_Area_RT info;
	info.A = -500;	//velocity of top plate
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.y_max = 0.03;
	info.y_min = -0.03;

	//Shock tube problem setting
	Vector direction = Vector(0,1,0);
	grid = GenGrid1D(N_cells, info.y_min, info.y_max, direction);
	
	//Set fluid properties	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.DisableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.005);
	model.SetSchemeOrder(2);

	//Bind computational grid
	model.BindGrid(grid);
	if(model.SchemeOrder > 1) model.EnableLimiter();
	else model.DisableLimiter();

	//Set initial conditions
	model.SetInitialConditions(FePuCollisionInit1D, &info);

	//free boundary condition
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition Symmetry(model);
	model.SetBoundaryCondition("left", NaturalBC);
	model.SetBoundaryCondition("right", Symmetry);

	//Load solution
	std::string outputSolutionFile = "collision1D";
	//model.LoadSolution("solution_to_load.txt");

	model.SaveToTechPlot("init.dat");

	int n_delta = (int)(model.totalTime/delta_time_to_save) + 1;
	double TimeToSave = n_delta*delta_time_to_save;
	//Total time
		for (int i = 0; i < 200000000; i++) {
		model.Step();
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_time_to_save;
			n_delta++;
		};
		if (model.totalTime > time) break;
	};

	//Save result to techplot
	model.SaveToTechPlot(outputSolutionFile + ".dat");
	return true;
};


//2D initializing of init conditions with constant speed of top plate (Pb)
ConservativeVariables FePlCollisionInit2D(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	if(CellPos.y <= 0) v.y = 0;
	else v.y = info.A;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};

//2D initializing of init conditions with constant speed of top plate (Pb) and distorted border
ConservativeVariables FePlCollision_1Mode_Init2D(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//common border
	double bound_pos = -info.delta*sin(PI*(0.5 + CellPos.x/(0.5*info.x_max)));
	//double bound_pos = -info.delta*sin(PI*(0.5 + CellPos.x/info.x_max));


	//density
	if(CellPos.y<=bound_pos) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	//if(CellPos.y <= 0) v.y = 0;
	if(CellPos.y <= bound_pos) v.y = 0;
	else v.y = info.A;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};

//2D initializing of init conditions with constant speed of top plate (Pb)
//and distorted density in two cells of top plate(постановка Фортовой)
ConservativeVariables FePlCollisionInit2DF(Vector CellPos, void* params) {
Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0 - ro*info.g*CellPos.y;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if(CellPos.y > 0) v.y = info.A;
	if((abs(CellPos.x)<4.0e-4)&&(abs(CellPos.y)<4.0e-4)&&(CellPos.y>0)) ro = 0.1*info.ro_top;
	
	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};

//2D initializing of init conditions for collision task with velocity distorted distribution
ConservativeVariables FePlCollision_Velocity_Distorted_Init2D(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;

	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//density
	if(CellPos.y<=0) ro = info.ro_bot;
	else ro = info.ro_top;

	//pressure
	p = info.P0;

	//velocity
	v.x = 0;
	v.z = 0;
	v.y = 0;
	if(CellPos.y > 0) v.y = info.A;
	if((abs(CellPos.x)<4.0e-4)&&(abs(CellPos.y)<4.0e-4)&&(CellPos.y>0)) v.y = 0.1*info.A;
	
	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;
};

//2D collision test
void RunCollisionTest2D(int N_x, int N_y, double delta_time_to_save, double Max_time)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pb discrete perturbation
	info.A = -5.0e2;
	info.g = 0.0;		//absolute value of accelaration
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	//info.ro_top = 11.34e3;	//Pb
	info.ro_top = 1.05*info.ro_bot;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();
	//model.DisableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(FePlCollision_Velocity_Distorted_Init2D, &info);
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "collision2D";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	int n_delta = (int)(model.totalTime/delta_time_to_save) + 1;
	double TimeToSave = n_delta*delta_time_to_save;
	
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
		};
		if (model.totalTime > Max_time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_time_to_save;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

//2D collision test with 500 m/s velosity of top (Pb) plate and Gausse acceleration
void RunCollisionTest2DGausseAcceleration(int N_x, int N_y, double delta_time_to_save, double Max_time)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation
	info.A = -5.0e2;
	info.g = 5.0e8;		//absolute value of accelaration
	info.gamma = 1.4;
	info.P0 = 1.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;
	info.delta = 0.0005;
	//Gausse distribution parameters
	double k = 0.01;	// k = exp(-0.5*t0*t0/dispertion)
	double t0 = 5.0e-7;	// expectation

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g = Vector(0, 0, 0);	//initialisation of acelaration depending on time

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();
	//model.DisableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(FePuCollisionInit1D, &info);
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", NaturalBC);
	model.SetBoundaryCondition("right", NaturalBC);
	model.SetBoundaryCondition("bottom", NaturalBC);
	model.SetBoundaryCondition("top", NaturalBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "collision2D";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	int n_delta = (int)(model.totalTime/delta_time_to_save) + 1;
	double TimeToSave = n_delta*delta_time_to_save;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		//compute acceleration
		g.y = -info.g*exp(log(k)*(model.totalTime/t0 - 1.0)*(model.totalTime/t0 - 1.0));
		model.SetUniformAcceleration(g);
		//Step
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
		};
		if (model.totalTime > Max_time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_time_to_save;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

//2D collision test with 500 m/s velosity of top (Pl) plate
//and density distorted in then one cell layer (постановка Фортовой)
void RunCollisionTest2DF(int N_x, int N_y)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pe discrete perturbation
	info.A = -5.0e2;	//velocity of top material
	info.g = 0.0;		
	info.gamma = 1.4;
	info.P0 = 2.0e11;
	info.ro_bot = 7.9e3;
	info.ro_top = 11.34e3;
	info.x_max = 0.01;
	info.x_min = -0.01;
	info.y_min = -0.01;
	info.y_max = 0.01;
	//set pointers for RTIprocessor 
	RTI_proc proc;
	proc.task2D();
	proc.SetParameters(info);
	proc.SetInitConditionFunction(FePlCollisionInit2DF);
	
	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);
	proc.BindGrid(grid, N_x, N_y, 0);
	proc.BindModel(model);

	//write settigs in file
	proc.WriteTaskSettings();

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(proc.GetPointerToFuncInitValues(), proc.GetInfoVar());
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", SymmetryBC);
	model.SetBoundaryCondition("right", SymmetryBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "collision2D";
	model.LoadSolution("solution_to_load.txt");

	//Run simulation
	double delta_t = 1e-6;
	double TimeToSave = delta_t;
	double Comp_Max_Time = 4.0e-5;
	int n_delta = 1;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		if((model.totalTime >= TimeToSave)&&(i==0)) {
			TimeToSave = ((int)(model.totalTime/delta_t) + 1)*delta_t;
			n_delta = (int)(model.totalTime/delta_t) + 1;
		};
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
		};
		if (model.totalTime > Comp_Max_Time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_t;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

#endif
#ifndef TURBO_RELEIGH_TAYLOR_PROCESSOR
#define TURBO_RELEIGH_TAYLOR_PROCESSOR

#include "tests.h"

//domain sizes and perturbation factor 
struct Rectangle_Area_RT {
	//coordinates of rectangle nodes
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	//amplitude of velocity perturbation
	double A;
	//amplituted of border disturbance
	double delta;
	//acceleration of field directed in -Y
	double g;
	//specific heat ratio
	double gamma;
	//zero altitude pressure
	double P0;
	//density bottom
	double ro_bot;
	//density top
	double ro_top;
};

//define type that means pointer to ConservativeVariables function of Vector and void*
typedef ConservativeVariables (*InitCondPointer)(Vector, void *);

//processor
class RTI_proc
{
private:
	//pointers to info structure and initial values function
	Rectangle_Area_RT* info;
	ConservativeVariables (*init_cond)(Vector, void *);

	//links to grid and model
	int Nx, Ny, Nz;
	Grid* grid;
	Model<Roe3DSolverPerfectGas>* model;

	//set of cells indexes for increment computations
	std::set<int> inc_set;

	//filenames
	std::string SettingsFileName;
	std::string IncrementFileName;

	//time interval for saving increment datum and time to save increment
	double DeltaTimeToSaveInc;
	double TimeToSaveInc;


public:
	//flag
	bool is3D;
	bool isIncrementRequired;

	//constructor by default
	RTI_proc() {
	isIncrementRequired = false;
	SettingsFileName = "Постановка.txt";
	IncrementFileName = "increment.dat";
	};

	//set filenames
	void SetSettingsFileName(char* filename) {
		SettingsFileName = filename;
	};
	void SetIncrementFileName(char* filename) {
		IncrementFileName = filename;
	};

	//set dimension
	void task3D() {
		is3D = true;
	};
	void task2D() {
		is3D = false;
	};

	//set appropriate pointers
	void SetParameters(Rectangle_Area_RT& _info) {
		info = &_info;
	};
	void SetInitConditionFunction(ConservativeVariables (&fun)(Vector, void *)) {
		init_cond = &fun;
	};

	//bind grid and model to proc structure
	void BindGrid(Grid &_grid, int N_x, int N_y, int N_z) {
		grid = &_grid;
		Nx = N_x;
		Ny = N_y;
		Nz = N_z;
	};
	void BindModel(Model<Roe3DSolverPerfectGas> &_m) {
		model = &_m;
	};

	//get pointers to InitValuesFunction and Info structure
	InitCondPointer GetPointerToFuncInitValues() {
		return init_cond;
	};
	Rectangle_Area_RT* GetInfoVar() {
		return info;
	};

	//write settings file
	void WriteTaskSettings() {
		std::ofstream ofs(SettingsFileName);
		ofs << "Grid: Nx = " << Nx << ", Ny = " << Ny << " and Nz = " << Nz << '\n';
		ofs << "info.A = " << info->A << '\n';
		ofs << "info.delta = " << info->delta << '\n';
		ofs << "info.g = " << info->g << '\n';
		ofs << "info.gamma = " << info->gamma << '\n';
		ofs << "info.P0 = " << info->P0 << '\n';
		ofs << "info.ro_bot = " << info->ro_bot << '\n';
		ofs << "info.ro_top = " << info->ro_top << '\n';
		ofs << "info.x_max = " << info->x_max << '\n';
		ofs << "info.x_min = " << info->x_min << '\n';
		ofs << "info.y_max = " << info->y_max << '\n';
		ofs << "info.y_min = " << info->y_min << '\n';
		ofs << "info.z_max = " << info->z_max << '\n';
		ofs << "info.z_min = " << info->z_min << '\n';
		ofs.close();
	}

	//processing functions for increment computing
	//get all central cells along y axis
	std::set<int> GetCentralAxialCells(double delta_x, double delta_z) {
		std::set<int> res;
		std::vector<Cell*> cells = grid->cells.getLocalNodes();
		for(int i=0; i<cells.size(); i++) {
			Vector P = cells[i]->CellCenter;
			if((P.x<0)||(P.z<0)) continue;
			if((P.x<=delta_x)&&(P.z<=delta_z)) res.insert(cells[i]->GlobalIndex);
		};
		return res;
	};
	void IncrementRequired() {
		isIncrementRequired = true;
		std::ofstream ofs(IncrementFileName);
		ofs.close();
		double delta_x = (info->x_max - info->x_min)/Nx;
		double delta_z;
		if(!is3D) delta_z = 1;
		else delta_z = (info->x_max - info->x_min)/Nz;
		inc_set = GetCentralAxialCells(0.75*delta_x, 0.75*delta_z);
	};
	//get y Position of target rho value in central line of cells (along y axis)
	double GetTargetRhoPos() {
		//by default
		double rho_target = 0.5*(info->ro_bot + info->ro_top);
		//positions of two cells with rho that is closest to target value
		double val_lower = -rho_target;
		double pos_l = info->y_min;
		double val_higher = -rho_target;
		double pos_h = info->y_max;

		//DEBUG
		std::ofstream ofs("srez.dat");
		for(int i : inc_set) {
			ofs << grid->cells[i].CellCenter.y << ' ' << model->U[i].ro << '\n';
		};ofs.close();

		//find two closest values
		for(int i : inc_set) {
			double rho = model->U[i].ro;
			double y = grid->cells[i].CellCenter.y;
			if(rho < rho_target) {
				if(pos_l < y) {
					pos_l = y;
					val_lower = rho;
				};
			} else {
				if(pos_h > y) {
					pos_h = y;
					val_higher = rho;
				};
			};
		};

		//compute approximate position of rho_target
		double res = pos_l;
		double der = (val_higher - val_lower)/(pos_h - pos_l);
		res += (rho_target - val_lower)/der;

		return res;
	};
	//set time to save increment and time interval of this saving
	void SetTimersForSavings(double total_time, double delta_t) {
		DeltaTimeToSaveInc = delta_t;
		TimeToSaveInc = ((int)(total_time/DeltaTimeToSaveInc) + 1)*delta_t;
	};
	//compute and write increment for during TimeStep
	void IncrementWrite(double total_time) {
		if(isIncrementRequired == true)
		{
			if(total_time > TimeToSaveInc) {
			TimeToSaveInc += DeltaTimeToSaveInc;
			std::ofstream ofs;
			ofs.open(IncrementFileName, std::ofstream::out | std::ofstream::app);
			ofs << total_time << '\t' << GetTargetRhoPos() << '\n';
			ofs.close();
			};
		};
	};
};




#endif
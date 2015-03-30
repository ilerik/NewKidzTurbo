#ifndef NewKidzTurbo_TestCases_EOSTests_TestCase1D_BarotropicEOS
#define NewKidzTurbo_TestCases_EOSTests_TestCase1D_BarotropicEOS

#include "TestCase.h"
#include "Gengrid1DUniform.h"
#include <memory>

namespace TestCasesEOS {

enum class BoundaryConditionType {
	FreeSurface,
	Natural,
	Symmetry,
	Wall
};

class TestCase1D_BarotropicEOS : public InitialConditions::InitialConditions {
private:
	//Data
	Godunov3DSolverPerfectGas _exactSolver;
	std::shared_ptr<ParallelManager> _MPIManager;
	std::unique_ptr<Kernel> _kernel;
	std::shared_ptr<Grid> _grid;
	Configuration _configuration;	

	//Generate grid
	void _prepareGrid(int nCells) {
		_grid = GenGrid1DUniform((_MPIManager), nCells, 0.0, 1.0, false);
	};

	//Generate initial conditions

	//Initial materials distribution
	virtual int getInitialGasModelIndex(const Cell& cell) override  {
		//Cell center
		double x = cell.CellCenter.x;		

		//Split in 2 halfs		
		if (x  < 0.5) {
			return 0; //Left
		} else {
			return 1; //Right
		};		

		return 0;
	};

	//Initial values distribution
	virtual std::vector<double> getInitialValues(const Cell& cell) override	{
		int nmat = getInitialGasModelIndex(cell); //get material index
		std::vector<double> initValues;		

		//Velocities
		double u = 0;
		double v = 0;
		double w = 0;

		//Internal energy
		double e = 0;		

		//Cell center
		double x = cell.CellCenter.x;		
				
		//Values
		double roE = 0;		
		double ro = 0;
		double pressure = 0.0;		
		if (nmat == 0) {
			//Light
			ro = 7900;
			u = 250;
			pressure = 2e14;
			e = 0;
		} else {
			//Heavy
			ro = 11400;
			u = -250;
			pressure = 2e14;
			e = 0;
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

	//Generate configuration
	GasModelConfiguration GetGasModelConfiguration(double E, double ro0, double p0) {
		GasModelConfiguration conf;

		conf.GasModelName = "BarotropicGasModel";
		conf.SetPropertyValue("YoungModulus", E);
		conf.SetPropertyValue("ReferencePressure", p0);
		conf.SetPropertyValue("ReferenceDensity", ro0);		
		
		return conf;		
	};

	BoundaryConditionConfiguration GetBoundaryConditionConfiguration(std::string materialName, BoundaryConditionType type) {
		BoundaryConditionConfiguration conf;

		if (type == BoundaryConditionType::Symmetry) {
			conf.BoundaryConditionType = BCType_t::BCSymmetryPlane;
			conf.MovementType = BoundaryConditionMovementType::Fixed;
			conf.MaterialName = materialName;
			return conf;
		};

		if (type == BoundaryConditionType::Natural) {
			conf.BoundaryConditionType = BCType_t::BCOutflowSupersonic;
			conf.MovementType = BoundaryConditionMovementType::Fixed;
			conf.MaterialName = materialName;
			return conf;			
		};

		if (type == BoundaryConditionType::FreeSurface) {			
			conf.BoundaryConditionType = BCType_t::BCGeneral;
			conf.MovementType = BoundaryConditionMovementType::FreeSurface;
			conf.MaterialName = materialName;
			return conf;
		};

		throw new Exception("Unspecified boundary condition type");
		return conf;
	};

	void _prepareConfiguration() {
		//Hardcode configuration for now
		_configuration.WorkingDirectory = "./results";
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Rieman solver settings
		_configuration.RiemannSolverConfiguration.riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::HLLC;

		//Availible gas models		
		//Left metal
		_configuration.AddGasModel("Light");
		_configuration.GasModelsConfiguration["Light"] = GetGasModelConfiguration(2e13, 7900, 2e14);

		//Right metal
		_configuration.AddGasModel("Heavy");
		_configuration.GasModelsConfiguration["Heavy"] = GetGasModelConfiguration(2e13, 11400, 2e14);

		//Boundary conditions		
		_configuration.BoundaryConditions["left"] = GetBoundaryConditionConfiguration("Light", BoundaryConditionType::Natural);
		_configuration.BoundaryConditions["right"] = GetBoundaryConditionConfiguration("Heavy", BoundaryConditionType::Natural);			
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.CFL = 0.1;
		_configuration.RungeKuttaOrder = 1;		

		//ALE settings
		_configuration.ALEConfiguration.MeshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDWnoRotation;
		//_configuration.ALEConfiguration.ALEMotionType = "Eulerian";	
		//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";
		_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";

		//Run settings
		_configuration.MaxIteration = 100000000;
		_configuration.MaxTime = 5e-6;
		_configuration.SaveSolutionSnapshotIterations = 0;
		_configuration.SaveSolutionSnapshotTime = 1e-6;

		//Gravity
		_configuration.g = Vector(0,0,0);				
	};

	void _prepareKernel() {
		//Initialize
		_kernel = std::unique_ptr<Kernel>(new Kernel());
		_kernel->Initilize(_MPIManager->comm());

		//Bind grid
		_kernel->BindGrid(_grid);

		//Bind configuration
		_kernel->VerboseOn();
		_kernel->BindConfiguration(_configuration);	
	};
	
public:
	void Init(int argc, char *argv[]) {
		//Init MPI
		_MPIManager = std::unique_ptr<ParallelManager>(new ParallelManager(argc, argv));
	};

	void Run(int nCells) {
		//Init kernel
		_prepareGrid(nCells);
		_prepareConfiguration();
		_prepareKernel();

		//Initialize calculation
		_kernel->InitCalculation();

		//Initial conditions
		_kernel->GenerateInitialConditions(this);
		//_kernel->SaveGrid("init.cgns");
		//_kernel->SaveSolution("init.cgns", "Solution");

		//Run computational cycle
		_kernel->RunCalculation();
		//_kernel->SaveGrid("result.cgns");
		//_kernel->SaveSolution("result.cgns", "Solution");

		//Finilize
		_kernel->FinalizeCalculation();
		_kernel->Finalize();

		_MPIManager->Barrier();
	};

	//Save result to tecplot
	void SaveSolutionTecplot() {		
		std::ofstream ofs("solution.dat");
		ofs<<"VARIABLES = \"x\" \"ro\"  \"p\"  \"u\" \"c\""<<std::endl;

		double L2error = 0.0;
		for (int i = 0; i<_grid->nCellsLocal; i++) {
			Cell& cell = _kernel->GetLocalCell(i);			
			double x = cell.CellCenter.x;			

			//Local gas model
			int nmat = _kernel->GetCellGasModelIndex(cell.GlobalIndex);
			GasModel* gasModelPtr = _gasModels[nmat].get();
		
			//Obtained solution
			std::vector<double> U = _kernel->GetCellValues(cell.GlobalIndex);
			double ro = U[0];
			double u = U[1] / U[0];
			double e = U[4] / ro - u*u/2.0;
			double p = 0;
			double c = 0;
			double Gr = 0;
			gasModelPtr->GetPressureAndSoundSpeed(U, p, c, Gr);

			//Output to file
			ofs<<x<<" ";
			ofs<<ro<<" ";			
			ofs<<p<<" ";			
			ofs<<u<<" ";
			ofs<<c<<" ";
			ofs<<std::endl;
		};		

		ofs.close();		
	};
};

};

#endif
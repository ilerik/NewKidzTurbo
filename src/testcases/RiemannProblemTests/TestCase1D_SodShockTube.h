#ifndef TURBO_TestCases_RiemannProblemTests_TestCase1D_SodShockTube
#define TURBO_TestCases_RiemannProblemTests_TestCase1D_SodShockTube

#include "TestCase.h"
#include "Gengrid1DUniform.h"
#include <memory>

namespace TestCases1D {

enum class BoundaryConditionType {
	FreeSurface,
	Natural,
	Symmetry,
	Wall
};

class TestCase1D_SodShockTube : public InitialConditions::InitialConditions {
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
			ro = 1.0;
			pressure = 1.0;
			e = pressure / (0.4 * ro);
		} else {
			//Heavy
			ro = 0.125;
			pressure = 0.1;
			e = pressure / (0.4 * ro);
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
	GasModelConfiguration GetGasModelConfiguration() {
		GasModelConfiguration conf;

		conf.GasModelName = "PerfectGasModel";
		conf.SetPropertyValue("IdealGasConstant", 8.3144621);
		conf.SetPropertyValue("SpecificHeatRatio", 1.4);
		conf.SetPropertyValue("SpecificHeatVolume", 1006.43 / 1.4);
		conf.SetPropertyValue("SpecificHeatPressure", 1006.43);
		
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
		_configuration.GasModelsConfiguration["Light"] = GetGasModelConfiguration();

		//Right metal
		_configuration.AddGasModel("Heavy");
		_configuration.GasModelsConfiguration["Heavy"] = GetGasModelConfiguration();

		//Boundary conditions		
		_configuration.BoundaryConditions["left"] = GetBoundaryConditionConfiguration("Light", BoundaryConditionType::Symmetry);
		_configuration.BoundaryConditions["right"] = GetBoundaryConditionConfiguration("Heavy", BoundaryConditionType::Symmetry);			
		
		//Solver settings					
		_configuration.SimulationType = TimeAccurate;
		_configuration.SpatialDiscretisation = SpatialDiscretisationType::ENO;
		_configuration.CFL = 0.1;
		_configuration.RungeKuttaOrder = 1;		

		//ALE settings
		_configuration.ALEConfiguration.MeshMovementAlgorithm = MeshMovement::MeshMovementAlgorithm::IDWnoRotation;
		_configuration.ALEConfiguration.ALEMotionType = "Eulerian";	
		//_configuration.ALEConfiguration.ALEMotionType = "Lagrangian";
		//_configuration.ALEConfiguration.ALEMotionType = "ALEMaterialInterfaces";

		//Run settings
		_configuration.MaxIteration = 100000000;
		_configuration.MaxTime = 0.2;
		_configuration.SaveSolutionSnapshotIterations = 1;
		_configuration.SaveSolutionSnapshotTime = 0;

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

	//Get error value
	double L2error() {
		_exactSolver.BindGasModels(_gasModels);
		GasModel::ConservativeVariables UL;
		GasModel::ConservativeVariables UR;
		UL.ro = 1.0;
		UL.rou = 0.0;
		UL.rov = 0.0;
		UL.row = 0.0;
		UL.roE = 1.0 / (0.4);

		UR.ro = 0.125;
		UR.rou = 0.0;
		UR.rov = 0.0;
		UR.row = 0.0;
		UR.roE = 0.1 / (0.4);

		Face f;
		f.FaceNormal = Vector(1.0, 0.0, 0.0);
		f.FaceCenter = Vector(0.0, 0.0, 0.0);
		f.FaceSquare = 1.0;

		std::ofstream ofs("solution.dat");
		ofs<<"VARIABLES = \"x\" \"ro\" \"roExact\" \"p\" \"pExact\" \"u\" \"uExact\""<<std::endl;

		double L2error = 0.0;
		for (int i = 0; i<_grid->nCellsLocal; i++) {
			Cell& cell = _kernel->GetLocalCell(i);
			double dx = -0.5;
			double x = cell.CellCenter.x + dx;
			double xL = _kernel->GetLocalFace(cell.Faces[0]).FaceCenter.x + dx;
			double xR = _kernel->GetLocalFace(cell.Faces[1]).FaceCenter.x + dx;
			double t = 0.2;
			std::vector<double> U1 = _exactSolver.SampleSolution( UL, UR, f, xL/t);
			std::vector<double> U2 = _exactSolver.SampleSolution( UL, UR, f, xR/t);

			//Exact second order solution
			std::vector<double> UExact = 0.5 * (U1 + U2);
			double roExact = UExact[0];
			double uExact = UExact[1] / UExact[0];
			double eExact = UExact[4] / roExact - uExact*uExact/2.0;
			double pExact = eExact * roExact / (0.4); 

			//Obtained solution
			std::vector<double> U = _kernel->GetCellValues(cell.GlobalIndex);
			double ro = U[0];
			double u = U[1] / U[0];
			double e = U[4] / ro - u*u/2.0;
			double p = e * ro / (0.4); 
			L2error += std::pow(ro - roExact, 2.0) * cell.CellVolume;

			//Output to file
			ofs<<x<<" ";
			ofs<<ro<<" ";
			ofs<<roExact<<" ";
			ofs<<p<<" ";
			ofs<<pExact<<" ";
			ofs<<u<<" ";
			ofs<<uExact<<" ";
			ofs<<std::endl;
		};
		L2error = std::sqrt(L2error);

		ofs.close();

		return L2error;
	};
};

};

#endif
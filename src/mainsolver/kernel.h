#ifndef TURBO_MAINSOLVER_KERNEL
#define TURBO_MAINSOLVER_KERNEL

#include "grid.h"
#include "mpi.h"
#include "logger.h"
#include "parmetis.h"
#include "CGNSReader.h"
#include "CGNSWriter.h"
#include "parallelHelper.h"
#include "perfomanceHelper.h"
#include "BoundaryConditions.h"
#include "InitialConditions.h"
#include "cmath"
#include <direct.h>
#include "geomfunctions.h"
#include "configuration.h"
#include "riemannsolvers.h"
#include "LomonosovFortovGasModel.h"
#include "BarotropicGasModel.h"
#include "PerfectGasModel.h"
#include "meshquality.h"
#include "ALEMethod.h"
#include "sources/GravitySource.h"
#include "mainsolver/SpatialDiscretisation/SpatialDiscretisation.h"

//Define error types
enum turbo_errt {
	TURBO_OK = 0,
	TURBO_ERROR = 1
};

//Step info
class StepInfo {
public:
	double Time;
	double TimeStep;	
	int Iteration;
	std::vector<double> Residual;

	//Perfomance times report
	double ConvectiveFluxesTime;
	double ComputeResidualTime;
	double ALETime;

};


//Forward declaration
class Kernel;

//Base class for storing activity executed after each computational step
class StepHistoryLogger {
protected:
	Kernel* _kernel;
	Logger* _logger;
	Grid* _grid;
	ParallelHelper* _parallelHelper;
	std::vector<GasModel*>* _gasModels;	
public:
	void setStepHistoryLoggerReferences(Kernel* kernel, Logger* logger, Grid* grid, ParallelHelper* parallelHelper, std::vector<GasModel*>* gasModels) {
		_kernel = kernel;
		_logger = logger;
		_grid = grid;
		_parallelHelper = parallelHelper;
		_gasModels = gasModels;
	};

	virtual void Init() = 0;
	virtual void SaveHistory() = 0;
	virtual void Finalize() = 0;
};

//Calculation kernel
class Kernel {
private:
	Kernel& operator=(const Kernel&); //non copyable

	//Logger
	std::string _logfilename;
	Logger _logger;

	//Helper for CGNS i\o
	GridLoading::CGNSReader _cgnsReader;
	PostProcessing::CGNSWriter _cgnsWriter;

	//History logger object
	StepHistoryLogger* _stepHistoryLogger;

	//Grid data
	std::shared_ptr<Grid> _gridPtr;

	//Configuration
	Configuration _configuration;

	//Spatial scheme
	SpatialDiscretisationType _spatialDiscretisation;

	//Gas model (equations of state)
	std::vector<std::shared_ptr<GasModel> > _gasModels;	

	//Sources (external forces, and so on)
	std::unique_ptr<GravitySource> _sources;

	//ALE data and logic
	ALEMethod _ALEmethod;

	//Solver (TO DO)
	SimulationType_t _simulationType; //Simulation type
	bool _isVerbose;				  //If standart output is on during calculation
	double CFL;
	int RungeKuttaOrder;
	double MaxIteration;
	double MaxTime;
	double CurrentTime;
	double NextSnapshotTime;
	double SaveSolutionSnapshotTime;
	int SaveSolutionSnapshotIterations;	

	//Riemann solver
	RiemannSolver* rSolver;

	//Step information
	StepInfo stepInfo;

	//Parallel run information		
	ParallelHelper _parallelHelper;
	int _rank;
	int _nProcessors;

	//Perfomance measurement helper
	//PerfomanceHelper _perfomanceHelper;
	Timer timerResidual;
	Timer timerALE;
	Timer timerConvective;

	//Boundary conditions
	std::map<int, int> _boundaryGasModelIndex;
	std::map<int, std::unique_ptr<BoundaryConditions::BoundaryCondition> > _boundaryConditions;

	//Internal storage		
	int nVariables; //number of variables in each cell
	std::vector<int> CellGasModel; //Index of gas model for cell constituents (localCell index -> GasModel index)
	std::vector<double> Values;	//Conservative variables
	std::vector<double> Residual; //Residual values	

	//Fluxes
	std::vector<std::vector<double>> FaceFluxes;
	std::vector<double> MaxWaveSpeed; //Maximal wave speed estimate for face
	std::vector<double> FacePressure; //Pressure at each face	

	//Gradients
	std::vector<Vector> GradientT;
	std::vector<Vector> GradientP;

public:	
	//Parallel helper
	inline ParallelHelper* getParallelHelper() {
		return &_parallelHelper;
	};

	//Step info
	inline StepInfo* getStepInfo() {
		return &stepInfo;
	};

	//History logger
	inline void setStepHistoryLogger(StepHistoryLogger* logger) {
		return;

		//Delete if exists
		if (_stepHistoryLogger != NULL) {
			_stepHistoryLogger->Finalize();
			delete _stepHistoryLogger;
		};

		//Properly attach object to kernel
		_stepHistoryLogger = logger;
		//_stepHistoryLogger->setStepHistoryLoggerReferences(this, &_logger, &_grid, &_parallelHelper, &_gasModels);

		//Synchronize
		_parallelHelper.Barrier();
	};

	//Initialize computational kernel
	turbo_errt Initilize(int *argc, char **argv[]) {
		MPI_Init(argc, argv);
		Initilize(MPI_COMM_WORLD);
	};

	turbo_errt Initilize(MPI_Comm comm) {
		//Initialize parallel MPI subsystem
		_parallelHelper.Init(comm);

		try {
			_rank = _parallelHelper.getRank();
			_nProcessors = _parallelHelper.getProcessorNumber();

			//Initialize perfomance watching subsystem
			//_perfomanceHelper.Init(_parallelHelper);

			//Initialize loggin subsystem
			_logfilename = "kernel"; //TO DO
			_logger.InitLogging(_parallelHelper, _logfilename);

			//Initialize cgns i\o subsystem
			_cgnsReader.Init(_logger, _parallelHelper);
			_cgnsWriter.Init(_logger, _parallelHelper);

			//No default history logger
			_stepHistoryLogger = NULL;
		} catch (...) {
			//Exception was thrown but shouldnt
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Kernel initialization failed");
			exit(0);
		};

		//Finish successfully
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Kernel initialization completed successfully");
		return TURBO_OK;
	};

	//Finilize kernel work
	turbo_errt Finalize() {
		_cgnsReader.Finalize();
		_cgnsWriter.Finalize();
		_parallelHelper.Finalize();			
		_logger.FinilizeLogging();
		//_perfomanceHelper.Finalize();
		return TURBO_OK;
	};

	turbo_errt InitCalculation() {
		//Change working directory
		if (_configuration.WorkingDirectory != "") {
			//const wchar_t *ptr = _configuration.WorkingDirectory.c_str();
			//_wchdir();
		};

		//Gas models setup
		//List of availible gas models
		_gasModels.resize(_configuration.GasModelsConfiguration.size());
		for (std::pair<std::string, GasModelConfiguration> pair : _configuration.GasModelsConfiguration) {
			int gmIndex = _configuration.GasModelNameToIndex[pair.first];
			std::string gmName = pair.second.GasModelName;
			if (gmName == "PerfectGasModel") {
				_gasModels[gmIndex] = std::unique_ptr<GasModel>(new PerfectGasModel(_logger));				
			};
			if (gmName == "LomonosovFortovGasModel") {
				_gasModels[gmIndex] = std::unique_ptr<GasModel>(new LomonosovFortovGasModel(_logger));
			};
			if (gmName == "BarotropicGasModel") {
				_gasModels[gmIndex] = std::unique_ptr<GasModel>(new BarotropicGasModel(_logger));
			};

			//Load configuration for gas model
			_gasModels[gmIndex]->loadConfiguration(pair.second);
		};

		//TO DO set shifts instead of global nVariables constant
		for (std::shared_ptr<GasModel>& gasModel : _gasModels) {			
			nVariables = gasModel->nConservativeVariables; //TO DO generalize for now equal number of conservative variables assumed
		};

		//Allocate memory for data structures
		Values.resize(nVariables * _gridPtr->nCellsLocal);
		Residual.resize(nVariables * _gridPtr->nCellsLocal);
		CellGasModel.resize(_gridPtr->nCellsLocal);
		FaceFluxes.resize(_gridPtr->nFaces);
		for (std::vector<double>& flux : FaceFluxes) {
			flux.resize(nVariables);
			for (double& v : flux) v = 0;
		};	
		MaxWaveSpeed.resize(_gridPtr->nFaces);		
		FacePressure.resize(_gridPtr->nFaces);

		//Solver settings
		bool isSolverSpecified = false;
		//Instantiate solver depending on user choice
		if (_configuration.RiemannSolverConfiguration.riemannSolverType == RiemannSolverConfiguration::RiemannSolverType::HLLC) {
			rSolver = new HLLCSolverGeneralEOS();
			isSolverSpecified = true;
		};
		if (_configuration.RiemannSolverConfiguration.riemannSolverType == RiemannSolverConfiguration::RiemannSolverType::Roe) {
			rSolver = new Roe3DSolverPerfectGas();
			isSolverSpecified = true;
		};
		if (_configuration.RiemannSolverConfiguration.riemannSolverType == RiemannSolverConfiguration::RiemannSolverType::Godunov) {
			rSolver = new Godunov3DSolverPerfectGas();
			isSolverSpecified = true;
		};
		if (!isSolverSpecified) {
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Riemann solver is not specified.");
			//Halt programm		
			FinalizeCalculation();
			Finalize();
			exit(0);
		};
		
		//Try to bind gas models to solver
		if (!rSolver->BindGasModels(_gasModels)) {
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Riemann solver doesn't support specified gas models.");
			//Halt programm		
			FinalizeCalculation();
			Finalize();
			exit(0);
		};
		//Try to load configuration for Riemann solver
		if (!rSolver->loadConfiguration(&_logger, _configuration.RiemannSolverConfiguration)) {
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Riemann solver configuration load unsuccesful.");
			//Halt programm		
			FinalizeCalculation();
			Finalize();
			exit(0);
		};

		//Spatial discretisation
		_spatialDiscretisation = _configuration.SpatialDiscretisation;

		//Simulation settings
		_simulationType = _configuration.SimulationType;
		CFL = _configuration.CFL;
		RungeKuttaOrder = _configuration.RungeKuttaOrder;		

		//ALE settings
		bool isALEMovementSet = false;
		_ALEmethod._moveHelper.meshMovementAlgorithm = _configuration.ALEConfiguration.MeshMovementAlgorithm;
		if (_configuration.ALEConfiguration.ALEMotionType == "Eulerian") {
			_ALEmethod.ALEMotionType = ALEMethod::ALEMotionType::PureEulerian;		
			isALEMovementSet = true;
		};
		if (_configuration.ALEConfiguration.ALEMotionType == "ALEMaterialInterfaces") {
			_ALEmethod.ALEMotionType = ALEMethod::ALEMotionType::ALEMaterialInterfaces;	
			isALEMovementSet = true;
		};		
		if (_configuration.ALEConfiguration.ALEMotionType == "Lagrangian") {		
			_ALEmethod.ALEMotionType = ALEMethod::ALEMotionType::PureLagrangian;	
			isALEMovementSet = true;
		};
		_ALEmethod.SetGrid(_gridPtr);		

		//Run settings
		MaxIteration = _configuration.MaxIteration;
		MaxTime = _configuration.MaxTime;
		SaveSolutionSnapshotIterations = _configuration.SaveSolutionSnapshotIterations;
		SaveSolutionSnapshotTime = _configuration.SaveSolutionSnapshotTime;	

		//Initialize start moment
		stepInfo.Time = 0.0; 
		NextSnapshotTime = stepInfo.Time;

		//Initial conditions
		//InitialConditions::InitialConditions ic;
		//GenerateInitialConditions(ic);

		//Boundary conditions
		if (InitBoundaryConditions() == turbo_errt::TURBO_ERROR) {
			//Halt programm		
			FinalizeCalculation();
			Finalize();
			exit(0);
		};

		//Initialize source terms
		_sources = std::unique_ptr<GravitySource>(new GravitySource(_configuration.g) );

		//Init parallel exchange
		InitParallelExchange();	

		//Initialize history logger
		if (_stepHistoryLogger != NULL) {
			_stepHistoryLogger->Init();
		};

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished calculation initialization");	
		return TURBO_OK;
	};

	turbo_errt RunCalculation() {
		//Start phase

		//Calculate snapshot times order of magnitude
		int snapshotTimePrecision = 0;
		if (SaveSolutionSnapshotTime > 0) {
			snapshotTimePrecision = 1 - std::floor(std::log10(SaveSolutionSnapshotTime));
		};

		//Calc loop start
		std::stringstream msg;	
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Calculation started!");
		//Create phases for time step calculations
		for (stepInfo.Iteration = 0; stepInfo.Iteration <= MaxIteration; stepInfo.Iteration++) {
			//Perform timestep
			ExplicitTimeStep();

			//Solution snapshots

			//Every few iterations
			if ((SaveSolutionSnapshotIterations != 0) && (stepInfo.Iteration % SaveSolutionSnapshotIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<"dataI"<<stepInfo.Iteration<<".cgns";
				SaveGrid(snapshotFileName.str());				
				SaveSolution(snapshotFileName.str(), "Solution");
			};

			//Every fixed time interval
			if ((SaveSolutionSnapshotTime > 0) && (NextSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<std::fixed;
				snapshotFileName.precision(snapshotTimePrecision);								
				snapshotFileName<<"dataT"<<stepInfo.Time<<".cgns";
				SaveGrid(snapshotFileName.str());				
				SaveSolution(snapshotFileName.str(), "Solution");

				//Adjust next snapshot time
				NextSnapshotTime += SaveSolutionSnapshotTime;
			};

			//Save history
			if (_stepHistoryLogger != NULL) _stepHistoryLogger->SaveHistory();

			//Output step information
			msg.clear();
			msg.str(std::string());
			msg<<"Iteration = "<<stepInfo.Iteration<<"; Total time = "<< stepInfo.Time << "; Time step = " <<stepInfo.TimeStep << "; RMSrou = "<<stepInfo.Residual[1]<<"\n";
			//msg<<"ALEtoFluxRation = "<<stepInfo.ALETime / stepInfo.ConvectiveFluxesTime <<"; ConvectiveFluxTime = "<<stepInfo.ConvectiveFluxesTime << "; ALETime = " << stepInfo.ALETime<<"\n";
			//msg<<"Computation time = "<<computationPhase->GetTotalTimeMilliseconds()<<std::endl;
			//msg<<"Snapshot time = "<<snapshotsPhase->GetTotalTimeMilliseconds()<<std::endl;
			//msg<<"Save history time = "<<saveHistoryPhase->GetTotalTimeMilliseconds()<<std::endl;
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());
			if (_parallelHelper.IsMaster() && _isVerbose) {
				std::cout<<msg.str();
			};

			//Convergence criteria
			if (stepInfo.Iteration == MaxIteration) {
				msg.clear();
				msg.str(std::string());
				msg<<"Maximal number of iterations reached.";
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, msg.str());
				break;
			};

			if (stepInfo.Time >= MaxTime) {
				msg.clear();
				msg.str(std::string());
				msg<<"Maximal time reached.";
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, msg.str());
				break;
			};

			//Synchronize
			_parallelHelper.Barrier();
		};

		//Synchronize
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Calculation finished!");	

		//Output total run time
		msg.str(std::string());
		/*msg<<runCalculationPhase->GetTotalTimeMilliseconds();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Total work time = " + msg.str() + " seconds.");	*/
		return TURBO_OK;
	};

	turbo_errt FinalizeCalculation() {
		//History logging
		if (_stepHistoryLogger != NULL) {
			_stepHistoryLogger->Finalize();
			delete _stepHistoryLogger;
		};

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished calculation finalization.");	
		return TURBO_OK;
	};

	//Load grid from CGNS file
	turbo_errt LoadGrid(std::string filename) {
		//_gridPtr = _cgnsReader.LoadGrid(filename);
		//PartitionGrid(_gridPtr);
		//GenerateGridGeometry(_gridPtr);
		return TURBO_OK;
	};	

	//Save grid to CGNS file
	turbo_errt SaveGrid(std::string filename)
	{
		_cgnsWriter.CreateFile(filename);
		_cgnsWriter.WriteGridToFile(_gridPtr);		
		_cgnsWriter.CloseFile();
		return TURBO_OK;
	}

	//Bind existing grid to kernel
	turbo_errt BindGrid(std::shared_ptr<Grid>& grid) {
		_gridPtr = grid;
		return TURBO_OK;
	};

	//Partion loaded grid
	turbo_errt PartitionGrid(Grid& grid) {
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Grid partitioning started");	

		//Partitioning settings
		idx_t *vwgt = NULL; //No vertex weights for now
		idx_t *adjwgt = NULL; //No edge weights for now
					
		/* This is used to indicate if the graph is weighted. wgtflag can take one of four values:
		0 No weights (vwgt and adjwgt are both NULL).
		1 Weights on the edges only (vwgt is NULL).
		2 Weights on the vertices only (adjwgt is NULL).
		3 Weights on both the vertices and edges. 	*/
		idx_t wgtflag = 0; 

		idx_t numflag = 0; //C-style numeration

		/* This is used to specify the number of weights that each vertex has. It is also the number of balance
		constraints that must be satisfied.	*/
		idx_t ncon = 1;

		idx_t nparts = _nProcessors; //Number of subdomains desired

		/* An array of size ncon  nparts that is used to specify the fraction of vertex weight that should
		be distributed to each sub-domain for each balance constraint. If all of the sub-domains are to be of
		the same size for every vertex weight, then each of the ncon  nparts elements should be set to
		a value of 1/nparts. If ncon is greater than 1, the target sub-domain weights for each sub-domain
		are stored contiguously (similar to the vwgt array). Note that the sum of all of the tpwgts for a
		give vertex weight should be one. */
		std::vector<real_t> tpwgts(ncon*nparts);
		for (int i = 0; i<nparts; i++) tpwgts[i] = 1.0/nparts;

		/* An array of size ncon that is used to specify the imbalance tolerance for each vertex weight, with 1
		being perfect balance and nparts being perfect imbalance. A value of 1.05 for each of the ncon
		weights is recommended. */
		std::vector<real_t> ubvec(ncon); // imbalance tolerance for each vertex weight,
		ubvec[0] = 1.02;	// recomended default value 

		//Algoritm options for displaing information
		idx_t options[3];
		options[0] = 0; //No information and default values

		idx_t edgecut; // Upon successful completion, the number of edges that are cut by the partitioning is written to this parameter.

		/* This is an array of size equal to the number of locally-stored vertices. Upon successful completion the
		partition vector of the locally-stored vertices is written to this array. (See discussion in Section 4.2.4). */		 
		std::vector<idx_t> part(grid.nCellsLocal);
		std::ostringstream msg;
		/*msg.str(std::string());
		msg<<"vdist[] = \n";
		for (int i = 0; i<=_nProcessors; i++) msg<<grid.vdist[i]<<" ";
		msg<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg.str(std::string());
		msg<<"xadj[] = \n";
		for (int i = 0; i<grid.xadj.size(); i++) msg<<grid.xadj[i]<<" ";
		msg<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg.str(std::string());
		msg<<"adjncy[] = \n";
		for (int i = 0; i<grid.adjncy.size(); i++) msg<<grid.adjncy[i]<<" ";
		msg<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	*/

		//Call partitioning function		
		MPI_Comm _comm = _parallelHelper.getComm();
		_parallelHelper.Barrier();
		if (_nProcessors < grid.nProperCells) { //Make sure that partitioning is needed
			int result = ParMETIS_V3_PartKway(&grid.vdist[0], &grid.xadj[0], grid.adjncy._Myfirst, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, &tpwgts[0], &ubvec[0], options, &edgecut, &part[0], &_comm);			
			if (result != METIS_OK) {
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "ParMETIS_V3_PartKway failed.");
				return TURBO_ERROR;
			};							
		};

		//Gather partitioning on every processor
		std::vector<int> recvcounts(_nProcessors);
		for (int i = 0; i<_nProcessors; i++) {
			recvcounts[i] = grid.vdist[i+1] - grid.vdist[i];
		};
		
		/*msg.clear();
		msg<<part.size()<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg<<"recvcounts[] = \n";
		for (int i = 0; i<_nProcessors; i++) msg<<recvcounts[i]<<" ";
		msg<<"\n";

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	*/
		_parallelHelper.Allgatherv( part, recvcounts, grid.cellsPartitioning);

		//Debug partitioning
		/*for (int i = 0; i<_nProcessors; i++) {
			for (int j = grid.vdist[i]; j< grid.vdist[i+1]; j++) grid.cellsPartitioning[j] = i;
		};		*/
		//for (int i = 10; i<51; i++) grid.cellsPartitioning[i] = 0;

		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good.");	

		//Otput result information		
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Partitioning edgecut = ", edgecut);		
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCells = ", grid.nCells);		
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nProperCells = ", grid.nProperCells);		

		grid.GenerateLocalCells(_parallelHelper.getRank(), grid.cellsPartitioning);

		//Otput result
		msg.str(std::string());
		msg<<"Number of local cells = "<<grid.nCellsLocal;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());

		_parallelHelper.part = grid.cellsPartitioning;

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Grid partitioning finished");	

		return TURBO_OK;
	}; 

	//Load required geometric grid data to appropriate data structures (provided we have partioned grid)
	turbo_errt GenerateGridGeometry(Grid& grid) {
		//return turbo_errt::TURBO_OK;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Geometry generation started");

		grid.GenerateLocalCells(_parallelHelper.getRank(), grid.cellsPartitioning);
		grid.GenerateLocalFaces(_parallelHelper.getRank());
		grid.UpdateGeometricProperties();

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCellsLocal = ", grid.nCellsLocal);
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCellsLocal + dummyLocal = ", grid.localCellIndexes.size());		
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating local geometry and faces");	
		return TURBO_OK;
	};

	turbo_errt InitParallelExchange() {		
		return TURBO_OK; //TO DO remove
		//Determine values required (TO DO Generalize) closest neighbours for now
		std::set<int> requiredValuesCells;
		std::set<int> requiredGradientsCells;
		for (Cell* cell : _gridPtr->localCells) {
			for (int c : cell->NeigbourCells) {
				if (_parallelHelper.getRank() != _gridPtr->cellsPartitioning[c]) {
					requiredValuesCells.insert(c);
					requiredGradientsCells.insert(c);				
				};
			};
		};

		//Set partitioning info
		_parallelHelper.SetCellsPartitioning(_gridPtr->cellsPartitioning);
		
		//Make requests	
		_parallelHelper.ClearRequests();
		
		//Debug info
		std::stringstream msg;
		msg.clear();
		msg<<"required cells :";
		for (int p : requiredValuesCells) msg<<" "<<p;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Cell values requests
		for (int cellID : requiredValuesCells) {
			_parallelHelper.RequestValues(cellID);
		};

		//Debug
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Request formation finished");	

		//Init exchange
		//_parallelHelper.InitExchange(_gasModels);

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Parallel data exchange initialized");	
		return TURBO_OK;
	};

	turbo_errt ParallelExchangeValues() {	
		return TURBO_OK; //TO DO remove

		//Exchange values
		//_parallelHelper.ExchangeValues(_grid, Values);
		
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Started : Exchange values.");	

		//Debug info
		std::stringstream msg;
		msg.clear();
		msg<<"toSendValuesNumberByProc :";
		for (int p : _parallelHelper.toSendValuesNumberByProc) msg<<"("<<p<<") ";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Allocate memory and fill datastructures
		int s = 0;
		int r = 0;
		std::vector<int> sdispl;
		std::vector<int> rdispl;
		std::vector<double> recvbuf;
		std::vector<double> sendbuf;		
		std::vector<int> sendcounts;
		std::vector<int> recvcounts;
		for (int i = 0; i<_nProcessors; i++) {
			sdispl.push_back(s);
			rdispl.push_back(r);
			//From current proc to i-th proc send values that proc requested			
			for (int j = 0; j<_parallelHelper.toSendValuesNumberByProc[i]; j++) {
				//Cell index
				int cellGlobalIndex = _parallelHelper.toSendValuesByProc[i][j];
				//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "cellGlobalIndex = ", cellGlobalIndex);	
				int cellLocalIndex = _gridPtr->cellsGlobalToLocal[cellGlobalIndex];
				//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "cellLocalIndex = ", cellLocalIndex);	
				//Write values from cell to send buffer
				for (int k = 0; k < nVariables; k++) {
					sendbuf.push_back(Values[cellLocalIndex * nVariables + k]);
				};
			};
			sendcounts.push_back(_parallelHelper.toSendValuesNumberByProc[i] * nVariables);
			recvcounts.push_back(_parallelHelper.toRecvValuesNumberByProc[i] * nVariables);			
			s += _parallelHelper.toSendValuesNumberByProc[i] * nVariables;
			r += _parallelHelper.toRecvValuesNumberByProc[i] * nVariables;			
		};

		//Debug		
		/*msg.clear();
		msg.clear();
		msg<<"Sendbuf values : ";
		for (double p : sendbuf) {
			msg<<p<<" ";
		};*/
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Info : Exchange values 1.");	

		//Make sure that sendbuf isnt empty
		sendbuf.push_back(0);

		//Allocate recieve buffer
		recvbuf.resize(r + 1);		

		//Exchange values
		MPI_Alltoallv(&sendbuf[0], &sendcounts[0], &sdispl[0], MPI_LONG_DOUBLE, &recvbuf[0], &recvcounts[0], &rdispl[0], MPI_LONG_DOUBLE, _parallelHelper.getComm());		

		//Debug		
	/*	msg.clear();
		msg<<"Recvbuf values : ";
		for (double p : recvbuf) {
			msg<<p<<" ";
		};
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	*/

		//Write recieved values to appropriate data structure
		int pointer = 0;
		for (int i = 0; i<_nProcessors; i++) {								
			for (int j = 0; j<_parallelHelper.toRecvValuesNumberByProc[i]; j++) {			
				//Cell index
				int cellGlobalIndex = _parallelHelper.toRecvValuesByProc[i][j];
				//Store values
				_parallelHelper.RequestedValues[cellGlobalIndex] = std::vector<double>(nVariables, 0);
				for (int k = 0; k < nVariables; k++) {
					_parallelHelper.RequestedValues[cellGlobalIndex][k] = recvbuf[pointer+k];
				};				
				//Increment pointer
				pointer += nVariables;
			};			
		};			

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Finished : Exchange values.");	
		return TURBO_OK;
	};

	turbo_errt ParallelExchangeGradients() {
		return TURBO_OK; //TO DO remove
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Started : Exchange gradients.");	

		//Debug info
		std::stringstream msg;
		msg.clear();
		msg<<"toSendValuesNumberByProc :";
		for (int p : _parallelHelper.toSendValuesNumberByProc) msg<<"("<<p<<") ";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Allocate memory and fill datastructures
		int s = 0;
		int r = 0;
		std::vector<int> sdispl;
		std::vector<int> rdispl;
		std::vector<double> recvbuf;
		std::vector<double> sendbuf;		
		std::vector<int> sendcounts;
		std::vector<int> recvcounts;
		for (int i = 0; i<_nProcessors; i++) {
			sdispl.push_back(s);
			rdispl.push_back(r);
			//From current proc to i-th proc send values that proc requested			
			for (int j = 0; j<_parallelHelper.toSendValuesNumberByProc[i]; j++) {
				//Cell index
				int cellGlobalIndex = _parallelHelper.toSendValuesByProc[i][j];
				_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "cellGlobalIndex = ", cellGlobalIndex);	
				int cellLocalIndex = _gridPtr->cellsGlobalToLocal[cellGlobalIndex];
				_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "cellLocalIndex = ", cellLocalIndex);	
				//Write values from cell to send buffer
				for (int k = 0; k < nVariables; k++) {
					sendbuf.push_back(Values[cellLocalIndex * nVariables + k]);
				};
			};
			sendcounts.push_back(_parallelHelper.toSendValuesNumberByProc[i] * nVariables);
			recvcounts.push_back(_parallelHelper.toRecvValuesNumberByProc[i] * nVariables);			
			s += _parallelHelper.toSendValuesNumberByProc[i] * nVariables;
			r += _parallelHelper.toRecvValuesNumberByProc[i] * nVariables;			
		};

		//Debug		
		msg.clear();
		msg.clear();
		msg<<"Sendbuf values : ";
		for (double p : sendbuf) {
			msg<<p<<" ";
		};
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Info : Exchange values 1.");	

		//Make sure that sendbuf isnt empty
		sendbuf.push_back(0);

		//Allocate recieve buffer
		recvbuf.resize(r + 1);		

		//Exchange values
		MPI_Alltoallv(&sendbuf[0], &sendcounts[0], &sdispl[0], MPI_LONG_DOUBLE, &recvbuf[0], &recvcounts[0], &rdispl[0], MPI_LONG_DOUBLE, _parallelHelper.getComm());		

		//Debug		
		msg.clear();
		msg<<"Recvbuf values : ";
		for (double p : recvbuf) {
			msg<<p<<" ";
		};
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Write recieved values to appropriate data structure
		int pointer = 0;
		for (int i = 0; i<_nProcessors; i++) {								
			for (int j = 0; j<_parallelHelper.toRecvValuesNumberByProc[i]; j++) {			
				//Cell index
				int cellGlobalIndex = _parallelHelper.toRecvValuesByProc[i][j];
				//Store values
				_parallelHelper.RequestedValues[cellGlobalIndex] = std::vector<double>(nVariables, 0);
				for (int k = 0; k < nVariables; k++) {
					_parallelHelper.RequestedValues[cellGlobalIndex][k] = recvbuf[pointer+k];
				};				
				//Increment pointer
				pointer += nVariables;
			};			
		};			

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Finished : Exchange gradients.");	
		return TURBO_OK;
	};

	turbo_errt BindConfiguration(Configuration& configuration) {
		_configuration = configuration;

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished : configuration loaded programmatically");	

		return turbo_errt::TURBO_OK;
	};

	turbo_errt ReadConfiguration(std::string fname) {
		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished : configuration loaded from file");	
		return TURBO_OK;
	};

	turbo_errt InitBoundaryConditions() {
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Initializing boundary conditions");	

		//Create Boundary conditions from configuration
		bool isUnspecifiedBC = false;
		for (const auto& bc : _gridPtr->patchesNames) {
			const std::string& bcName = bc.first;			
			int bcMarker = bc.second;			
			if ( _configuration.BoundaryConditions.find(bcName) == _configuration.BoundaryConditions.end()) {
				//Boundary condition unspecified
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unspecified boundary condition :" + bcName);	
				isUnspecifiedBC = true;	

				//Skip missing boundary condition
				continue;
			};
			
			
			//Initialize boundary conditions
			bool bcTypeCheckPassed = false;
			BoundaryConditionConfiguration& bcConfig = _configuration.BoundaryConditions[bcName];

			std::string mName = bcConfig.MaterialName;
			if (_configuration.GasModelNameToIndex.find(mName) == std::end(_configuration.GasModelNameToIndex)) {
				//Material name unspecified
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unspecified material name " + mName + ", for boundary condition :" + bcName);	

				//Synchronize
				_parallelHelper.Barrier();
				return turbo_errt::TURBO_ERROR;
			};
			int materialIndex = _configuration.GasModelNameToIndex[mName];

			BCType_t bcType = bcConfig.BoundaryConditionType;
			if (bcType == BCType_t::BCSymmetryPlane) {				
				_boundaryConditions[bcMarker] = std::unique_ptr<BoundaryConditions::BoundaryCondition>(new BoundaryConditions::BCSymmetryPlane());
				bcTypeCheckPassed = true;
			};
			if (bcType == BCType_t::BCOutflowSupersonic) {
				_boundaryConditions[bcMarker] = std::unique_ptr<BoundaryConditions::BoundaryCondition>(new BoundaryConditions::BCNatural());
				bcTypeCheckPassed = true;
			};
			/*if (bcType == BCType_t::BCInflowSupersonic) {
				_boundaryConditions[bcMarker] = std::unique_ptr<BoundaryConditions::BoundaryCondition>(new BoundaryConditions::BCInflowSupersonic());
				bcTypeCheckPassed = true;
			};*/
			if (bcType == BCType_t::BCGeneral) {
				_boundaryConditions[bcMarker] = std::unique_ptr<BoundaryConditions::BoundaryCondition>(new BoundaryConditions::BCGeneral());
				bcTypeCheckPassed = true;
			};

			if (!bcTypeCheckPassed) {
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "BCType is not supported.");

				//Synchronize
				_parallelHelper.Barrier();
				return turbo_errt::TURBO_ERROR;			
			} else {
				//Attach needed data structures to boundary condition class
				_boundaryConditions[bcMarker]->setGrid(_gridPtr);
				_boundaryConditions[bcMarker]->setGasModel(_gasModels);
				_boundaryConditions[bcMarker]->setMaterialIndex(materialIndex);

				//Load specific configuration parameters
				_boundaryConditions[bcMarker]->loadConfiguration(bcConfig);

				//Set and check mesh movement type (TO DO check for consistency with ALE settings)
				_boundaryConditions[bcMarker]->movementType = bcConfig.MovementType;
			};
		};

		if (isUnspecifiedBC) {
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unspecified boundary condition.");

			//Synchronize
			_parallelHelper.Barrier();
			return turbo_errt::TURBO_ERROR;
		};

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished initialization of boundary conditions");	
		return turbo_errt::TURBO_OK;
	};

	turbo_errt ReadInitialConditions(std::string solutionName) {		
		//Read physical data
		//std::vector<double> density(_gridPtr->nCellsLocal);
		//std::vector<double> velocityX(_gridPtr->nCellsLocal);
		//std::vector<double> velocityY(_gridPtr->nCellsLocal);
		//std::vector<double> velocityZ(_gridPtr->nCellsLocal);		
		//_cgnsReader.ReadField(_gridPtr, solutionName, "Density", density);
		//_cgnsReader.ReadField(_gridPtr, solutionName, "VelocityX", velocityX);
		//_cgnsReader.ReadField(_gridPtr, solutionName, "VelocityY", velocityY);
		//_cgnsReader.ReadField(_gridPtr, solutionName, "VelocityZ", velocityZ);
		//_cgnsReader.ReadField(_gridPtr, solutionName, "VelocityZ", velocityZ);

		////Fill data structures
		//nVariables = 5;
		//Values.resize(nVariables * _gridPtr->nCellsLocal);
		//for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
		//	double ro = density[i];
		//	double rou = velocityX[i] * ro;
		//	double rov = velocityY[i] * ro;
		//	double row = velocityZ[i] * ro;
		//	Values[i + 0] = ro;
		//	Values[i + 1] = rou;
		//	Values[i + 2] = rov;
		//	Values[i + 3] = row;
		//};

		////Synchronize		
		//_parallelHelper.Barrier();
		//_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished reading initial conditions");	
		return TURBO_OK;
	};

	turbo_errt GenerateInitialConditions(InitialConditions::InitialConditions* initConditions) {
		//Attach data structures
		initConditions->setGrid(_gridPtr);
		initConditions->setGasModel(_gasModels);
				
		//For each local cell write initial conditions
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {			
			Cell* c = _gridPtr->localCells[i];

			//Get initial values
			std::vector<double> values = initConditions->getInitialValues(*c);						
			for (int j = 0; j<nVariables; j++) {				
				Values[i * nVariables + j] = values[j];
			};

			//Get material distribution
			CellGasModel[i] = initConditions->getInitialGasModelIndex(*c);			
		};

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating initial conditions");	
		return TURBO_OK;
	};

	turbo_errt FillInitialConditions(std::function<std::vector<double>(Vector)> FVariables, std::function<int(Vector)> FMaterial) {
		//For each local cell write initial conditions
		for (int i = 0; i < _gridPtr->nCellsLocal; i++) {
			Cell* cell = _gridPtr->localCells[i];

			//Get initial values
			std::vector<double> values = FVariables(cell->CellCenter);						
			for (int j = 0; j<nVariables; j++) {				
				Values[i * nVariables + j] = values[j];
			};

			//Get material distribution
			CellGasModel[i] = FMaterial(cell->CellCenter);
		};

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating initial conditions");	
		return TURBO_OK;
	};

	turbo_errt SaveSolution(std::string fname, std::string solutionName) {	
		//Open file
		_cgnsWriter.OpenFile(fname);
		int nv = nVariables;
		std::vector<std::string> storedFields = _gasModels[0]->GetStoredFieldsNames();

		//Write simulation type

		//Write iterative data
		double iter = stepInfo.Iteration;
		double time = stepInfo.Time;		
		if (_simulationType == TimeAccurate) {
			_cgnsWriter.WriteIterativeData(time, iter);
		};

		//Write solution node		
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "1");	
		_cgnsWriter.WriteSolution(_gridPtr, solutionName);
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "2");	

		//Write physical quantities
		std::vector<double> buffer(_gridPtr->nCellsLocal);
		std::vector<double> buffer2(_gridPtr->nCellsLocal);
		std::vector<double> buffer3(_gridPtr->nCellsLocal);
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "3");	

		//Stored fields
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "storedFields.size() = ", storedFields.size());	
		for (int fieldIndex = 0; fieldIndex < storedFields.size(); fieldIndex++ ) {
			std::string fieldName = storedFields[fieldIndex];			
			//Write to buffer
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Field name " + fieldName);	
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCellsLocal = ", _gridPtr->nCellsLocal);				
			for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
				int nmat = GetCellGasModelIndex(_gridPtr->localCells[i]->GlobalIndex);
				GasModel::ConservativeVariables U(&Values[i * nv]);
				std::vector<double> storedValues = _gasModels[nmat]->GetStoredValuesFromConservative(U);				
				buffer[i] = storedValues[fieldIndex];
			};
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "2");	
			_cgnsWriter.WriteField(_gridPtr, solutionName, fieldName, buffer);		
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "3");	
		};

		//Additional fields (TO DO) generalize
		//VelocityX
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double rou = Values[i * nv + 1];
			double u = rou / ro;
			buffer[i] = u;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "VelocityX", buffer); 

		//VelocityY
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double rov = Values[i * nv + 2];
			double v = rov / ro;
			buffer[i] = v;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "VelocityY", buffer); 

		//VelocityZ
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double row = Values[i * nv + 3];
			double w = row / ro;
			buffer[i] = w;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "VelocityZ", buffer); 

		//Pressure and sound speed and Gruneisen coef
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			int nmat = CellGasModel[i];
			double P = 0;
			double C = 0;
			double Gr = 0;
			_gasModels[nmat]->GetPressureAndSoundSpeed(&Values[i * nv + 0], P, C, Gr);
			buffer[i] = P;
			buffer2[i] = C;
			buffer3[i] = Gr;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Pressure", buffer);
		_cgnsWriter.WriteField(_gridPtr, solutionName, "SoundSpeed", buffer2);
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Gruneisen", buffer3);

		//Internal energy
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double u = Values[i * nv + 1] / ro;
			double v = Values[i * nv + 2] / ro;
			double w = Values[i * nv + 3] / ro;
			double E = Values[i * nv + 4] / ro;
			double e = E - (u*u + v*v + w*w) / 2.0;	
			buffer[i] = e;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "EnergyInternal", buffer);

		//Temperature
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			int nmat = CellGasModel[i];
			double T = _gasModels[nmat]->GetTemperature(&Values[i * nv + 0]);
			buffer[i] = T;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Temperature", buffer);

		//Phase information
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			int nmat = CellGasModel[i];
			GasModel::MediumPhase phase = _gasModels[nmat]->GetPhase(&Values[i * nv + 0]);
			int indicator = phase == GasModel::MediumPhase::AboveMeltingPoint;
			buffer[i] = indicator;
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Phase", buffer);

		//Material index
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			buffer[i] = GetCellGasModelIndex(i);
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Material", buffer);

		//Quality
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			buffer[i] = MeshQuality::Anisotropy(_gridPtr, *_gridPtr->localCells[i]);
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Anisotropy", buffer);
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			buffer[i] = MeshQuality::LinearSize(_gridPtr, *_gridPtr->localCells[i]);
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "LinearSize", buffer);
		for (int i = 0; i<_gridPtr->nCellsLocal; i++) {
			buffer[i] = MeshQuality::LinearSizeRatio(_gridPtr, *_gridPtr->localCells[i]);
		};
		_cgnsWriter.WriteField(_gridPtr, solutionName, "LinearSizeRatio", buffer);

		//Write partitioning to solution
		std::vector<double> part(_gridPtr->nCellsLocal, _parallelHelper.getRank());		
		_cgnsWriter.WriteField(_gridPtr, solutionName, "Processor", part); 

		//Close file
		_cgnsWriter.CloseFile();

		return TURBO_OK;
	};	

	//Reconstruct solution for each face, left state and right state
	void ComputeSolutionReconstruction(const std::vector<double>& cellValues, std::vector<std::vector<double> >& cellValuesL, std::vector<std::vector<double> >& cellValuesR) {
		int nv = nVariables;
		int nFaces = _gridPtr->localFaces.size();
		if (cellValuesL.size() != nFaces) cellValuesL.resize(nFaces);
		if (cellValuesR.size() != nFaces) cellValuesR.resize(nFaces);

		//In every cell proper cell reconstruct solution
		for (int localCellIndex = 0; localCellIndex < _gridPtr->nCellsLocal; localCellIndex++) {
			Cell* cell = _gridPtr->localCells[localCellIndex];			  			
			CellSpatialDiscretisation cellSpatial(cell->GlobalIndex, _gridPtr, _spatialDiscretisation);			

			//Compute stencil
			int rootMaterial = GetCellGasModelIndex(cell->GlobalIndex);

			std::function<bool(Cell&)> isGoodCell = [&](Cell& c) {
				//return (!c.IsDummy);
				return (!c.IsDummy) && (GetCellGasModelIndex(c.GlobalIndex) == rootMaterial);
			};

			cellSpatial.CalculateStencil(isGoodCell);

			//Gather cell values according to stencil
			std::vector< int > stencilMaterials;
			std::vector< std::vector<double> > stencilValues;				
			for (int index : cellSpatial.stencilIndexes()) {	
				stencilMaterials.push_back( GetCellGasModelIndex(index) );
				stencilValues.push_back(std::vector<double>(cellValues.begin() + nv * index, cellValues.begin() + nv * (index + 1))); 
			};										

			//Reconstruct solution
			cellSpatial.ReconstructSolution(stencilMaterials, stencilValues);

			//Interpolate to every face
			for (int faceIndex : cell->Faces) {
				Face& face = _gridPtr->localFaces[faceIndex];
				std::vector<double> U = cellSpatial.GetSolutionAtFace(face);
				
				if (face.FaceCell_1 == localCellIndex) {
					cellValuesL[face.GlobalIndex] = U;
				} else {
					cellValuesR[face.GlobalIndex] = U;
				};
			};
		};

		//Apply boundary conditions for dummy cells
		for (int localCellIndex = _gridPtr->nCellsLocal; localCellIndex < _gridPtr->nCells; localCellIndex++) {
			//Get cell and face info
			Cell* cell = _gridPtr->localCells[localCellIndex];			  			
			Face& face = _gridPtr->localFaces[cell->Faces[0]];
			int nmat = GetCellGasModelIndex(cell->GlobalIndex);
			int faceLocalIndex = face.GlobalIndex;			

			//Boundary face material correction			
			nmat = GetCellGasModelIndex(face.FaceCell_1); 

			//Obtain dummy cell state
			cellValuesR[faceLocalIndex] = _boundaryConditions[cell->BCMarker]->getDummyValues(nmat, cellValuesL[faceLocalIndex], *cell);
		};

		//Handle periodic connectivity issue
		assert(_nProcessors == 1);
		for (Face& face : _gridPtr->localFaces) {
			//Additional check
			if (cellValuesL[face.GlobalIndex].size() == 0) throw Exception("Reconstruction incomplete");

			//If right state is empty check for periodic counterpart
			if (cellValuesR[face.GlobalIndex].size() == 0) {
				//Find corresponding nodes
				std::set<int> identicalNodes;
				for (int faceNode : face.FaceNodes) {
					identicalNodes.insert(faceNode);
					identicalNodes.insert(std::begin(_gridPtr->periodicNodesIdentityList[faceNode]), std::end(_gridPtr->periodicNodesIdentityList[faceNode]));
				};

				//Find corresponding face on right neighbour cell
				int facePairIndex = -1;
				bool isFound = false;
				Cell& pairCell = _gridPtr->Cells[face.FaceCell_2];
				for (int faceIndex : pairCell.Faces) {
					Face& pairFace = _gridPtr->localFaces[faceIndex];
					isFound = true;

					//Check if all face nodes are right
					for (int nodeIndex : pairFace.FaceNodes) {
						if (identicalNodes.find(nodeIndex) == identicalNodes.end()) {
							isFound = false;
							break;
						};
					};

					//Break if we found right face
					if (isFound) {
						facePairIndex = faceIndex;
						break;
					};
				};

				if (isFound && (cellValuesL[facePairIndex].size() != 0)) {
					//Everything went fine
					cellValuesR[face.GlobalIndex] = cellValuesL[facePairIndex];
				} else {
					throw Exception("Reconstruction incomplete");
				};
			};
		};

	}; //Compute solution reconstruction

	//Compute convective flux and max wave propagation speed throught each face
	void ComputeConvectiveFluxes(std::vector<std::vector<double>>& fluxes, std::vector<double>& maxWaveSpeed, std::vector<std::vector<double> >& cellValuesL, std::vector<std::vector<double> >& cellValuesR, std::vector<double>& ALEindicators) {
		timerConvective.Resume();

		//Nullify all fluxes
		for (std::vector<double>& flux : fluxes) {
			for (double& v : flux) v = 0;
		};					

		//Exchange requested cell values
		ParallelExchangeValues();
		std::stringstream msg;
		/*msg.clear();		
		msg<<"Obtained values : ";
		for (std::pair<int, std::vector<double>> p : _parallelHelper.RequestedValues) {
			msg<<"cell="<<p.first<<"values=("<<p.second[0]<<","<<p.second[1]<<","
				<<p.second[2]<<","<<p.second[3]<<","<<p.second[4]<<")";
		};
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());*/	


		//Compute convective flux for each cell face and apply boundary conditions								
		#pragma omp for
		for (int i = 0; i<_gridPtr->localFaces.size(); i++) {
			Face& f = _gridPtr->localFaces[i];											
			std::vector<double> flux;
			int ALEIndicator = ALEindicators[i];

			//Solution reconstruction
			int nmatL = GetCellGasModelIndex(f.FaceCell_1);
			int nmatR = GetCellGasModelIndex(f.FaceCell_2);
			GasModel::ConservativeVariables UL = cellValuesL[f.GlobalIndex];
			GasModel::ConservativeVariables UR = cellValuesR[f.GlobalIndex];

			//Determine face velocity
			Vector faceVelocityAverage = Vector(0,0,0);
			if (_ALEmethod.ALEMotionType != ALEMethod::ALEMotionType::PureEulerian) {			
				faceVelocityAverage = _ALEmethod.facesVelocity[f.GlobalIndex];
			};

			//Compute flux		
			RiemannProblemSolutionResult result = rSolver->Solve(nmatL, UL, nmatR, UR, f, faceVelocityAverage);			

			//Store interface velocity and pressure for ALE
			_ALEmethod.facesPressure[f.GlobalIndex] = FacePressure[f.GlobalIndex] = result.Pressure;			

			//Store wave speeds			
			maxWaveSpeed[f.GlobalIndex] = result.MaxEigenvalue;

			//Compute flux exactly for Lagrangian interfaces				
			if (ALEIndicator == 1.0) {
				//if (IsBoundaryFace(f)) result.Pressure = 1e5;
				result.Fluxes[0] = 0;
				result.Fluxes[1] = f.FaceNormal.x * result.Pressure;
				result.Fluxes[2] = f.FaceNormal.y * result.Pressure;
				result.Fluxes[3] = f.FaceNormal.z * result.Pressure;
				result.Fluxes[4] = (faceVelocityAverage * f.FaceNormal) * result.Pressure;
			};

			//Store flux
			fluxes[f.GlobalIndex] = result.Fluxes;								

		};

		//Synchronyze
		_parallelHelper.Barrier();
		timerConvective.Pause();
	};	

	//Compute residual for each cell		
	void ComputeResidual(std::vector<double>& residual, std::vector<double>& cellValues) {	
		//Reconstruct solution at faces
		std::vector<std::vector<double> > cellValuesL;
		std::vector<std::vector<double> > cellValuesR;
		ComputeSolutionReconstruction(cellValues, cellValuesL, cellValuesR); 

		//Compute ALE indicator values
		timerALE.Resume();
		std::vector<double> ALEindicators(_gridPtr->nFaces, 0);
		for (Face& f : _gridPtr->localFaces) {
			int nmatL = GetCellGasModelIndex(f.FaceCell_1);
			int nmatR = GetCellGasModelIndex(f.FaceCell_2);
			if (!IsDummyCell(f.FaceCell_2)) {
				//Determine ALE movement indicators for inner cells
				if (_ALEmethod.ALEMotionType == ALEMethod::ALEMotionType::PureLagrangian) {
					ALEindicators[f.GlobalIndex] = 1; // all move
				};
				if (_ALEmethod.ALEMotionType == ALEMethod::ALEMotionType::ALEMaterialInterfaces) { 
					ALEindicators[f.GlobalIndex] = 0;
					if (nmatL != nmatR) { //move material interface
						ALEindicators[f.GlobalIndex] = 1;
					};				
				};
				if (_ALEmethod.ALEMotionType == ALEMethod::ALEMotionType::PureEulerian) {
					ALEindicators[f.GlobalIndex] = 0;
				};
			} else {
				//Determine ALE movement indicators for boundaries
				BoundaryConditionMovementType movementType = GetCellBoundaryCondition(f.FaceCell_2)->movementType;
				if (movementType == BoundaryConditionMovementType::Fixed) {
					ALEindicators[f.GlobalIndex] = 0;
				};
				if (movementType == BoundaryConditionMovementType::FreeSurface) {
					ALEindicators[f.GlobalIndex] = 1;
				};
			};
		};

		//ALE step determine mesh motion
		if (_ALEmethod.ALEMotionType != ALEMethod::ALEMotionType::PureEulerian) {
			//Compute face velocities
			_ALEmethod.facesVelocity.clear();
			for (Face& f : _gridPtr->localFaces) {
				int ALEindicator = ALEindicators[f.GlobalIndex];
				if (ALEindicator == 0) continue; //Skip face that don't participate in obligatory motion
				int nmatL = GetCellGasModelIndex(f.FaceCell_1);
				int nmatR = GetCellGasModelIndex(f.FaceCell_2);
				if (IsBoundaryFace(f)) {
					nmatR = nmatL;
				};
				GasModel::ConservativeVariables UL = GasModel::ConservativeVariables(GetCellValues(f.FaceCell_1));
				GasModel::ConservativeVariables UR = GasModel::ConservativeVariables(GetCellValues(f.FaceCell_2));															
				Vector faceVelocity = _ALEmethod.ComputeFaceVelocity(_gasModels[nmatL], UL, _gasModels[nmatR], UR, f, ALEindicator);
				_ALEmethod.facesVelocity[f.GlobalIndex] = faceVelocity;
			};

			//Compute moving nodes velocities
			_ALEmethod.ComputeMovingNodesVelocities();

			//Impose boundary movement
			for (Face& f : _gridPtr->localFaces) {
				if (!_gridPtr->IsBoundaryFace(f)) continue;
				BoundaryConditionMovementType movementType = GetFaceBoundaryCondition(f)->movementType;
				if (movementType == BoundaryConditionMovementType::Fixed) {
					//Each node must stick to its place
					for (int nodeIndex : f.FaceNodes) _ALEmethod.nodesVelocity[nodeIndex] = Vector(0,0,0);
				};
			};

			//Compute velocities of all other nodes
			_ALEmethod.ComputeFreeNodesVelocities();

			//Compute new face velocities
			for (Face& f : _gridPtr->localFaces) {
				_ALEmethod.facesVelocity[f.GlobalIndex] = _ALEmethod.ComputeFaceVelocityByNodes(f);
			};
		};
		timerALE.Pause();		

		//Compute convective fluxes and max wave speeds		
		ComputeConvectiveFluxes(FaceFluxes, MaxWaveSpeed, cellValuesL, cellValuesR, ALEindicators);
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Convective fluxes calculated");

		////Compute gradients
		//if (IsGradientRequired) ComputeGradients();

		////Compute viscous fluxes
		//ComputeViscousFluxes();	
		
		//Distribute residual to each cell		
		for ( int cellIndex = 0; cellIndex<_gridPtr->nCellsLocal; cellIndex++ )
		{
			//Get cell reference
			Cell& cell = GetLocalCell(cellIndex);

			//Nullify residual first
			for (int i = 0; i < nVariables; i++) residual[cellIndex*nVariables + i] = 0;

			//For each face distribute
			double sweptVolumeTotal = 0;
			std::vector<int>& nFaces = cell.Faces;
			for (int nFaceIndex : nFaces)
			{
				//Get face reference
				Face& face = GetLocalFace(nFaceIndex);
				
				//Determine direction
				int fluxDirection = (face.FaceCell_1 == cell.GlobalIndex) ? 1 : -1;		

				//Distribute face flux to cell
				std::vector<double> fluxc = FaceFluxes[nFaceIndex];
				for (int j = 0; j<nVariables; j++) {
					residual[cellIndex * nVariables + j] +=  (fluxc[j]) * face.FaceSquare * fluxDirection;
				};		

				//Mesh deformation contribution (GCL part statisfied explicitly)
				if (_ALEmethod.ALEMotionType != ALEMethod::ALEMotionType::PureEulerian) {
					//Compute volume swept by current face
					double sweptVolume = _ALEmethod.facesVelocity[face.GlobalIndex] * face.FaceNormal * face.FaceSquare * fluxDirection;
					sweptVolumeTotal += sweptVolume;

					//Obtain cell values at moving face
					std::vector<double>& values =  (face.FaceCell_1 == cell.GlobalIndex) ? cellValuesL[cellIndex] : cellValuesR[cellIndex];		
				};
			};

			//Mesh deformation contribution (GCL part statisfied explicitly)
			if (_ALEmethod.ALEMotionType != ALEMethod::ALEMotionType::PureEulerian) {
				for (int j = 0; j<nVariables; j++) {
					residual[cellIndex * nVariables + j] += cellValues[cellIndex * nVariables + j] * sweptVolumeTotal;
				};			
			}
		};

		

		//Compute source term
		for ( int cellIndex = 0; cellIndex<_gridPtr->nCellsLocal; cellIndex++ )
		{
			Cell* cell = _gridPtr->localCells[cellIndex];						
			std::vector<double> val(&cellValues[cellIndex * nVariables], &cellValues[cellIndex * nVariables] + nVariables);
			std::vector<double> RSource = _sources->GetResidual(val);
			for (int i = 0; i < nVariables; i++) residual[cellIndex*nVariables + i] -= cell->CellVolume * RSource[i];
		};

 		return;
	};

	//Compute spectral radii estimate for each cell
	void ComputeSpectralRadius(std::map<int, double>& spectralRadius, std::vector<double>& maxWaveSpeed, std::vector<double>& cellValues) {
		spectralRadius.clear();		
		for (int cellIndex = 0; cellIndex<_gridPtr->nCellsLocal; cellIndex++)
		{
			Cell* cell = _gridPtr->localCells[cellIndex];			
			spectralRadius[cellIndex] = 0;
			std::vector<int>& nFaces = cell->Faces;
			for (int nFaceIndex : nFaces)
			{
				//Blazek f. 6.21
				Face& face = _gridPtr->localFaces[nFaceIndex];			
				spectralRadius[cellIndex] +=  maxWaveSpeed[nFaceIndex] * face.FaceSquare;
			};
		};
	};

	//Compute local time step for each cell utilizing spectral radius estimates
	void ComputeLocalTimeStep(std::map<int, double>& localTimeStep, std::map<int, double>& spectralRadius) {
		localTimeStep.clear();
		for (int cellIndex = 0; cellIndex<_gridPtr->nCellsLocal; cellIndex++)
		{
			Cell* cell = _gridPtr->localCells[cellIndex];			
			double sR = spectralRadius[cellIndex];
			localTimeStep[cellIndex] = CFL * cell->CellVolume / sR; //Blazek f. 6.20
		}
	};	

	//Explicit time step
	void ExplicitTimeStep() {
		std::stringstream msg;

		//Compute residual	
		timerResidual.Resume();
		ComputeResidual(Residual, Values);
		timerResidual.Pause();

		//Synchronize
		_parallelHelper.Barrier();

		//Determine time step as global minimum over local time steps		
		std::map<int, double> spectralRadius;
		ComputeSpectralRadius(spectralRadius, MaxWaveSpeed, Values);

		//Synchronize
		_parallelHelper.Barrier();

		std::map<int, double> localTimeStep;
		ComputeLocalTimeStep(localTimeStep, spectralRadius);

		//Synchronize
		_parallelHelper.Barrier();

		stepInfo.TimeStep = MaxTime - stepInfo.Time; 
		if (SaveSolutionSnapshotTime > 0) stepInfo.TimeStep = min(NextSnapshotTime - stepInfo.Time, stepInfo.TimeStep);
		for (std::pair<int, double> p : localTimeStep)
		{
			double& timeStep = p.second;
			if (timeStep < stepInfo.TimeStep) stepInfo.TimeStep = timeStep;
		}
		stepInfo.TimeStep = _parallelHelper.Min(stepInfo.TimeStep);		

		//Synchronize
		_parallelHelper.Barrier();

		//Runge-Kutta explicit time stepping
		int nStages = RungeKuttaOrder;
		std::vector<double> alpha;//{ 0.0833, 0.2069, 0.4265, 1.000 };
		if (nStages == 1) {
			alpha.push_back(1.0);
		};
		if (nStages == 2) {
			alpha.push_back(0.5);
			alpha.push_back(1.0);
		};
		if (nStages == 4) {
			//Fluent coefficients second order
			/*alpha.push_back(0.25);
			alpha.push_back(0.3333);
			alpha.push_back(0.5);
			alpha.push_back(1.0);*/

			//Second order 4 stage optimized cooefficients Blazek table 6.1
			alpha.push_back(0.1084);
			alpha.push_back(0.2602);
			alpha.push_back(0.5052);
			alpha.push_back(1.0);
		};
		double sumAlpha = 0;
		for (double a : alpha) sumAlpha += a;

		stepInfo.Residual = std::vector<double>(5, 0.0);
		
		//Time stepping
		for (int stage = 0; stage<nStages; stage++) {
			double stageTimeStep = stepInfo.TimeStep * alpha[stage] / sumAlpha;

			timerALE.Resume();
			//ALE step mesh transformation
			if (_ALEmethod.ALEMotionType != ALEMethod::ALEMotionType::PureEulerian) {						
				//Move mesh
				_ALEmethod.MoveMesh(stageTimeStep);	

				//Regenerate geometric entities
				_gridPtr->UpdateGeometricProperties();
			};
			timerALE.Pause();

			for ( int cellIndex = 0; cellIndex<_gridPtr->nCellsLocal; cellIndex++ )
			{
				Cell* cell = _gridPtr->localCells[cellIndex];			
				//Update values
				for (int i = 0; i < nVariables; i++) {
					Values[cellIndex*nVariables + i] += Residual[cellIndex*nVariables + i] * (- stageTimeStep / cell->CellVolume);
					stepInfo.Residual[i] += pow(Residual[cellIndex*nVariables + i], 2);
				};
			};

			//Synchronize
			_parallelHelper.Barrier();

			if (stage == nStages - 1) break;

			//Compute new residual
			timerResidual.Resume();
			ComputeResidual(Residual, Values);
			timerResidual.Pause();

			//Synchronize
			_parallelHelper.Barrier();
		};	

		//Compute new result and residual
		//Compute RMS residual		
		for (int i = 0; i<nVariables; i++) {
			stepInfo.Residual[i] = _parallelHelper.SumDouble(stepInfo.Residual[i]);
			stepInfo.Residual[i] = sqrt(stepInfo.Residual[i]);
		};

		//Advance total time
		stepInfo.Time += stepInfo.TimeStep;
		
		//Update perfomance statistics
		stepInfo.ComputeResidualTime = timerResidual.ElapsedTimeMilliseconds();
		stepInfo.ALETime = timerALE.ElapsedTimeMilliseconds();
		stepInfo.ConvectiveFluxesTime = timerConvective.ElapsedTimeMilliseconds();

		timerALE.Reset();
		timerConvective.Reset();
		timerResidual.Reset();

		//Synchronize
		_parallelHelper.Barrier();
	};	
	
	////Compute scalar function gradient in each cell	
	//void ComputeFunctionGradient( std::vector<Vector>& grads, std::vector<double>& values, double (*func)(const std::vector<double>&) ) {
	//	//Allocate memory for gradients
	//	grads.resize(_gridPtr->nCellsLocal);		

	//	//For each cell compute gradient of given function
	//	std::vector<Vector> nPoints;
	//	std::vector<double> nValues;
	//	for (int i = 0; i<_gridPtr->nCellsLocal; i++) {			
	//		Cell& cell = *_gridPtr->localCells[i];			

	//		//Determine required set of point and values
	//		nPoints.clear();			
	//		nValues.clear();			

	//		//Add all neighbours
	//		for (int j = 0; j<cell.NeigbourCells.size(); j++) {				
	//			Cell& nCell = _gridPtr->Cells[cell.NeigbourCells[j]];
	//			std::vector<double> nU;				
	//			if (nCell.IsDummy) {
	//			};
	//			if (face.isExternal) {
	//				nU = GetDummyCellValues(values[cell.GlobalIndex], face);
	//			} else {
	//				if (face.FaceCell_1 == cell.GlobalIndex) {
	//					nU = values[face.FaceCell_2];
	//				} else {
	//					nU = values[face.FaceCell_1];
	//				};
	//			};
	//			Vector nPoint;
	//			if (face.isExternal) {
	//				nPoint = 2 * (face.FaceCenter - cell.CellCenter) + cell.CellCenter;
	//			} else {
	//				if (face.FaceCell_1 == cell.GlobalIndex) {
	//					nPoint = _gridPtr->cells[face.FaceCell_2].CellCenter;
	//				} else {
	//					nPoint = _gridPtr->cells[face.FaceCell_1].CellCenter;
	//				};
	//			};
	//			double nValue = (this->*func)(nU);
	//			nPoints.push_back(nPoint);
	//			nValues.push_back(nValue);
	//		};

	//		double cellValue =  (this->*func)(U[cell.GlobalIndex]);
	//		grads[cell.GlobalIndex] = ComputeGradientByPoints(cell.CellCenter, cellValue, nPoints, nValues);
	//	};	

	//	return;
	//};

	//Return number of local cells (stored on this process)
	inline int GetLocalCellsNumber() {
		return _gridPtr->nCellsLocal;
	};

	inline bool IsLocalCell(int globalIndex) {
		return _gridPtr->GetCellPart(globalIndex) == _parallelHelper.getRank();
	};

	inline bool IsDummyCell(int globalIndex) {
		return _gridPtr->IsDummyCell(globalIndex);
	};

	inline bool IsBoundaryFace(Face& face) {
		return _gridPtr->IsBoundaryFace(face);
	};	

	//Get face boundary condition reference
	inline BoundaryConditions::BoundaryCondition* GetFaceBoundaryCondition(Face& face) {
		if (!IsBoundaryFace(face)) throw 1; //TO DO check
		return GetCellBoundaryCondition(face.FaceCell_2);
	};

	//Get cell boundary condition reference
	inline BoundaryConditions::BoundaryCondition* GetCellBoundaryCondition(int globalIndex) {
		if (!IsDummyCell(globalIndex)) throw 1; //TO DO check
		Cell& cell = _gridPtr->Cells[globalIndex];						
		return _boundaryConditions[cell.BCMarker].get();
	};	

	//Get cell gas model index
	inline int GetCellGasModelIndex(int globalIndex, int localIndex = -1) {		
		if (IsLocalCell(globalIndex)) {
			if (IsDummyCell(globalIndex)) {
				//If dummy cell compute on the go				
				Cell& cell = _gridPtr->Cells[globalIndex];
				int nCellIndex = _gridPtr->cellsGlobalToLocal[cell.NeigbourCells[0]]; //Obtain neighbour

				//IF dummy get from boundary condition
				int nmat = _boundaryConditions[cell.BCMarker]->_nmat;
				//int nmat = CellGasModel[nCellIndex]; //Neighbour cell for now
				//nmat = 0; //Air for now
				return nmat;
			} else {
				//If proper cell return part of Values array				
				if (localIndex == -1) {
					//If local index is unknown determine it
					//Lower perfomance WARNING
					localIndex = _gridPtr->cellsGlobalToLocal[globalIndex];										
				};
				int nmat = CellGasModel[localIndex];
				return nmat;
			};
		} else {
			//Return result of interprocessor exchange
			return _parallelHelper.RequestedGasModelIndex[globalIndex];
		};
	};

	//Get cell values
	inline std::vector<double> GetCellValues(int globalIndex, int localIndex = -1) {
		/*std::stringstream msg;
		msg.str("");
		msg<<"Global Index = "<<globalIndex;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str() );*/
		if (IsLocalCell(globalIndex)) {
			if (IsDummyCell(globalIndex)) {
				//If dummy cell compute on the go
				/*msg.str("");
				msg<<"Cells = "<<globalIndex;
				_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str() );*/
				Cell& cell = _gridPtr->Cells[globalIndex];				
				//Obtain neighbour
				int nmat = GetCellGasModelIndex(cell.NeigbourCells[0]);
				//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "!" );				
				return _boundaryConditions[cell.BCMarker]->getDummyValues(nmat, Values, cell); //TO DO now boundary always 0 gas model
			} else {
				//If proper cell return part of Values array				
				if (localIndex == -1) {
					//If local index is unknown determine it
					//Lower perfomance WARNING
					localIndex = _gridPtr->cellsGlobalToLocal[globalIndex];					
					//msg.str("");					
					//msg<<"Local Index = "<<localIndex;
					//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str() );
				};
				int nmat = CellGasModel[localIndex];
				return std::vector<double>(&Values[localIndex * _gasModels[nmat]->nConservativeVariables], &Values[localIndex * _gasModels[nmat]->nConservativeVariables] + _gasModels[nmat]->nConservativeVariables);
			};
		} else {
			//Return result of interprocessor exchange
			return _parallelHelper.RequestedValues[globalIndex];
		};
	};

	//Switch verbose state
	inline void VerboseOn() { _isVerbose = true; };
	inline void VerboseOff() { _isVerbose = false; };	

	//Grid elements accessors
	inline Cell& GetLocalCell(int localIndex) {
		return *_gridPtr->localCells[localIndex];
	};

	inline Face& GetLocalFace(int localIndex) {
		return _gridPtr->localFaces[localIndex];
	};

	//Medium properties accessors
	inline double GetCellPressure(int localIndex) {
		int nmat = CellGasModel[localIndex];
		double p = _gasModels[nmat]->GetPressure(GasModel::ConservativeVariables(&Values[nVariables * localIndex]));
		return p;
	};

}; //Kernel

#endif
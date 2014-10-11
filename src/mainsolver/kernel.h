#ifndef TURBO_MAINSOLVER_KERNEL
#define TURBO_MAINSOLVER_KERNEL

#include "grid.h"
#include "mpi.h"
#include "logger.h"
#include "parmetis.h"
#include "CGNSReader.h"
#include "CGNSWriter.h"
#include "parallelHelper.h"
#include "BoundaryConditions.h"
#include "InitialConditions.h"
//#include "BoundaryCondition.h"
//#include "BCSymmetryPlane.h"
#include "configuration.h"

//Define error types
enum turbo_errt {
	TURBO_OK = 0,
	TURBO_ERROR = 1
};

//Calculation kernel
class Kernel {
private:
	//Logger
	std::string _logfilename;
	Logger _logger;

	//Helper for CGNS i\o
	GridLoading::CGNSReader _cgnsReader;
	PostProcessing::CGNSWriter _cgnsWriter;

	//Grid data
	Grid _grid;

	//Configuration
	Configuration _configuration;	

	//Gas model
	GasModel _gasModel;		

	//Solver (TO DO)		
	SimulationType_t _simulationType; //Simulation type
	double CFL;
	double RungeKuttaOrder;
	double MaxIteration;
	double MaxTime;
	double CurrentTime;
	double SaveSolutionSnapshotTime;
	int SaveSolutionSnapshotIterations;
	//Roe3DSolverPerfectGas rSolver;
	Godunov3DSolverPerfectGas rSolver;
	StepInfo stepInfo;

	//Parallel run information		
	ParallelHelper _parallelHelper;
	int _rank;
	int _nProcessors;

	//Boundary conditions
	std::map<int, BoundaryConditions::BoundaryCondition*> _boundaryConditions;

	//Internal storage		
	int nVariables; //number of variables in each cell
	std::vector<double> Values;	//Conservative variables
	std::vector<double> Residual; //Residual values	

	//Fluxes
	std::vector<std::vector<double>> FaceFluxes;
	std::vector<double> MaxWaveSpeed;

	//Gradients
	std::vector<Vector> GradientT;

public:
	//Parallel helper
	ParallelHelper* getParallelHelper() {
		return &_parallelHelper;
	};	

	//Initialize computational kernel
	turbo_errt Initilize(int *argc, char **argv[]) {		
		try {
			//Initialize parallel MPI subsystem
			_parallelHelper.Init(argc, argv);
			_rank = _parallelHelper.getRank();
			_nProcessors = _parallelHelper.getProcessorNumber();

			//Initialize loggin subsistem
			_logfilename = "kernel.log"; //TO DO
			_logger.InitLogging(_parallelHelper, _logfilename);

			//Initialize cgns i\o subsystem
			_cgnsReader.Init(_logger, _parallelHelper);
			_cgnsWriter.Init(_logger, _parallelHelper);
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
		return TURBO_OK;
	};

	turbo_errt InitCalculation() {
		//Gas model parameters
		//_gasModel.loadConfiguration(_configuration);
		nVariables = _gasModel.nConservativeVariables = 5;

		//Allocate memory for data structures
		Values.resize(nVariables * _grid.nCellsLocal);
		Residual.resize(nVariables * _grid.nCellsLocal);
		FaceFluxes.resize(_grid.nFaces);
		for (std::vector<double>& flux : FaceFluxes) {
			flux.resize(nVariables);
			for (double& v : flux) v = 0;
		};	
		MaxWaveSpeed.resize(_grid.nFaces);
		
		//Gas model setup
		_gasModel = GasModel();
		_gasModel.loadConfiguration(_configuration);

		//Solver settings
		rSolver.SetGamma(_gasModel.Gamma);
		rSolver.SetHartenEps(0.0);
		//rSolver.SetOperatingPressure(0.0);

		//Initial conditions
		InitialConditions::InitialConditions ic;
		GenerateInitialConditions(ic);

		//Boundary conditions
		InitBoundaryConditions();

		//Init parallel exchange
		InitParallelExchange();

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished calculation initialization");	
		return TURBO_OK;
	};

	turbo_errt RunCalculation() {
		//Start timer
		double workTime = 0.0;
		clock_t start, stop;
		/* Start timer */
		assert((start = clock())!=-1);

		//Calc loop
		std::stringstream msg;	
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Calculation started!");
		for (stepInfo.Iteration = 0; stepInfo.Iteration <= MaxIteration; stepInfo.Iteration++) {
			ExplicitTimeStep();

			//Output step information
			//if (_parallelHelper.IsMaster()) {
			msg.clear();
			msg.str(std::string());
			msg<<"Iteration = "<<stepInfo.Iteration<<"; Total time = "<< stepInfo.Time << "; Time step = " <<stepInfo.TimeStep << "; RMSro = "<<stepInfo.Residual[0]<<"\n";
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());
			//};

			//Debug
			//Output values
			/*msg.clear();
			msg.str(std::string());
			msg<<"Values:\n";
			for (int i = 0; i<_grid.nCellsLocal; i++) {
				Cell* c = _grid.localCells[i];				
				msg<<"ID = "<<c->GlobalIndex<<", ( ";
				for (int j = 0; j<_gasModel.nConservativeVariables; j++) msg<<Values[i * _gasModel.nConservativeVariables + j]<<" ";
				msg<<"\n";
			};
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());*/

			//Solution snapshots
			if ((SaveSolutionSnapshotIterations != 0) && (stepInfo.Iteration % SaveSolutionSnapshotIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<"dataI"<<stepInfo.Iteration<<".cgns";
				SaveGrid(snapshotFileName.str());				
				SaveSolution(snapshotFileName.str(), "Solution");
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
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Calculation finished!");	

		/* Stop timer */
		stop = clock();
		workTime = (double) (stop-start)/CLOCKS_PER_SEC;
		msg.clear();
		msg.str(std::string());
		msg<<workTime;
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Total work time = " + msg.str() + " seconds.");	
		return TURBO_OK;
	};

	turbo_errt FinalizeCalculation() {
		//Free memory				
		//Boundary conditions
		for (std::pair<int, BoundaryConditions::BoundaryCondition*> p : _boundaryConditions) {
			delete (p.second);
		};

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished calculation finalization.");	
		return TURBO_OK;
	};

	//Load grid from CGNS file
	turbo_errt LoadGrid(std::string filename) {
		_grid = _cgnsReader.LoadGrid(filename);
		PartitionGrid();
		GenerateGridGeometry();
		return TURBO_OK;
	};	

	//Save grid to CGNS file
	turbo_errt SaveGrid(std::string filename)
	{
		_cgnsWriter.CreateFile(filename);
		_cgnsWriter.WriteGridToFile(_grid);		
		_cgnsWriter.CloseFile();
		return TURBO_OK;
	}

	//Bind existing grid to kernel
	turbo_errt BindGrid(Grid& grid) {
		_grid = grid;
		PartitionGrid();
		GenerateGridGeometry();
		return TURBO_OK;
	};

	//Partion loaded grid
	turbo_errt PartitionGrid() {
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
		std::vector<idx_t> part(_grid.nCellsLocal);
		std::ostringstream msg;
		msg<<"vdist[] = \n";
		for (int i = 0; i<=_nProcessors; i++) msg<<_grid.vdist[i]<<" ";
		msg<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg<<"xadj[] = \n";
		for (int i = 0; i<_grid.xadj.size(); i++) msg<<_grid.xadj[i]<<" ";
		msg<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Call partitioning function		
		MPI_Comm _comm = _parallelHelper.getComm();
		_parallelHelper.Barrier();
		int result = ParMETIS_V3_PartKway(&_grid.vdist[0], &_grid.xadj[0], &_grid.adjncy[0], vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, &tpwgts[0], &ubvec[0], options, &edgecut, &part[0], &_comm);			
		if (result != METIS_OK) {
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "ParMETIS_V3_PartKway failed.");
			return TURBO_ERROR;
		};							

		//Gather partitioning on every processor
		std::vector<int> recvcounts(_nProcessors);
		for (int i = 0; i<_nProcessors; i++) {
			recvcounts[i] = _grid.vdist[i+1] - _grid.vdist[i];
		};
		
		/*msg.clear();
		msg<<part.size()<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg<<"recvcounts[] = \n";
		for (int i = 0; i<_nProcessors; i++) msg<<recvcounts[i]<<" ";
		msg<<"\n";

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	*/
		_parallelHelper.Allgatherv( part, recvcounts, _grid.cellsPartitioning);

		//Debug partitioning
		/*for (int i = 0; i<_nProcessors; i++) {
			for (int j = _grid.vdist[i]; j< _grid.vdist[i+1]; j++) _grid.cellsPartitioning[j] = i;
		};		*/
		//for (int i = 10; i<51; i++) _grid.cellsPartitioning[i] = 0;

		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good.");	

		//Otput result information		
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Partitioning edgecut = ", edgecut);		
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCells = ", _grid.nCells);		
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nProperCells = ", _grid.nProperCells);		
		
		//Add dummy cells
		for (int i = _grid.nProperCells; i<_grid.nCells; i++) {
			Cell* cell = &_grid.Cells[i];
			int neigbour = cell->NeigbourCells[0];
			int p = _grid.cellsPartitioning[neigbour];
			_grid.cellsPartitioning.push_back(p);			
		};

		//Otput cells partitioning
		/*msg.clear();
		msg<<"Partitioning = [ ";
		for (int p : _grid.cellsPartitioning) {
			msg<<p<<" ";
		};
		msg<<"]\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());		*/


		//Extract local cells indexes		
		int count = 0;
		_grid.localCells.clear();
		_grid.localCellIndexes.clear();
		for (int i = 0; i<_grid.nCells; i++) {
			if ((_grid.cellsPartitioning[i] == _rank) && (!_grid.Cells[i].IsDummy)) {
				_grid.localCells.push_back(&_grid.Cells[i]);
				_grid.localCellIndexes.push_back(i);				
				//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.Cells[i].getInfo());
			};
		};		
		_grid.nCellsLocal = _grid.localCellIndexes.size(); //Without dummy cells
		for (int i = 0; i<_grid.nCells; i++) {
			if ((_grid.cellsPartitioning[i] == _rank) && (_grid.Cells[i].IsDummy)) {
				_grid.localCells.push_back(&_grid.Cells[i]);
				_grid.localCellIndexes.push_back(i);
				//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.Cells[i].getInfo());
			};
		};		

		//Otput result
		msg.str(std::string());
		msg<<"Number of local cells = "<<_grid.nCellsLocal;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());															

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Grid partitioning finished");	

		return TURBO_OK;
	}; 

	//Load required geometric grid data to appropriate data structures
	//Provided we have partioned grid
	turbo_errt GenerateGridGeometry() {
		//return turbo_errt::TURBO_OK;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Geometry generation started");

		//Create local cells
		_grid.localCells.resize(_grid.localCellIndexes.size());		
		for (int i = 0; i < _grid.localCellIndexes.size(); i++) {	
			int cellGlobalIndex = _grid.localCellIndexes[i];
			int cellLocalIndex = i;
			_grid.localCells[i] = &_grid.Cells[cellGlobalIndex];	
			_grid.cellsGlobalToLocal[cellGlobalIndex] = cellLocalIndex;
		};

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCellsLocal = ", _grid.nCellsLocal);
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCellsLocal + dummyLocal = ", _grid.localCellIndexes.size());		

		//Create faces		
		int faceIndex = 0;
		std::map<std::set<int>, int> faces;	//Face nodes to face index		
		for (int i = 0; i < _grid.localCellIndexes.size(); i++) {
			//Generate all faces for the cell
			Cell* cell = _grid.localCells[i];		
			std::vector<Face> newFaces;
			if (!cell->IsDummy) {
				newFaces = _grid.ObtainFaces(cell);
			} else {
				Face boundaryFace;
				boundaryFace.CGNSType = cell->CGNSType;
				boundaryFace.FaceNodes = cell->Nodes;
				boundaryFace.BCMarker = cell->BCMarker;
				newFaces.clear();
				newFaces.push_back(boundaryFace);
			};
			for (int j = 0; j<newFaces.size(); j++) {
				Face& face = newFaces[j];				
				std::set<int> nodes(face.FaceNodes.begin(), face.FaceNodes.end());
				std::pair<std::set<int>, int> newFaceInfo(nodes, 0);
				std::pair<std::map<std::set<int>, int>::iterator, bool> result = faces.insert( newFaceInfo );				
				//Add face and save connectivity info
				if ( result.second ) {	
					face.GlobalIndex = faceIndex;
					result.first->second = face.GlobalIndex;
					_grid.localFaces.push_back(face);
					_grid.localFaces[faceIndex].isExternal = true;
					//_grid.localFaces[faceIndex].FaceCell_1 = i;
					_grid.localFaces[faceIndex].FaceCell_1 = cell->GlobalIndex;
					_grid.localCells[i]->Faces.push_back(face.GlobalIndex);
					faceIndex++;
				} else {					
					int fIndex = result.first->second;
					if (cell->GlobalIndex == 98) {
						_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "cell 98 found, face index = ",fIndex);
					};
					_grid.localFaces[fIndex].isExternal = false;
					if (cell->IsDummy) {
						_grid.localFaces[fIndex].FaceCell_2 = cell->GlobalIndex;
					} else {
						//_grid.localFaces[faceIndex].FaceCell_2 = i;
						_grid.localFaces[fIndex].FaceCell_2 = cell->GlobalIndex;
					};
					_grid.localCells[i]->Faces.push_back(fIndex);
				};
			};
		};

		//Finish generating external interprocessor faces
		std::stringstream msg;
		msg.clear();
		int nExternal = 0;
		for (Face& face : _grid.localFaces) if (face.isExternal) {
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External faceID ", face.GlobalIndex);				
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External face FaceCell_1 ", face.FaceCell_1);
			nExternal++;
			int faceNodes = face.FaceNodes.size();
			msg<<"faceNodes = "<<faceNodes<<"\n";
			Cell& cell = _grid.Cells[face.FaceCell_1];
			for (int nCellID : cell.NeigbourCells) {
				msg<<"nCellID = "<<faceNodes<<"\n";
				Cell& nCell = _grid.Cells[nCellID];
				int nCommonNodes = 0;
				std::set<int> cellNodes;
				for (int cellNodeID : nCell.Nodes) {
					cellNodes.insert(cellNodeID);
					msg<<"cellNodeID = "<<cellNodeID<<"\n";
				};
				for (int nodeID : face.FaceNodes) {
					msg<<"nodeID = "<<nodeID<<"\n";
					if (cellNodes.find(nodeID) != cellNodes.end()) nCommonNodes++;
				};
				if (nCommonNodes == faceNodes) {
					//We found second cell
					face.FaceCell_2 = nCellID;
					break;
				};
			};
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External face FaceCell_2 ", face.FaceCell_2);
		};

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Debug
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External faces ", nExternal);	

		//Update cells geometric properties
		for (int i = 0; i < _grid.nCellsLocal; i++) {
			_grid.ComputeGeometricProperties(_grid.localCells[i]);
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.localCells[i]->getInfo());
		};		

		//Generate face geometric properties		
		_grid.nFaces = faceIndex;		
		for (Face& face : _grid.localFaces) {
			int index = face.GlobalIndex;						
			_grid.ComputeGeometricProperties(&_grid.localFaces[index]);

			//Orient normal
			Vector cellCenter = _grid.Cells[_grid.localFaces[index].FaceCell_1].CellCenter;
			Vector faceCenter = _grid.localFaces[index].FaceCenter;
			if (((cellCenter - faceCenter) * _grid.localFaces[index].FaceNormal) > 0) {
				_grid.localFaces[index].FaceNormal *= -1;
			};

			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, face.getInfo());
		};


		//Compute dummy cell geometric properties
		for (int i = _grid.nCellsLocal; i < _grid.localCells.size(); i++) {			
			int neighbour = _grid.localCells[i]->NeigbourCells[0];
			Cell& cell = _grid.Cells[neighbour];
			_grid.localCells[i]->CellVolume = cell.CellVolume;

			//Reflect cell center over boundary face plane
			Face& face = _grid.localFaces[ _grid.localCells[i]->Faces[0]];
			Vector dR = ((cell.CellCenter - face.FaceCenter) * face.FaceNormal) * face.FaceNormal / face.FaceNormal.mod();
			Vector dummyCenter = cell.CellCenter - 2 * dR;
			_grid.localCells[i]->CellCenter = dummyCenter;
			//s_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.localCells[i]->getInfo());
		};

		//Synchronize
		_parallelHelper.Barrier();
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Geometry generation finished");
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating local geometry and faces");	
		return TURBO_OK;
	};

	turbo_errt InitParallelExchange() {
		//Determine values required (TO DO Generalize) closest neighbours for now
		std::set<int> requiredValuesCells;
		std::set<int> requiredGradientsCells;
		for (Cell* cell : _grid.localCells) {
			for (int c : cell->NeigbourCells) {
				if (_parallelHelper.getRank() != _grid.cellsPartitioning[c]) {
					requiredValuesCells.insert(c);
					requiredGradientsCells.insert(c);				
				};
			};
		};

		//Set partitioning info
		_parallelHelper.SetCellsPartitioning(_grid.cellsPartitioning);
		
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
		_parallelHelper.InitExchange();

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Parallel data exchange initialized");	
		return TURBO_OK;
	};

	turbo_errt ParallelExchangeValues() {		
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
				_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "cellGlobalIndex = ", cellGlobalIndex);	
				int cellLocalIndex = _grid.cellsGlobalToLocal[cellGlobalIndex];
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
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Finished : Exchange values.");	
		return TURBO_OK;
	};

	turbo_errt ParallelExchangeGradients() {
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
				int cellLocalIndex = _grid.cellsGlobalToLocal[cellGlobalIndex];
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

	turbo_errt ReadConfiguration(std::string fname) {
		//Hardcode configuration for now
		_configuration.InputCGNSFile = "";
		_configuration.OutputCGNSFile = "result.cgns";

		//Ideal gas law
		_configuration.GasModel = CaloricallyPerfect;
		_configuration.IdealGasConstant = 8.3144621;
		_configuration.SpecificHeatRatio = 1.4;
		_configuration.SpecificHeatVolume = 1006.43 / 1.4;
		_configuration.SpecificHeatPressure = 1006.43;

		//Solver settings
		_simulationType = TimeAccurate;
		CFL = 0.5;
		RungeKuttaOrder = 1;
		MaxIteration = 10000;
		MaxTime = 0.2;
		SaveSolutionSnapshotIterations = 1;

		//Initialize start moment
		stepInfo.Time = 0.0;

		//Boundary conditions
		//_configuration.BoundaryConditions["top_left"].BoundaryConditionType = BCOutflow;
		//_configuration.BoundaryConditions["top_left"].SetPropertyValue("StaticPressure", 1905);
		_configuration.BoundaryConditions["top"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["bottom"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["left"].BoundaryConditionType = BCType_t::BCSymmetryPlane;
		_configuration.BoundaryConditions["right"].BoundaryConditionType = BCType_t::BCSymmetryPlane;

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating local geometry and faces");	
		return TURBO_OK;
	};

	turbo_errt InitBoundaryConditions() {
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Initializing boundary conditions");	

		//Create Boundary conditions from configuration
		for (const auto& bc : _grid.patchesNames) {
			const std::string& bcName = bc.first;			
			int bcMarker = bc.second;
			if ( _configuration.BoundaryConditions.find(bcName) == _configuration.BoundaryConditions.end()) {
				//Boundary condition unspecified
				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unspecified boundary condition");	

				//Synchronize
				_parallelHelper.Barrier();		
				return turbo_errt::TURBO_ERROR;
			}

			//Initialize boundary conditions
			BoundaryConditionConfiguration& bcConfig = _configuration.BoundaryConditions[bcName];
			BCType_t bcType = bcConfig.BoundaryConditionType;
			if (bcType == BCType_t::BCSymmetryPlane) {				
				BoundaryConditions::BCSymmetryPlane* bc = new BoundaryConditions::BCSymmetryPlane();
				_boundaryConditions[bcMarker] = bc; 
			};

			//Attach needed data structures to boundary condition class
			_boundaryConditions[bcMarker]->setGrid(_grid);
			_boundaryConditions[bcMarker]->setGasModel(_gasModel);
		};

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished initialization of boundary conditions");	
		return turbo_errt::TURBO_OK;
	};

	turbo_errt ReadInitialConditions(std::string solutionName) {		
		//Read physical data
		std::vector<double> density(_grid.nCellsLocal);
		std::vector<double> velocityX(_grid.nCellsLocal);
		std::vector<double> velocityY(_grid.nCellsLocal);
		std::vector<double> velocityZ(_grid.nCellsLocal);		
		_cgnsReader.ReadField(_grid, solutionName, "Density", density);
		_cgnsReader.ReadField(_grid, solutionName, "VelocityX", velocityX);
		_cgnsReader.ReadField(_grid, solutionName, "VelocityY", velocityY);
		_cgnsReader.ReadField(_grid, solutionName, "VelocityZ", velocityZ);
		_cgnsReader.ReadField(_grid, solutionName, "VelocityZ", velocityZ);

		//Fill data structures
		nVariables = 5;
		Values.resize(nVariables * _grid.nCellsLocal);
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = density[i];
			double rou = velocityX[i] * ro;
			double rov = velocityY[i] * ro;
			double row = velocityZ[i] * ro;
			Values[i + 0] = ro;
			Values[i + 1] = rou;
			Values[i + 2] = rov;
			Values[i + 3] = row;
		};

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished reading initial conditions");	
		return TURBO_OK;
	};

	turbo_errt GenerateInitialConditions(InitialConditions::InitialConditions& initConditions) {
		//Attach data structures
		initConditions.setGrid(_grid);
		initConditions.setGasModel(_gasModel);

		//For each local cell write initial conditions
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			Cell* c = _grid.localCells[i];
			for (int j = 0; j<_gasModel.nConservativeVariables; j++) {
				std::vector<double> values = initConditions.getInitialValues(*c);
				Values[i * _gasModel.nConservativeVariables + j] = values[j];
			};
		};

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating initial conditions");	
		return TURBO_OK;
	};

	turbo_errt SaveSolution(std::string fname, std::string solutionName) {	
		//Open file
		_cgnsWriter.OpenFile(fname);
		int nv = _gasModel.nConservativeVariables;
		std::vector<std::string> storedFields = _gasModel.GetStoredFieldsNames();

		//Write simulation type

		//Write iterative data
		double iter = stepInfo.Iteration;
		double time = stepInfo.Time;		
		if (_simulationType == TimeAccurate) {
			_cgnsWriter.WriteIterativeData(time, iter);
		};

		//Write solution node		
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "1");	
		_cgnsWriter.WriteSolution(_grid, solutionName);
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "2");	

		//Write physical quantities
		std::vector<double> buffer(_grid.nCellsLocal);
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "3");	

		//Stored fields
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "storedFields.size() = ", storedFields.size());	
		for (int fieldIndex = 0; fieldIndex < storedFields.size(); fieldIndex++ ) {
			std::string fieldName = storedFields[fieldIndex];			
			//Write to buffer
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Field name " + fieldName);	
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "nCellsLocal = ", _grid.nCellsLocal);				
			for (int i = 0; i<_grid.nCellsLocal; i++) {
				GasModel::ConservativeVariables U(&Values[i * nv]);
				std::vector<double> storedValues = _gasModel.GetStoredValuesFromConservative(U);				
				buffer[i] = storedValues[fieldIndex];
			};
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "2");	
			_cgnsWriter.WriteField(_grid, solutionName, fieldName, buffer);		
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "3");	
		};

		//Additional fields (TO DO) generalize
		//VelocityX
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double rou = Values[i * nv + 1];
			double u = rou / ro;
			buffer[i] = u;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "VelocityX", buffer); 

		//VelocityY
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double rov = Values[i * nv + 2];
			double v = rov / ro;
			buffer[i] = v;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "VelocityY", buffer); 

		//VelocityZ
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double row = Values[i * nv + 3];
			double w = row / ro;
			buffer[i] = w;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "VelocityZ", buffer); 

		//Pressure		
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i * nv + 0];
			double rou = Values[i * nv + 1];
			double rov = Values[i * nv + 2];
			double row = Values[i * nv + 3];
			double roE = Values[i * nv + 4];
			double u = rou / ro;
			double v = rov / ro;
			double w = row / ro;
			double P = (roE - ro*(w*w + v*v + u*u)/2.0) / (_gasModel.Gamma - 1.0);
			buffer[i] = P;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "Pressure", buffer); 

		//Write partitioning to solution
		std::vector<double> part(_grid.nCellsLocal, _parallelHelper.getRank());		
		_cgnsWriter.WriteField(_grid, solutionName, "Processor", part); 

		//Close file
		_cgnsWriter.CloseFile();

		return TURBO_OK;
	};	

	//Compute convective flux and max wave propagation speed throught each face
	void ComputeConvectiveFluxes(std::vector<std::vector<double>>& fluxes, std::vector<double>& maxWaveSpeed, std::vector<double>& cellValues) {		
		//Compute gradients for second order reconstruction
		/*if (IsSecondOrder) {
			ComputeFunctionGradient(gradCellsRo, U, &Model<RiemannSolver>::GetDensity);
			ComputeFunctionGradient(gradCellsRoU, U, &Model<RiemannSolver>::GetRoU);
			ComputeFunctionGradient(gradCellsRoV, U, &Model<RiemannSolver>::GetRoV);
			ComputeFunctionGradient(gradCellsRoW, U, &Model<RiemannSolver>::GetRoW);
			ComputeFunctionGradient(gradCellsRoE, U, &Model<RiemannSolver>::GetRoE);
		};*/

		//Apply boundary conditions (TO DO) decide where to store dummy values
		for (int i = 0; i < _grid.dummyLocalCells.size(); i++) {
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
		for (int i = 0; i<_grid.localFaces.size(); i++) {
			Face& f = _grid.localFaces[i];			
			
			//Debug
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "FaceID = ", f.GlobalIndex);	
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "FaceCell1 = ", cellIndexLeft);	
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "FaceCell2 = ", cellIndexRight);	
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "isExternal = ", f.isExternal);
			
			std::vector<double> flux;
			GasModel::ConservativeVariables UL;
			GasModel::ConservativeVariables UR;
			//TO DO check for indexes LOCAL vs GLOBAL for cells
			//UL = GasModel::ConservativeVariables(&cellValues[cellIndexLeft * nVariables]);
			/*msg.str("Reconstruct cell values.");
			msg<<"cLeft.GlobalIndex = "<< cLeftGlobalIndex << "\n";
			msg<<"cRight.GlobalIndex = "<< cRightGlobalIndex << "\n";
			msg<<"cellIndexLeft = "<< cellIndexLeft << "\n";
			msg<<"cellIndexRight = "<< cellIndexRight << "\n";
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str() );*/
			UL = GasModel::ConservativeVariables(GetCellValues(f.FaceCell_1));
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "UL reconstructed");
			UR = GasModel::ConservativeVariables(GetCellValues(f.FaceCell_2));
			//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "UR reconstructed");

			//GasModel::ConservativeVariables UL_ = GasModel::ConservativeVariables(&cellValues[cellIndexLeft * nVariables]);
			//GasModel::ConservativeVariables UR_;
			//if (!f.isExternal) {			 
			//	//If it is local face and both cells local				
			//	if (cRight.IsDummy) {
			//		//Dummy cell
			//		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Dummy cell");	
			//		std::vector<BoundaryConditions::BoundaryConditionResultType> resultType = _boundaryConditions[cRight.BCMarker]->bcResultTypes;				
			//		UR_ = GasModel::ConservativeVariables(_boundaryConditions[cRight.BCMarker]->getDummyValues(UL, f));
			//	} else {
			//		//Ordinary cell
			//		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Local cell");	
			//		UR_ = GasModel::ConservativeVariables(&cellValues[cellIndexRight * nVariables]);
			//	};
			//} else {
			//	//If it is interprocessor face obtain values from parallel helper
			//	//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Requested cell");	
			//	int globalIndex = cRight.GlobalIndex;
			//	UR_ = _parallelHelper.RequestedValues[globalIndex];
			//};

			//Compute flux
			flux = rSolver.ComputeFlux( UL,  UR, f);	
			/*msg.str("");		
			msg<<"Flux for face = "<<f.GlobalIndex<<" Cell1 = "<<f.FaceCell_1<<" Cell2 = "<<f.FaceCell_2<<" computed. Flux = ("
				<<flux[0]<<" , "<<flux[1]<<" , "<<flux[2]<<" , "<<flux[3]<<" , "<<flux[4]<<") FaceNormal = ("
				<<f.FaceNormal.x<<","<<f.FaceNormal.y<<","<<f.FaceNormal.z<<")";
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());*/

			//Store wave speeds
			maxWaveSpeed[f.GlobalIndex] = rSolver.MaxEigenvalue;

			//Store flux
			fluxes[f.GlobalIndex] = flux;
			
			//Vector dRL = f.FaceCenter - _grid.cells[f.FaceCell_1].CellCenter;
			//if (f.isExternal) {				
			//	//Apply boundary conditions
			//	UL = cellValues[f.FaceCell_1]; 						
			//	if (!IsSecondOrder) {
			//		//Constant reconstruction
			//	} else {
			//		//Try linear reconstruction
			//		dRL = Vector(0,0,0);
			//		UL.ro = UL.ro + gradCellsRo[f.FaceCell_1] * dRL;
			//		UL.rou = UL.rou + gradCellsRoU[f.FaceCell_1] * dRL;
			//		UL.rov = UL.rov + gradCellsRoV[f.FaceCell_1] * dRL;
			//		UL.row = UL.row + gradCellsRoW[f.FaceCell_1] * dRL;
			//		UL.roE = UL.roE + gradCellsRoE[f.FaceCell_1] * dRL;
			//	};
			//	flux = _boundaryConditions[f.BCMarker]->ComputeConvectiveFlux(f, UL);				
			//} else {
			//	UL = cellValues[f.FaceCell_1];
			//	UR = cellValues[f.FaceCell_2];	
			//	if (!IsSecondOrder) {
			//		//Constant reconstruction
			//	} else {
			//		//Try linear reconstruction					
			//		UL.ro = UL.ro + gradCellsRo[f.FaceCell_1] * dRL;
			//		UL.rou = UL.rou + gradCellsRoU[f.FaceCell_1] * dRL;
			//		UL.rov = UL.rov + gradCellsRoV[f.FaceCell_1] * dRL;
			//		UL.row = UL.row + gradCellsRoW[f.FaceCell_1] * dRL;
			//		UL.roE = UL.roE + gradCellsRoE[f.FaceCell_1] * dRL;
			//	
			//		Vector dRR = f.FaceCenter - _grid.cells[f.FaceCell_2].CellCenter;					
			//		UR.ro = UR.ro + gradCellsRo[f.FaceCell_2] * dRR;
			//		UR.rou = UR.rou + gradCellsRoU[f.FaceCell_2] * dRR;
			//		UR.rov = UR.rov + gradCellsRoV[f.FaceCell_2] * dRR;
			//		UR.row = UR.row + gradCellsRoW[f.FaceCell_2] * dRR;
			//		UR.roE = UR.roE + gradCellsRoE[f.FaceCell_2] * dRR;				
			//	};
			//	flux = rSolver.ComputeFlux( UL,  UR, f);
			//};				
					
		};

		//Synchronyze
		_parallelHelper.Barrier();
	};	

	////Compute residual for each cell		
	void ComputeResidual(std::vector<double>& residual, std::vector<double>& cellValues) {		
		//Compute convective fluxes and max wave speeds
		for (std::vector<double>& flux : FaceFluxes) {
			for (double& v : flux) v = 0;
		};
		ComputeConvectiveFluxes(FaceFluxes, MaxWaveSpeed, cellValues);
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Convective fluxes calculated");

		////Compute gradients
		//if (IsGradientRequired) ComputeGradients();

		////Compute viscous fluxes
		//ComputeViscousFluxes();	

		//Compute residual for each cell		
		for ( int cellIndex = 0; cellIndex<_grid.nCellsLocal; cellIndex++ )
		{
			Cell* cell = _grid.localCells[cellIndex];			
			for (int i = 0; i < nVariables; i++) residual[cellIndex*nVariables + i] = 0;
			std::vector<int>& nFaces = cell->Faces;
			for (int nFaceIndex : nFaces)
			{
				Face& face = _grid.localFaces[nFaceIndex];
				int fluxDirection = (face.FaceCell_1 == cell->GlobalIndex) ? 1 : -1;		
				std::vector<double> fluxc = FaceFluxes[nFaceIndex];
				//std::vector<double> fluxv = vfluxes[nFaceIndex];
				for (int j = 0; j<nVariables; j++) {
					residual[cellIndex * nVariables + j] +=  (fluxc[j]) * face.FaceSquare * fluxDirection;
				};
				
			};
		};

		return;
	};

	//Compute spectral radii estimate for each cell
	void ComputeSpectralRadius(std::map<int, double>& spectralRadius, std::vector<double>& maxWaveSpeed, std::vector<double>& cellValues) {
		spectralRadius.clear();		
		for (int cellIndex = 0; cellIndex<_grid.nCellsLocal; cellIndex++)
		{
			Cell* cell = _grid.localCells[cellIndex];			
			spectralRadius[cellIndex] = 0;
			std::vector<int>& nFaces = cell->Faces;
			for (int nFaceIndex : nFaces)
			{
				//Blazek f. 6.21
				Face& face = _grid.localFaces[nFaceIndex];			
				spectralRadius[cellIndex] +=  maxWaveSpeed[nFaceIndex] * face.FaceSquare;
			};
		};
	};

	//Compute local time step for each cell utilizing spectral radius estimates
	void ComputeLocalTimeStep(std::map<int, double>& localTimeStep, std::map<int, double>& spectralRadius) {
		localTimeStep.clear();
		for (int cellIndex = 0; cellIndex<_grid.nCellsLocal; cellIndex++)
		{
			Cell* cell = _grid.localCells[cellIndex];			
			double sR = spectralRadius[cellIndex];
			localTimeStep[cellIndex] = CFL * cell->CellVolume / sR; //Blazek f. 6.20
		}
	};	

	//Explicit time step
	void ExplicitTimeStep() {
		std::stringstream msg;

		//Compute residual		
		ComputeResidual(Residual, Values);

		//Synchronize
		_parallelHelper.Barrier();

		//Determine time step as global minimum over local time steps		
		std::map<int, double> spectralRadius;
		ComputeSpectralRadius(spectralRadius, MaxWaveSpeed, Values);

		/*msg.clear();
		msg.str(std::string());
		msg<<"MaxWaveSpeed :\n";
		int ind = 0;
		for (auto p : MaxWaveSpeed) {
			msg<<"Face = "<<ind++<<" MaxWaveSpeed = "<<p<<"\n";
		};
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	*/

		//Synchronize
		_parallelHelper.Barrier();

		std::map<int, double> localTimeStep;
		ComputeLocalTimeStep(localTimeStep, spectralRadius);
		
		/*msg.clear();
		msg.str(std::string());
		msg<<"Local time step :\n";
		for (auto p : localTimeStep) {
			msg<<"cell= "<<p.first<<" timeStep= "<<p.second<<"\n";
		};*/
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Synchronize
		_parallelHelper.Barrier();

		stepInfo.TimeStep = MaxTime - stepInfo.Time; 
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
		
		//Time stepping
		for (int stage = 0; stage<nStages-1; stage++) {			
			for ( int cellIndex = 0; cellIndex<_grid.nCellsLocal; cellIndex++ )
			{
				Cell* cell = _grid.localCells[cellIndex];			
				//Update values
				for (int i = 0; i < nVariables; i++) {
					Values[cellIndex*nVariables + i] += Residual[cellIndex*nVariables + i] * (-stepInfo.TimeStep / cell->CellVolume) * alpha[stage];
				};
			};

			//Compute new residual
			ComputeResidual(Residual, Values);
			//Synchronize
			_parallelHelper.Barrier();
		};

		//Compute new result and residual
		stepInfo.Residual = std::vector<double>(5, 0.0);		
		for ( int cellIndex = 0; cellIndex<_grid.nCellsLocal; cellIndex++ )
		{
			Cell* cell = _grid.localCells[cellIndex];			
			//Update values
			for (int i = 0; i < nVariables; i++) {
				Values[cellIndex*nVariables + i] += Residual[cellIndex*nVariables + i] * (-stepInfo.TimeStep / cell->CellVolume) * alpha[nStages-1];
				stepInfo.Residual[i] += pow(Residual[cellIndex*nVariables + i], 2);
			};			
		};

		//Compute RMS residual		
		for (int i = 0; i<nVariables; i++) {
			stepInfo.Residual[i] = _parallelHelper.SumDouble(stepInfo.Residual[i]);
			stepInfo.Residual[i] = sqrt(stepInfo.Residual[i]);
		};

		//Advance total time
		stepInfo.Time += stepInfo.TimeStep;

		//Synchronize
		_parallelHelper.Barrier();
	};	

	//Implicit time step
	void ImplicitTimeStep() {
	};
	
	//Compute scalar function gradient in each cell	
	//void ComputeFunctionGradient( std::vector<Vector>& grads, std::vector<double>& values, double (*func)(const std::vector<double>&) ) {
	//	//Allocate memory for gradients
	//	grads.resize(_grid.nCellsLocal);		

	//	//For each cell compute gradient of given function
	//	std::vector<Vector> nPoints;
	//	std::vector<double> nValues;
	//	for (int i = 0; i<_grid.nCellsLocal; i++) {			
	//		Cell& cell = *_grid.localCells[i];			

	//		//Determine required set of point and values
	//		nPoints.clear();			
	//		nValues.clear();			

	//		//Add all neighbours
	//		for (int j = 0; j<cell.NeigbourCells.size(); j++) {				
	//			Cell& nCell = _grid.Cells[cell.NeigbourCells[j]];
	//			std::vector<double> nU;
	//			if (nCell.
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
	//					nPoint = _grid.cells[face.FaceCell_2].CellCenter;
	//				} else {
	//					nPoint = _grid.cells[face.FaceCell_1].CellCenter;
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

	inline bool IsLocalCell(int globalIndex) {
		return _grid.cellsPartitioning[globalIndex] == _parallelHelper.getRank();
	};

	inline bool IsDummyCell(int globalIndex) {
		return _grid.Cells[globalIndex].IsDummy;
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
				Cell& cell = _grid.Cells[globalIndex];				
				//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "!" );
				return _boundaryConditions[cell.BCMarker]->getDummyValues(Values, cell);
			} else {
				//If proper cell return part of Values array				
				if (localIndex == -1) {
					//If local index is unknown determine it
					//Lower perfomance WARNING
					localIndex = _grid.cellsGlobalToLocal[globalIndex];					
					//msg.str("");					
					//msg<<"Local Index = "<<localIndex;
					//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str() );
				};
				return std::vector<double>(&Values[localIndex * _gasModel.nConservativeVariables], &Values[localIndex * _gasModel.nConservativeVariables] + _gasModel.nConservativeVariables);
			};
		} else {
			//Return result of interprocessor exchange
			return _parallelHelper.RequestedValues[globalIndex];
		};
	};
};

#endif
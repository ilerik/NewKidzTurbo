#ifndef TURBO_MAINSOLVER_KERNEL
#define TURBO_MAINSOLVER_KERNEL
#include "grid.h"
#include "mpi.h"
#include "logger.h"
#include "parmetis.h"
#include "CGNSReader.h"
#include "CGNSWriter.h"
#include "parallelHelper.h"

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

	//Task parameters
	int _rungeKuttaOrder;	

	//Parallel run information		
	ParallelHelper _parallelHelper;
	int _rank;
	int _nProcessors;

	//Internal storage		
	int nVariables; //number of variables in each cell
	std::vector<double> Values;	//Conservative variables
	std::vector<double> Residual; //Residual values
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
		
		msg.clear();
		msg<<part.size()<<"\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg<<"recvcounts[] = \n";
		for (int i = 0; i<_nProcessors; i++) msg<<recvcounts[i]<<" ";
		msg<<"\n";

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	
		_parallelHelper.Allgatherv( part, recvcounts, _grid.cellsPartitioning);

		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good.");	

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
		msg.clear();
		msg<<"Partitioning = [ ";
		for (int p : _grid.cellsPartitioning) {
			msg<<p<<" ";
		};
		msg<<"]\n";
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());		


		//Extract local cells indexes		
		_grid.localCells.clear();
		_grid.localCellIndexes.clear();
		for (int i = 0; i<_grid.nCells; i++) {
			if ((_grid.cellsPartitioning[i] == _rank) && (!_grid.Cells[i].IsDummy)) {
				_grid.localCells.push_back(&_grid.Cells[i]);
				_grid.localCellIndexes.push_back(i);
				_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.Cells[i].getInfo());
			};
		};		
		_grid.nCellsLocal = _grid.localCellIndexes.size(); //Without dummy cells
		for (int i = 0; i<_grid.nCells; i++) {
			if ((_grid.cellsPartitioning[i] == _rank) && (_grid.Cells[i].IsDummy)) {
				_grid.localCells.push_back(&_grid.Cells[i]);
				_grid.localCellIndexes.push_back(i);
				_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.Cells[i].getInfo());
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
			_grid.localCells[i] = &_grid.Cells[_grid.localCellIndexes[i]];			
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
					_grid.localFaces[faceIndex].FaceCell_1 = i;
					_grid.localCells[i]->Faces.push_back(face.GlobalIndex);
					faceIndex++;
				} else {
					int fIndex = result.first->second;
					_grid.localFaces[fIndex].isExternal = false;
					_grid.localFaces[fIndex].FaceCell_2 = i;
					_grid.localCells[i]->Faces.push_back(fIndex);
				};
			};
		};

		//Update cells geometric properties
		for (int i = 0; i < _grid.nCellsLocal; i++) {
			_grid.ComputeGeometricProperties(_grid.localCells[i]);
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.localCells[i]->getInfo());
		};		

		//Generate face geometric properties		
		_grid.nFaces = faceIndex;		
		for (Face& face : _grid.localFaces) {
			int index = face.GlobalIndex;						
			_grid.ComputeGeometricProperties(&_grid.localFaces[index]);

			//Orient normal
			Vector cellCenter = _grid.localCells[_grid.localFaces[index].FaceCell_1]->CellCenter;
			Vector faceCenter = _grid.localFaces[index].FaceCenter;
			if (((cellCenter - faceCenter) * _grid.localFaces[index].FaceNormal) > 0) {
				_grid.localFaces[index].FaceNormal *= -1;
			};

			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, face.getInfo());
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
			_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, _grid.localCells[i]->getInfo());
		};

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Geometry generation finished");
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating local geometry and faces");	
		return TURBO_OK;
	};

	turbo_errt InitParallelExchange() {
		//Determine values required

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Parallel data exchange initialized");	
		return TURBO_OK;
	};



	turbo_errt ReadBoundaryConditions() {

		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating local geometry and faces");	
		return TURBO_OK;
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

	turbo_errt SaveSolution(std::string fname, std::string solutionName) {	
		//Open file
		_cgnsWriter.OpenFile(fname);

		//Write physical quantities
		//Density
		std::vector<double> buffer(_grid.nCellsLocal);
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i + 0];
			buffer[i] = ro;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "Density", buffer); 

		//VelocityX
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i + 0];
			double rou = Values[i + 1];
			double u = rou / ro;
			buffer[i] = u;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "VelocityX", buffer); 

		//VelocityY
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i + 0];
			double rov = Values[i + 2];
			double v = rov / ro;
			buffer[i] = v;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "VelocityY", buffer); 

		//VelocityZ
		for (int i = 0; i<_grid.nCellsLocal; i++) {
			double ro = Values[i + 0];
			double row = Values[i + 3];
			double w = row / ro;
			buffer[i] = w;
		};
		_cgnsWriter.WriteField(_grid, solutionName, "VelocityZ", buffer); 

		//Write partitioning to solution
		std::vector<double> part(_grid.nCellsLocal, _parallelHelper.getRank());		
		_cgnsWriter.WriteField(_grid, solutionName, "Processor", part); 

		//Close file
		_cgnsWriter.CloseFile();

		return TURBO_OK;
	};

	//Compute residual
	void ComputeResidual(double *R, double *U){
		
	};

	//Compute convective flux through face
	void ComputeConvectiveFlux() {
	};	

	//Explicit time step
	void ExplicitTimeStep() {
	};

	//Implicit time step
	void ImplicitTimeStep() {
	};

	//Medium properties and features


};

#endif
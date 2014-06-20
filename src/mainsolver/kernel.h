#ifndef TURBO_MAINSOLVER_KERNEL
#define TURBO_MAINSOLVER_KERNEL
#include "grid.h"
#include "mpi.h"
#include "logger.h"
#include "parmetis.h"
#include "CGNSReader.h"
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

public:
	//Initialize computational kernel
	turbo_errt Initilize(int *argc, char **argv[]) {		
		try {
			//Initialize parallel MPI subsystem
			_parallelHelper.Init(argc, argv);
			_rank = _parallelHelper.getRank();
			_nProcessors = _parallelHelper.getProcessorNumber();

			//Initialize loggin subsistem
			_logfilename = "kernel.log"; //TO DO
			_logger.InitLogging(_logfilename, _rank);

			//Initialize cgns i\o subsystem
			_cgnsReader.Init(_logger, _parallelHelper);
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
		_parallelHelper.Finalize();			
		return TURBO_OK;
	};

	//Load grid from CGNS file
	turbo_errt LoadGrid(std::string filename) {
		_grid = _cgnsReader.LoadGrid(filename, _parallelHelper);
		PartitionGrid();
		//_grid.UpdateGeometricProperties();
		return TURBO_OK;
	};

	//Load grid topology info
	turbo_errt LoadGridTopologyAndInfo(std::string filename) {		
	//	//Open cgns file
	//	_grid.gridInfo.fileLocation = filename;
	//	int cgfile, mode = CG_MODE_READ;		
	//	_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Opening cgns file"+filename+"");		
	//	if (cg_open (filename.c_str(), mode, &cgfile)) {
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Cannot open grid file.");		
	//		cg_error_exit();
	//	};
	//			
	//	//Determine the number of bases in the grid. This example assumes 
	//	//one base. However it is allowed to have multiple bases.    
	//	int nBases;
	//	if(cg_nbases(cgfile, &nBases)!= CG_OK) {			
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Cannot read number of bases.");		
	//		cg_error_exit();
	//	};
	//	if (nBases != 1) {
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Cannot process grid. Must be one base.");
	//		return TURBO_ERROR;
	//	};
	//	int base = 1;	
	//	_grid.gridInfo.MainBaseIndex = base;

	//	//Check the cell and physical dimensions of the base.
	//	int physDim;
	//	int cellDim;
	//	char cgnsName[255];
	//	if(cg_base_read(cgfile, base, cgnsName, &cellDim, &physDim) != CG_OK) {
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Cannot read base info.");
	//		cg_error_exit();
	//		return TURBO_ERROR;
	//	};		
	//	_grid.gridInfo.MainBaseName = std::string(cgnsName);
	//	_grid.gridInfo.CellDimensions = cellDim;
	//	_grid.gridInfo.GridDimensions = physDim;


	//	//Read the number of zones in the grid.
	//	//This example assumes one zone.
	//	int nZones;
	//	if(cg_nzones(cgfile, base, &nZones) != CG_OK) {
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Cannot read number of zones.");
	//		cg_error_exit();     
	//	};
	//	if(nZones != 1) {
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "This CGNS loader assumes one zone.");
	//		return TURBO_ERROR;
	//	}
	//	int zone = 1;
	//	_grid.gridInfo.MainZoneIndex = 1;

	//	//Check the zone type. This should be Unstructured.
	//	ZoneType_t zoneType;
	//	if(cg_zone_type(cgfile, base, zone, &zoneType) != CG_OK) {
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Cannot read zone type.");
	//		cg_error_exit(); 		
	//	};
	//	if(zoneType != Unstructured) {			
	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unstructured zone expected.");
	//		return TURBO_ERROR;
	//	};

	//	//Determine the number of vertices and volume elements in this */
	//	//zone (and thus in the grid, because one zone is assumed). */
	//	char zoneName[255];
	//	cgsize_t sizes[3];
	//	if(cg_zone_read(cgfile, base, zone, zoneName, sizes) != CG_OK) cg_error_exit();		
	//	int nVertices = sizes[0];
	//	int nCells = sizes[1];

	//	//Save to grid info
	//	_grid.gridInfo.MainZoneName = std::string(zoneName);
	//	_grid.gridInfo.nNodes = nVertices;
	//	_grid.gridInfo.nCells = nCells;

	//	//Determine the number and names of the coordinates.
	//	int nCoords;
	//	if(cg_ncoords(cgfile, base, zone, &nCoords) != CG_OK)
	//		cg_error_exit();	

	//	char name[255];
	//	DataType_t dataType;
	//	if(cg_coord_info(cgfile, base, zone, 1, &dataType, name) != CG_OK)
	//		cg_error_exit();			

	//	//Load shared grid data 
	//	//Connectivity info, etc.

	//	//Read nodes coordinates
	//	/* Read the x-coordinates. The y and z-coordinates can be read */
	//	/* similarly. Just replace CoordinateX by CoordinateY and */
	//	/* CoordinateZ respectively. This assumes Cartesian coordinates */
	//	/* in double precision. Note that CGNS starts the numbering at */
	//	/* 1 even if C is used. */
	//
	//	int one = 1;	
	//	std::vector<double> coorX(nVertices);
	//	if(cg_coord_read(cgfile, base, zone, "CoordinateX", RealDouble, &one,
	//					&nVertices, &coorX[0]) != CG_OK) cg_error_exit();	
	//
	//	//Check if we have Y coordinate
	//	std::vector<double> coorY(nVertices);
	//	for (int i = 0; i<nVertices; i++) coorY[i] = 0;
	//	if (nCoords > 1) {
	//		if(cg_coord_read(cgfile, base, zone, "CoordinateY", RealDouble, &one,
	//					&nVertices, &coorY[0]) != CG_OK) cg_error_exit();	
	//	};
	//
	//	//Check if its not 3D case and make all Z equal zero
	//	std::vector<double> coorZ(nVertices);
	//	for (int i = 0; i<nVertices; i++) coorZ[i] = 0;
	//	if (nCoords > 2) {
	//		if(cg_coord_read(cgfile, base, zone, "CoordinateZ", RealDouble, &one,
	//					&nVertices, &coorZ[0]) != CG_OK) cg_error_exit();	
	//	};

	//	//Create grid nodes
	//	_grid.localNodes.resize(nVertices);
	//	for (int i = 0; i<nVertices; i++) {
	//		Node newNode;
	//		newNode.GlobalIndex = i+1;
	//		newNode.P.x = coorX[i];
	//		newNode.P.y = coorY[i];
	//		newNode.P.z = coorZ[i];
	//		_grid.localNodes[i] = newNode;
	//	};
	//	
	//	//Get sections information		
	//	int nSections;
	//	if(cg_nsections(cgfile, base, zone, &nSections) != CG_OK)
	//		cg_error_exit();

	//	// Find volume elements section also check all element types
	//	/* Loop over the number of sections and read the element */
	//	/* connectivities. As CGNS starts the numbering at 1 the */
	//	/* for-loop starts at 1 as well. */		
	//	std::vector<int> conn;
	//	int nNodes;
	//	int nConnNodes;		
	//	int currentCellNumber = 0;
	//	_grid.gridInfo.GlobalIndexToNumber.clear();
	//	_grid.gridInfo.NumberToGlobalIndex.clear();
	//	std::map<int, std::vector<int> > cgnsCellIndexToElementNodes;
	//	std::map<int, ElementType_t > cgnsCellIndexToElementType;
	//	std::map<int, std::vector<int> > cgnsFaceIndexToElementNodes;
	//	std::map<int, ElementType_t > cgnsFaceIndexToElementType;

	//	for(int sec=1; sec<=nSections; sec++)
	//	{
	//		/* Determine the element type and set the pointer for the */
	//		/* connectivity accordingly. */
	//		char secName[255];
	//		ElementType_t type;
	//		cgsize_t eBeg;
	//		cgsize_t eEnd;
	//		int nBdry;
	//		int parentFlag;
	//		if(cg_section_read(cgfile, base, zone, sec, secName, &type,
	//						&eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
	//			cg_error_exit();
	//		int nElements = (eEnd-eBeg+1);
	//		
	//		//Read zone elements and their types and connectivity info
	//		if (type == MIXED) {			
	//			int connDataSize;				
	//			int parentData;
	//			cg_ElementDataSize(cgfile, base, zone, sec, &connDataSize);
	//			std::vector<int> connectivityInfo(connDataSize);
	//			cg_elements_read(cgfile, base, zone, sec, &connectivityInfo[0], &parentData);
	//			nConnNodes = connDataSize - nElements;				
	//			conn.resize(nConnNodes);

	//			//Parse each element connectivity info
	//			int connP = 0;
	//			int connCounter = 0;
	//			for (int ind = eBeg; ind<=eEnd; ind++) {										
	//				ElementType_t elementType = static_cast<ElementType_t>(connectivityInfo[connP++]);
	//				if (!IsSupprotedCGNSEelementType(elementType)) {
	//					_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unsupported element encountered.");
	//					return TURBO_ERROR;	
	//				};
	//				int elememtNodeCount = 0;
	//				cg_npe(elementType, &elememtNodeCount);
	//				for (int i = 0; i<elememtNodeCount; i++) {
	//					conn[connCounter + i] = connectivityInfo[connP++];
	//				};
	//				connCounter++;
	//			};
	//			//if mixed
	//		} else {
	//			switch (type)
	//			{
	//			case QUAD_4: 
	//				nNodes = NPE_QUAD_4;			
	//				break;	
	//			case BAR_2:
	//				nNodes = NPE_BAR_2;			
	//				break;	
	//			case HEXA_8:
	//				nNodes = NPE_HEXA_8;
	//				break;
	//			default:				
	//				_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unsupported element encountered.");
	//				return TURBO_ERROR;				
	//			}
	//		};
	//		

	//		//Create cells if its main section (contains elements with dim = physDim)
	//		//For now we assume that all volume elements are in one section so following condition must hold
	//		//TO DO Additional checks
	//		bool isVolumeElementsSection;
	//		if ((GetCGNSElementDimensions(type) == _grid.gridInfo.CellDimensions)) {
	//			isVolumeElementsSection = true;
	//			_grid.gridInfo.CellsSections.push_back(sec);
	//		} else {
	//			isVolumeElementsSection = false;
	//			_grid.gridInfo.BoundarySections.push_back(sec);
	//		};			
	//		
	//		//Read section connectivity info			
	//		nConnNodes = nNodes*nElements;
	//		conn.resize(nConnNodes); 			
 //  
	//		/* Read the connectivity. Again, the node numbering of the */
	//		/* connectivities start at 1. If internally a starting index */
	//		/* of 0 is used (typical for C-codes) 1 must be substracted */
	//		/* from the connectivities read. */   		
	//		if(cg_elements_read(cgfile, base, zone, sec, &conn[0], NULL) != CG_OK)
	//			cg_error_exit();	

	//		//Add all elements read to approriate structures for cells
	//		/*if (!isVolumeElementsSection) continue;			
	//		for (int ind = eBeg; ind<=eEnd; ind++) {
	//			_grid.gridInfo.GlobalIndexToNumber[ind] = currentCellNumber;
	//			_grid.gridInfo.NumberToGlobalIndex[currentCellNumber] = ind;
	//			numberToElementType[currentCellNumber] = type;
	//			numberToElementNodes[currentCellNumber].resize(nNodes);
	//			for (int j = 0; j<nNodes; j++) numberToElementNodes[currentCellNumber][j] = (conn[currentCellNumber*nNodes + j]);
	//			currentCellNumber++;
	//		};*/
	//	};

	//	//Distribute cells evenly between processors		
	//	int vProc = nCells / _nProcessors;
	//	int vLeft = nCells % _nProcessors;
	//	_grid.vdist.resize(_nProcessors + 1);
	//	_grid.vdist[0] = 0;
	//	for (int i = 0; i<_nProcessors; i++) {
	//		_grid.vdist[i+1] = _grid.vdist[i] + vProc;
	//		if (vLeft > 0) {
	//			vdist[i+1]++;
	//			vLeft--;
	//		};
	//	};
	//	vdist[_nProcessors] = nCells;
	//	
	//	//TO DO Generalize for different types of volume elements
	//	//Convert connectivity info to graph
	//	_nLocalCells = vdist[_rank + 1] - vdist[_rank];
	//	idx_t *eptr;
	//	idx_t *eind;
	//	idx_t ne = _nLocalCells;
	//	idx_t nn = 0;
	//	idx_t ncommon = _grid.gridInfo.CellDimensions;
	//	idx_t numflag = 0;		

	//	//Fill in connectivity data structures		
	//	/*eptr = new idx_t[ne+1];			
	//	eptr[0] = 0;		
	//	int localI = 0;
	//	for (int i = vdist[_rank]; i < vdist[_rank + 1]; i++)
	//	{		
	//		int nCellNodes = numberToElementNodes[i].size();
	//		nConnNodes += nCellNodes;
	//		eptr[localI+1] = eptr[localI] + nCellNodes;			
	//		localI++;
	//	};
	//	nn = eptr[ne];
	//	eind = new idx_t[nn];
	//	localI = 0;
	//	for (int i = vdist[_rank]; i < vdist[_rank + 1]; i++)
	//	{		
	//		int nCellNodes = numberToElementNodes[i].size();
	//		for (int j = 0; j<nCellNodes; j++) {			
	//			eind[eptr[localI] + j] = numberToElementNodes[i][j];
	//		};
	//		localI++;
	//	};*/

	//	//Call graph generation function
	//	//if (_nProcessors == 1) {
	//	//	idx_t* _xadj;
	//	//	idx_t* _adjncy;
	//	//	int result = METIS_MeshToDual(&ne, &nn, eptr, eind, &ncommon, &numflag, &_xadj, &_adjncy);			
	//	//	if (result != METIS_OK) {
	//	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "METIS_MeshToDual failed.");
	//	//		return TURBO_ERROR;
	//	//	};
	//	//} else {
	//	//	int result = ParMETIS_V3_Mesh2Dual(&vdist[0], eptr, eind, &numflag, &ncommon, &&xadj[0], &&adjncy[0], &_comm);
	//	//	if (result != METIS_OK) {
	//	//		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "ParMETIS_V3_Mesh2Dual failed.");
	//	//		return TURBO_ERROR;
	//	//	};
	//	//};

	//	//Free memory
	//	free(eptr);
	//	free(eind);

	//	/* Determine the number of boundary conditions for this zone. */

	//	int nBocos;
	//	if(cg_nbocos(cgfile, base, zone, &nBocos) != CG_OK)
	//		cg_error_exit();

	//	/* Loop over the number of boundary conditions. */

	//	for(int boco=1; boco<=nBocos; boco++)
	//	{
	//		/* Read the info for this boundary condition. */
	//		char bocoName[255];
	//		BCType_t bocoType;
	//		PointSetType_t ptsetType;
	//		int normalIndex;
	//		cgsize_t nBCElem;
	//		cgsize_t normListFlag;
	//		DataType_t normDataType;
	//		int nDataSet;
	//		if(cg_boco_info(cgfile, base, zone, boco, bocoName, &bocoType, &ptsetType,
	//						&nBCElem, &normalIndex, &normListFlag, &normDataType,
	//						&nDataSet) != CG_OK) cg_error_exit();

	//		/* Read the element ID's. */

	//		std::vector<cgsize_t> BCElemRead(nBCElem);
	//		if(cg_boco_read(cgfile, base, zone, boco, &BCElemRead[0], NULL) != CG_OK) cg_error_exit();

	//		/* And much more to make it fit into the */
	//		/* internal data structures. */

	//		//Assign to all faces boundary marker
	//		bool notBoundary = false;
	//		//for (int i = 0; i<nBCElem; i++) {
	//		//	std::set<int> idx = indexToElementNodes[BCElemRead[i]];
	//		//	if (faces.find(idx) != faces.end()) {
	//		//		//Element is a face
	//		//		int globalInd = faces[idx];
	//		//		_grid.faces[globalInd].BCMarker = boco;
	//		//		if (!_grid.faces[globalInd].isExternal) throw Exception("Boundary face must be marked isExternal");
	//		//	} else {
	//		//		//Its not a face so this is not valid Boundary condition
	//		//		//??
	//		//		notBoundary = true;
	//		//	};
	//		//};

	//		//Add new patch to grid
	//		if (!notBoundary) _grid.addPatch(bocoName, boco);		
	//	}

	//	//Make all boundary data consistent
	//	//_grid.ConstructAndCheckPatches();

	//	//Close CGNS file
	//	_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Closing CGNS file...");
	//	cg_close(cgfile);

	//	//Synchronize
	//	MPI_Barrier(_comm);

	//	return TURBO_OK;
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
		real_t* tpwgts = new real_t[ncon*nparts];
		for (int i = 0; i<nparts; i++) tpwgts[i] = 1.0/nparts;

		/* An array of size ncon that is used to specify the imbalance tolerance for each vertex weight, with 1
		being perfect balance and nparts being perfect imbalance. A value of 1.05 for each of the ncon
		weights is recommended. */
		real_t* ubvec = new real_t[ncon]; // imbalance tolerance for each vertex weight,
		ubvec[0] = 1.02;	// recomended default value 

		//Algoritm options for displaing information
		idx_t options[3];
		options[0] = 0; //No information and default values

		idx_t edgecut; // Upon successful completion, the number of edges that are cut by the partitioning is written to this parameter.

		/* This is an array of size equal to the number of locally-stored vertices. Upon successful completion the
		partition vector of the locally-stored vertices is written to this array. (See discussion in Section 4.2.4). */
		std::vector<idx_t> part(_grid.nProperCells);

		//Call partitioning function		
		MPI_Comm _comm = _parallelHelper.getComm();
		int result = ParMETIS_V3_PartKway(&_grid.vdist[0], &_grid.xadj[0], &_grid.adjncy[0], vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, &part[0], &_comm);			
		if (result != METIS_OK) {
			_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "ParMETIS_V3_PartKway failed.");
			return TURBO_OK;
		};							

		//Gather partitioning on every processor
		std::vector<int> recvcounts(_nProcessors);
		std::vector<int> displs(_nProcessors);
		displs[0] = 0;
		for (int i = 0; i<_nProcessors; i++) {
			recvcounts[i] = _grid.vdist[i+1] - _grid.vdist[i];
			if (i>0) displs[i] = displs[i-1] + recvcounts[i-1];
		};

		_grid.cellsPartitioning.resize(_grid.nProperCells);
		MPI_Allgatherv(&part[0], _grid.nProperCells, IDX_T, &_grid.cellsPartitioning[0], &recvcounts[0], &displs[0], IDX_T, _comm);

		//Otput result information
		std::ostringstream msg;
		msg<<"Partitioning edgecut = "<<edgecut;
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, msg.str());		
		
		//Add dummy cells
		for (int i = _grid.nProperCells; i<_grid.nCells; i++) {
			Cell* cell = _grid.localCells[i];
			int neigbour = cell->NeigbourCells[0];
			int p = _grid.cellsPartitioning[neigbour];
			_grid.cellsPartitioning.push_back(p);
		};

		//Extract local cells indexes		
		_grid.localCellIndexes.clear();
		for (int i = 0; i<_grid.nCells; i++) {
			if ((_grid.cellsPartitioning[i] == _rank) && (!_grid.Cells[i].IsDummy)) _grid.localCellIndexes.push_back(i);
		};		
		_grid.nCellsLocal = _grid.localCellIndexes.size(); //Without dummy cells
		for (int i = 0; i<_grid.nCells; i++) {
			if ((_grid.cellsPartitioning[i] == _rank) && (_grid.Cells[i].IsDummy)) _grid.localCellIndexes.push_back(i);
		};		

		//Otput result
		msg.str(std::string());
		msg<<"Number of local cells = "<<_grid.nCellsLocal;
		_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());

		//Free memory
		free(tpwgts);
		free(ubvec);									

		//Synchronize		
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Grid partitioning finished");	

		return TURBO_OK;
	}; 

	//Load required geometric grid data to appropriate data structures
	//Provided we have pationed grid
	turbo_errt GenerateGridGeometry() {
		//Create local cells
		//_grid.localCells.resize(_grid.nCellsLocal);		
		//for (int i = 0; i < _grid.localCellIndexes.size(); i++) {					
		//	_grid.localCells[i] = _grid.Cells[_grid.localCellIndexes[i]];
		//};

		////Update cells geometric properties
		//for (int i = 0; i < _grid.nCellsLocal; i++) {
		//	_grid.ComputeGeometricProperties(_grid.localCells[i]);
		//};

		////Compute dummy cell geometric properties

		////Create faces
		//int faceIndex = 0;
		//std::set<Face> faces;
		//_grid.localFaces.clear();
		//for (int i = 0; i < _grid.nCellsLocal; i++) {
		//	//Generate all faces for the cell
		//	Cell& cell = _grid.localCells[i];
		//	std::vector<Face> newFaces = _grid.ObtainFaces(cell);
		//	for (int j = 0; j<newFaces.size(); j++) {
		//		Face& face = newFaces[j];
		//		std::pair<std::set<Face>::iterator, bool> result = faces.insert(face);
		//		//Add face and save connectivity info
		//		if ( result.second ) {
		//			face.GlobalIndex = faceIndex++;
		//			face.FaceCell_1 = cell.GlobalIndex;
		//			face.isExternal = true;
		//			_grid.localFaces.push_back(face);
		//		} else {
		//			face.FaceCell_2 = cell.GlobalIndex;
		//			face.isExternal = false;
		//		};
		//	};
		//};

		//Generate face geometric properties


		//Synchronize
		_parallelHelper.Barrier();
		_logger.WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished generating local geometry and faces");	
		return TURBO_OK;
	};

	turbo_errt InitParallelExchange() {

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

	turbo_errt ReadInitialConditions() {
		//Synchronize
		//MPI_Barrier(_comm);
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

};

#endif
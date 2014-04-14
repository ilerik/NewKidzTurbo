#include "grid.h"
#include "mpi.h"
#include "logger.h"
#include "parmetis.h"

//Define error types
#define err_t int
#define TURBO_OK 0
#define TURBO_ERROR 0

////Calculation kernel
class Kernel {
public:
	//Logger
	std::string _logfilename;
	Logger _logger;

	//Grid data
	Grid _grid;

	//Task parameters
	int _rungeKuttaOrder;	

	//Parallel run information
	int _nProcessors;
	int _rank; 
	int _comm;	//global communicator
	int _nLocalCells;	

	//Internal storage		
	double* Uconservative;	//Conservative variables
	double* Residual;		//Residual

	//METIS datastructures
	idx_t *vdist;
	idx_t *xadj;
	idx_t *adjncy;		

	//Initialize computational kernel
	err_t Initilize(int *argc, char **argv[]) {
		_comm = MPI_COMM_WORLD;
		MPI_Init(argc, argv);		
		MPI_Comm_size(_comm, &_nProcessors);
		MPI_Comm_rank(_comm, &_rank);

		_logfilename = "kernel.log";
		_logger.InitLogging(_logfilename, _rank);
		return TURBO_OK;
	};

	//Finilize kernel work
	err_t Finalize() {
		MPI_Finalize();	
		delete[] vdist;
		delete[] xadj;
		delete[] adjncy;		
		return TURBO_OK;
	};

	//Load grid topology info
	err_t LoadGridTopologyAndInfo(std::string filename) {		
		//Open cgns file
		_grid.gridInfo.fileLocation = filename;
		int cgfile, mode = CG_MODE_READ;		
		_logger.WriteMessage(GLOBAL, INFORMATION, "Opening cgns file"+filename+"");		
		if (cg_open (filename.c_str(), mode, &cgfile)) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Cannot open grid file.");		
			cg_error_exit();
		};
				
		//Determine the number of bases in the grid. This example assumes 
		//one base. However it is allowed to have multiple bases.    
		int nBases;
		if(cg_nbases(cgfile, &nBases)!= CG_OK) {			
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Cannot read number of bases.");		
			cg_error_exit();
		};
		if (nBases != 1) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Cannot process grid. Must be one base.");
			return 1;
		};
		int base = 1;	
		_grid.gridInfo.MainBaseIndex = base;

		//Check the cell and physical dimensions of the base.
		int physDim;
		int cellDim;
		char cgnsName[255];
		if(cg_base_read(cgfile, base, cgnsName, &cellDim, &physDim) != CG_OK) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Cannot read base info.");
			cg_error_exit();
			return 1;
		};		
		_grid.gridInfo.MainBaseName = std::string(cgnsName);
		_grid.gridInfo.CellDimensions = cellDim;
		_grid.gridInfo.GridDimensions = physDim;


		//Read the number of zones in the grid.
		//This example assumes one zone.
		int nZones;
		if(cg_nzones(cgfile, base, &nZones) != CG_OK) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Cannot read number of zones.");
			cg_error_exit();     
		};
		if(nZones != 1) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "This CGNS loader assumes one zone.");
			return 1;
		}
		int zone = 1;
		_grid.gridInfo.MainZoneIndex = 1;

		//Check the zone type. This should be Unstructured.
		ZoneType_t zoneType;
		if(cg_zone_type(cgfile, base, zone, &zoneType) != CG_OK) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Cannot read zone type.");
			cg_error_exit(); 		
		};
		if(zoneType != Unstructured) {			
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Unstructured zone expected.");
			return 1;
		};

		//Determine the number of vertices and volume elements in this */
		//zone (and thus in the grid, because one zone is assumed). */
		char zoneName[255];
		cgsize_t sizes[3];
		if(cg_zone_read(cgfile, base, zone, zoneName, sizes) != CG_OK) cg_error_exit();		
		int nVertices = sizes[0];
		int nCells = sizes[1];

		//Save to grid info
		_grid.gridInfo.MainZoneName = std::string(zoneName);
		_grid.gridInfo.nNodes = nVertices;
		_grid.gridInfo.nCells = nCells;

		//Determine the number and names of the coordinates.
		int nCoords;
		if(cg_ncoords(cgfile, base, zone, &nCoords) != CG_OK)
			cg_error_exit();	

		char name[255];
		DataType_t dataType;
		if(cg_coord_info(cgfile, base, zone, 1, &dataType, name) != CG_OK)
			cg_error_exit();			

		//Load shared grid data 
		//Connectivity info, etc.

		//Read nodes coordinates
		/* Read the x-coordinates. The y and z-coordinates can be read */
		/* similarly. Just replace CoordinateX by CoordinateY and */
		/* CoordinateZ respectively. This assumes Cartesian coordinates */
		/* in double precision. Note that CGNS starts the numbering at */
		/* 1 even if C is used. */
	
		int one = 1;	
		std::vector<double> coorX(nVertices);
		if(cg_coord_read(cgfile, base, zone, "CoordinateX", RealDouble, &one,
						&nVertices, &coorX[0]) != CG_OK) cg_error_exit();	
	
		//Check if we have Y coordinate
		std::vector<double> coorY(nVertices);
		for (int i = 0; i<nVertices; i++) coorY[i] = 0;
		if (nCoords > 1) {
			if(cg_coord_read(cgfile, base, zone, "CoordinateY", RealDouble, &one,
						&nVertices, &coorY[0]) != CG_OK) cg_error_exit();	
		};
	
		//Check if its not 3D case and make all Z equal zero
		std::vector<double> coorZ(nVertices);
		for (int i = 0; i<nVertices; i++) coorZ[i] = 0;
		if (nCoords > 2) {
			if(cg_coord_read(cgfile, base, zone, "CoordinateZ", RealDouble, &one,
						&nVertices, &coorZ[0]) != CG_OK) cg_error_exit();	
		};

		//Create grid nodes
		_grid.localNodes.resize(nVertices);
		for (int i = 0; i<nVertices; i++) {
			Node newNode;
			newNode.GlobalIndex = i+1;
			newNode.P.x = coorX[i];
			newNode.P.y = coorY[i];
			newNode.P.z = coorZ[i];
			_grid.localNodes[i] = newNode;
		};
		
		//Get sections information		
		int nSections;
		if(cg_nsections(cgfile, base, zone, &nSections) != CG_OK)
			cg_error_exit();

		//Find volume elements section also check all element types
		/* Loop over the number of sections and read the element */
		/* connectivities. As CGNS starts the numbering at 1 the */
		/* for-loop starts at 1 as well. */		
		std::vector<int> conn;
		int nNodes;	
		int nConnNodes;		
		int currentCellNumber = 0;
		_grid.gridInfo.GlobalIndexToNumber.clear();
		_grid.gridInfo.NumberToGlobalIndex.clear();
		std::map<int, std::vector<int> > numberToElementNodes;
		std::map<int, ElementType_t > numberToElementType;

		for(int sec=1; sec<=nSections; sec++)
		{
			/* Determine the element type and set the pointer for the */
			/* connectivity accordingly. */
			char secName[255];
			ElementType_t type;
			cgsize_t eBeg;
			cgsize_t eEnd;
			int nBdry;
			int parentFlag;
			if(cg_section_read(cgfile, base, zone, sec, secName, &type,
							&eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
				cg_error_exit();
			int nElements = (eEnd-eBeg+1);
			
			switch (type)
			{
			case QUAD_4: 
				nNodes = 4;			
				break;	
			case BAR_2:
				nNodes = 2;			
				break;	
			case HEXA_8:
				nNodes = 8;
				break;

			default:				
				_logger.WriteMessage(GLOBAL, FATAL_ERROR, "Unsupported element encountered.");
				return 1;				
			}

			//Create cells if its main section (contains elements with dim = physDim)
			//For now we assume that all volume elements are in one section so following condition must hold
			//TO DO Additional checks
			bool isVolumeElementsSection = false;
			if ((GetElementDimensions(type) == _grid.gridInfo.CellDimensions)) {
				isVolumeElementsSection = true;
				_grid.gridInfo.CellsSections.push_back(sec);
			} else {
				isVolumeElementsSection = false;
				_grid.gridInfo.BoundarySections.push_back(sec);
			};			
			
			//Read section connectivity info			
			nConnNodes = nNodes*nElements;
			conn.resize(nConnNodes); 			
   
			/* Read the connectivity. Again, the node numbering of the */
			/* connectivities start at 1. If internally a starting index */
			/* of 0 is used (typical for C-codes) 1 must be substracted */
			/* from the connectivities read. */
   		
			if(cg_elements_read(cgfile, base, zone, sec, &conn[0], NULL) != CG_OK)
				cg_error_exit();	

			//Add all elements read to approriate structures for cells
			if (!isVolumeElementsSection) continue;			
			for (int ind = eBeg; ind<=eEnd; ind++) {
				_grid.gridInfo.GlobalIndexToNumber[ind] = currentCellNumber;
				_grid.gridInfo.NumberToGlobalIndex[currentCellNumber] = ind;
				numberToElementType[currentCellNumber] = type;
				numberToElementNodes[currentCellNumber].resize(nNodes);
				for (int j = 0; j<nNodes; j++) numberToElementNodes[currentCellNumber][j] = (conn[currentCellNumber*nNodes + j]);
				currentCellNumber++;
			};
		};

		//Distribute cells evenly between processors		
		int vProc = nCells / _nProcessors;
		int vLeft = nCells % _nProcessors;
		vdist = new idx_t[_nProcessors + 1];
		vdist[0] = 0;
		for (int i = 0; i<_nProcessors; i++) {
			vdist[i+1] = vdist[i] + vProc;
			if (vLeft > 0) {
				vdist[i+1]++;
				vLeft--;
			};
		};
		vdist[_nProcessors] = nCells;
		
		//TO DO Generalize for different types of volume elements
		//Convert connectivity info to graph
		_nLocalCells = vdist[_rank + 1] - vdist[_rank];
		idx_t *eptr;
		idx_t *eind;
		idx_t ne = _nLocalCells;
		idx_t nn = 0;
		idx_t ncommon = _grid.gridInfo.CellDimensions;
		idx_t numflag = 0;		

		//Fill in connectivity data structures		
		eptr = new idx_t[ne+1];			
		eptr[0] = 0;		
		int localI = 0;
		for (int i = vdist[_rank]; i < vdist[_rank + 1]; i++)
		{		
			int nCellNodes = numberToElementNodes[i].size();
			nConnNodes += nCellNodes;
			eptr[localI+1] = eptr[localI] + nCellNodes;			
			localI++;
		};
		nn = eptr[ne];
		eind = new idx_t[nn];
		localI = 0;
		for (int i = vdist[_rank]; i < vdist[_rank + 1]; i++)
		{		
			int nCellNodes = numberToElementNodes[i].size();
			for (int j = 0; j<nCellNodes; j++) {			
				eind[eptr[localI] + j] = numberToElementNodes[i][j];
			};
			localI++;
		};

		//Call graph generation function
		if (_nProcessors == 1) {
			int result = METIS_MeshToDual(&ne, &nn, eptr, eind, &ncommon, &numflag, &xadj, &adjncy);			
			if (result != METIS_OK) {
				_logger.WriteMessage(GLOBAL, FATAL_ERROR, "METIS_MeshToDual failed.");
				return 1;
			};
		} else {
			int result = ParMETIS_V3_Mesh2Dual(vdist, eptr, eind, &numflag, &ncommon, &xadj, &adjncy, &_comm);
			if (result != METIS_OK) {
				_logger.WriteMessage(GLOBAL, FATAL_ERROR, "ParMETIS_V3_Mesh2Dual failed.");
				return 1;
			};
		};

		//Free memory
		free(eptr);
		free(eind);

		/* Determine the number of boundary conditions for this zone. */

		int nBocos;
		if(cg_nbocos(cgfile, base, zone, &nBocos) != CG_OK)
			cg_error_exit();

		/* Loop over the number of boundary conditions. */

		for(int boco=1; boco<=nBocos; boco++)
		{
			/* Read the info for this boundary condition. */
			char bocoName[255];
			BCType_t bocoType;
			PointSetType_t ptsetType;
			int normalIndex;
			cgsize_t nBCElem;
			cgsize_t normListFlag;
			DataType_t normDataType;
			int nDataSet;
			if(cg_boco_info(cgfile, base, zone, boco, bocoName, &bocoType, &ptsetType,
							&nBCElem, &normalIndex, &normListFlag, &normDataType,
							&nDataSet) != CG_OK) cg_error_exit();

			/* Read the element ID's. */

			std::vector<cgsize_t> BCElemRead(nBCElem);
			if(cg_boco_read(cgfile, base, zone, boco, &BCElemRead[0], NULL) != CG_OK) cg_error_exit();

			/* And much more to make it fit into the */
			/* internal data structures. */

			//Assign to all faces boundary marker
			bool notBoundary = false;
			//for (int i = 0; i<nBCElem; i++) {
			//	std::set<int> idx = indexToElementNodes[BCElemRead[i]];
			//	if (faces.find(idx) != faces.end()) {
			//		//Element is a face
			//		int globalInd = faces[idx];
			//		_grid.faces[globalInd].BCMarker = boco;
			//		if (!_grid.faces[globalInd].isExternal) throw Exception("Boundary face must be marked isExternal");
			//	} else {
			//		//Its not a face so this is not valid Boundary condition
			//		//??
			//		notBoundary = true;
			//	};
			//};

			//Add new patch to grid
			if (!notBoundary) _grid.addPatch(bocoName, boco);		
		}

		//Make all boundary data consistent
		//_grid.ConstructAndCheckPatches();

		//Close CGNS file
		_logger.WriteMessage(GLOBAL, INFORMATION, "Closing CGNS file...");
		cg_close(cgfile);

		//Synchronize
		MPI_Barrier(_comm);

		return 0;
	};

	//Partion loaded grid
	err_t PartitionGrid() {		
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
		idx_t *part = new idx_t[_nLocalCells];

		//Call partitioning function		
		int result = ParMETIS_V3_PartKway(vdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &_comm);			
		if (result != METIS_OK) {
			_logger.WriteMessage(GLOBAL, FATAL_ERROR, "ParMETIS_V3_PartKway failed.");
			return 1;
		};							

		//Gather partitioning on every processor
		int *recvcounts = new int[_nProcessors];
		int *displs = new int[_nProcessors];
		displs[0] = 0;
		for (int i = 0; i<_nProcessors; i++) {
			recvcounts[i] = vdist[i+1] - vdist[i];
			if (i>0) displs[i] = displs[i-1] + recvcounts[i-1];
		};

		int* cellsPartition = new int[_grid.gridInfo.nCells];
		MPI_Allgatherv(part, _nLocalCells, IDX_T, cellsPartition, recvcounts, displs, IDX_T, _comm);

		_grid.cellsPartitioning.resize(_grid.gridInfo.nCells);
		for (int i = 0; i<_grid.gridInfo.nCells; i++) _grid.cellsPartitioning[i] = cellsPartition[i];

		//Extract local cells indexes
		std::vector<int> localCells;
		for (int i = 0; i<_grid.gridInfo.nCells; i++) {
			if (_grid.cellsPartitioning[i] == _rank) localCells.push_back(i);
		};

		//Free memory
		free(tpwgts);
		free(ubvec);
		free(part);
		free(recvcounts);
		free(displs);
		free(cellsPartition);

		//Otput result information
		std::ostringstream msg;
		msg<<"Partitioning edgecut = "<<edgecut;
		_logger.WriteMessage(GLOBAL, INFORMATION, msg.str());		
		msg.str(std::string());
		msg<<"Number of cells = "<<localCells.size();
		_logger.WriteMessage(LOCAL, INFORMATION, msg.str());

		//Synchronize		
		MPI_Barrier(_comm);
	}; 

	//Load required geometric grid data
	err_t LoadGridGeometry() {
		//Create local cells

		//


		//Synchronize
		MPI_Barrier(_comm);
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
#ifndef TURBO_GridLoading_CGNSReader
#define TURBO_GridLoading_CGNSReader

//Inspired by Coolfluid3 code implementation
//https://github.com/coolfluid/coolfluid3/blob/master/cf3/mesh/CGNS/Reader.cpp

#include "stdafx.h"

#include "grid.h"
#include "cgnslib.h"
#include "logger.h"
#include <string>
#include <vector>
#include "parallelHelper.h"

namespace GridLoading {

#define CGNS_CHAR_MAX 1024
#define CGNS_VERT_IDX 0
#define CGNS_CELL_IDX 1
#define CGNS_BVRT_IDX 2

class CGNSReader {	
private:
	//boost::shared_ptr<Logger> _logger;	
	//boost::shared_ptr<ParallelHelper> _parallelHelper;	
	Logger* _logger;
	ParallelHelper* _parallelHelper;
protected:
	struct CGNS_Size
	{
		int nbCoords;
		int nbVertices;
		int nbElements;
		int nbBoundaryVertices;
	} m_size;
  
	bool m_isCoordinatesCreated;
	bool m_uniqueBase;
  
	struct CGNS_File
	{
		int idx;
		int nbBases;
	} m_file;
  
	struct CGNS_Base
	{
		int idx;
		int cell_dim;
		int phys_dim;
		std::string name;
		bool unique;
		int nbZones;
	} m_base;
  
	struct CGNS_Zone
	{
		int idx;
		bool unique;
		std::string name;
		ZoneType_t type;
		int total_nbVertices;
		int nbVertices[3];
		int nbElements;
		int nbBdryVertices;
		int coord_dim;
		int nbGrids;
		int nbSols;
		int nbSections;
		int nbBocos;		
	} m_zone;
  
	struct CGNS_Section
	{
		int idx;
		bool unique;
		std::string name;
		ElementType_t type;
		int eBegin;
		int eEnd;
		int nbBdry;
		int parentFlag;
		int elemNodeCount;
		int elemDataSize;
		int parentData;
		int elemStartIdx;
		int elemEndIdx;
	} m_section;

	struct CGNS_Boco
	{
		int idx;
		bool unique;
		std::string name;
		BCType_t boco_type;  // e.g. BCDirichlet, BCSubsonicInflow, ...
		PointSetType_t ptset_type; // PointList / PointRange / ElementList / ElementRange
		int nBC_elem;
		int normalIndex;
		int normalListFlag;
		DataType_t normalDataType;
		int nDataSet;
	} m_boco; 	

	//Internal numbering mapping from CGNS index to ElementType_t and list of element nodes
	std::unordered_map<int, ElementType_t> elementType;
	std::unordered_map<int, std::vector<int> > elementNodes;

	//Cells information
	int nCells; //Total number of cells in grid
	std::unordered_set<int> cells;//CGNS indexes of cell elements
	std::unordered_map<int, int> _cellCToG; // CGNS index to global index
	std::unordered_map<int, int> _cellGToC; // global index to CGNS index

	//Faces information
	std::unordered_set<int> faces;//CGNS indexes of face elements	
	std::unordered_map<int, int> _faceCToG; // CGNS index to global index
	std::unordered_map<int, int> _faceGToC; // global index to CGNS index

	//Boundary conditions
	int nbBCs;
	std::vector<std::string> bcNames;
	std::vector<std::unordered_set<int>> bcElements;

	//Grid vertices information
	int nNodes; //Total number of nodes in grid
	std::unordered_map<int, int> _nodeCToG; // CGNS index to global index
	std::unordered_map<int, int> _nodeGToC; // global index to CGNS index
	std::vector<Node> nodes; //List of nodes 

private:
	void open_file(std::string filename, int mode = CG_MODE_READ) {
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Opening cgns file"+filename+"");				
		CALL_CGNS(cg_open(filename.c_str(), mode, &m_file.idx));
	};

	void read_coordinates()
	{
		// read coordinates
		int one = 1;
		std::vector<double> xCoord;
		std::vector<double> yCoord;
		std::vector<double> zCoord;
		switch (m_zone.coord_dim)
		{
		case 3:
			zCoord.resize(m_zone.total_nbVertices);
			CALL_CGNS(cg_coord_read(m_file.idx,m_base.idx,m_zone.idx, "CoordinateZ", RealDouble, &one, &m_zone.total_nbVertices, &zCoord[0]));
		case 2:
			yCoord.resize(m_zone.total_nbVertices);
			CALL_CGNS(cg_coord_read(m_file.idx,m_base.idx,m_zone.idx, "CoordinateY", RealDouble, &one, &m_zone.total_nbVertices, &yCoord[0]));
		case 1:
			xCoord.resize(m_zone.total_nbVertices);
			CALL_CGNS(cg_coord_read(m_file.idx,m_base.idx,m_zone.idx, "CoordinateX", RealDouble, &one, &m_zone.total_nbVertices, &xCoord[0]));
		}
		
		for (int i=0; i<m_zone.total_nbVertices; ++i)
		{
			Node newNode;
			switch (m_zone.coord_dim)
			{
				case 3:
					newNode.P.z = zCoord[i];
				case 2:
					newNode.P.y = yCoord[i];
				case 1:
					newNode.P.x = xCoord[i];							
			}
			_nodeCToG[i+1] = nNodes;
			_nodeGToC[nNodes] = i+1;
			newNode.GlobalIndex = nNodes++;
			nodes.push_back(newNode);
		}
		
	}

	void read_section() {
		char section_name_char[CGNS_CHAR_MAX];

		// read section information
		CALL_CGNS(cg_section_read(m_file.idx, m_base.idx, m_zone.idx, m_section.idx, section_name_char, &m_section.type,
								&m_section.eBegin, &m_section.eEnd, &m_section.nbBdry, &m_section.parentFlag));
		m_section.name=section_name_char;

		// replace whitespace by underscore
		boost::algorithm::replace_all(m_section.name," ","_");
		boost::algorithm::replace_all(m_section.name,".","_");
		boost::algorithm::replace_all(m_section.name,":","_");
		boost::algorithm::replace_all(m_section.name,"/","_");	

		//elements stored in section
		int nElements = m_section.eEnd - m_section.eBegin + 1;

		if (m_section.type == CGNS_ENUMV( MIXED )) // Different element types, Can also be faces
		{		
			int connDataSize;				
			int parentData;
			CALL_CGNS(cg_ElementDataSize(m_file.idx, m_base.idx, m_zone.idx, m_section.idx, &connDataSize));
			std::vector<int> connectivityInfo(connDataSize);
			CALL_CGNS(cg_elements_read(m_file.idx, m_base.idx, m_zone.idx, m_section.idx, &connectivityInfo[0], &parentData));			
						
			//Parse each element connectivity info
			int connP = 0;			
			for (int ind = m_section.eBegin; ind<=m_section.eEnd; ind++) {
				//Check element type
				ElementType_t eType = static_cast<ElementType_t>(connectivityInfo[connP++]);
				if (!IsSupprotedCGNSEelementType(eType)) {
					_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unsupported element encountered.");
					exit(0);	
				};				

				//Determine if it is a cell or a face
				bool isCell = GetCGNSElementDimensions(eType) == m_base.cell_dim;
				if (isCell) {
					cells.insert(ind);
				} else {
					faces.insert(ind);
				};

				//Store connectivity info
				int elememtNodeCount = 0;
				CALL_CGNS(cg_npe(eType, &elememtNodeCount)); //obtain number of nodes for given element type
				elementType[ind] = eType;
				elementNodes[ind] = std::vector<int>(elememtNodeCount);
				for (int i = 0; i<elememtNodeCount; i++) {
					elementNodes[ind][i] = connectivityInfo[connP++];
				};				
			};			
		} // if mixed
		else // Single element type in this section
		{
			//Check element type			
			if (!IsSupprotedCGNSEelementType(m_section.type)) {
				_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Unsupported element encountered.");
				exit(0);	
			};

			//Determine if it is a cell or a face
			bool isCell = GetCGNSElementDimensions(m_section.type) == m_base.cell_dim;

			//obtain number of nodes for given element type
			int elememtNodeCount = 0;
			CALL_CGNS(cg_npe(m_section.type, &elememtNodeCount)); 
				
			//Read all connectivity
			int connDataSize;				
			int parentData;
			CALL_CGNS(cg_ElementDataSize(m_file.idx, m_base.idx, m_zone.idx, m_section.idx, &connDataSize));
			std::vector<int> connectivityInfo(connDataSize);
			CALL_CGNS(cg_elements_read(m_file.idx, m_base.idx, m_zone.idx, m_section.idx, &connectivityInfo[0], &parentData));			
						
			//Parse each element connectivity info
			int connP = 0;			
			for (int ind = m_section.eBegin; ind<=m_section.eEnd; ind++) {				
				//Add to apropriate list
				if (isCell) {
					cells.insert(ind);
				} else {
					faces.insert(ind);
				};				

				//Store connectivity info				
				elementType[ind] = m_section.type;						
				elementNodes[ind] = std::vector<int>(elememtNodeCount);
				for (int i = 0; i<elememtNodeCount; i++) {
					elementNodes[ind][i] = connectivityInfo[connP++];
				};		
			};	//for ind	
		}; //if single
	};

	void read_boco() {
		// Read the info for this boundary condition.
		char boco_name_char[CGNS_CHAR_MAX];
		CALL_CGNS(cg_boco_info(m_file.idx, m_base.idx, m_zone.idx, m_boco.idx, boco_name_char, &m_boco.boco_type, &m_boco.ptset_type,
					&m_boco.nBC_elem, &m_boco.normalIndex, &m_boco.normalListFlag, &m_boco.normalDataType, &m_boco.nDataSet));
		m_boco.name = boco_name_char;
  
		// replace whitespace by underscore		
		boost::algorithm::replace_all(m_boco.name," ","_");
		boost::algorithm::replace_all(m_boco.name,".","_");
		boost::algorithm::replace_all(m_boco.name,":","_");
		boost::algorithm::replace_all(m_boco.name,"/","_");

		// Read boundary condition types

		// Read the element ID's
		int* boco_elems = new int [m_boco.nBC_elem];
		void* NormalList(NULL);
		CALL_CGNS(cg_boco_read(m_file.idx, m_base.idx, m_zone.idx, m_boco.idx, boco_elems, NormalList));

		//Process boundary elements dependion on representation type
		switch (m_boco.ptset_type)
		{
		case ElementRange : // all bc elements are within a range given by 2 global element numbers
		{			
			bcNames.push_back(m_boco.name);
			bcElements.push_back(std::unordered_set<int>());
			for (int ind=boco_elems[0]; ind<=boco_elems[1]; ind++)
			{
				bcElements[nbBCs].insert(ind);
			};
			nbBCs++;
			break;
		}
		case ElementList : // all bc elements are listed as global element numbers
		{
			bcNames.push_back(m_boco.name);
			bcElements.push_back(std::unordered_set<int>());						
			for (int i=0; i<m_boco.nBC_elem; ++i)
			{
				int ind = boco_elems[i];        			
				bcElements[nbBCs].insert(ind);
			}
			nbBCs++;
			break;
		}
		case PointRange :  // see default
		case PointList :   // see default
		default :
			_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR,"CGNS: no boundary with pointset_type (PointRange or PointList) supported."); 
			exit(0);
		};
	};

	void read_zone() {
		// get zone type (CGNS_ENUMV( Structured ) or CGNS_ENUMV( Unstructured ))
		CALL_CGNS(cg_zone_type(m_file.idx,m_base.idx,m_zone.idx,&m_zone.type));

		//// For now only CGNS_ENUMV( Unstructured ) and CGNS_ENUMV( Structured ) zone types are supported
		if (m_zone.type != CGNS_ENUMV( Unstructured ))
		{
			_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Only CGNS_ENUMV( Unstructured ) zone type is supported");
			exit(0);
		};

		// Read zone size and name
		if (m_zone.type == CGNS_ENUMV( Unstructured ))
		{
			cgsize_t size[3][1];
			char zone_name_char[CGNS_CHAR_MAX];
			CALL_CGNS(cg_zone_read(m_file.idx,m_base.idx,m_zone.idx,zone_name_char,size[0]));
			m_zone.name = zone_name_char;
			boost::algorithm::replace_all(m_zone.name," ","_");
			boost::algorithm::replace_all(m_zone.name,".","_");
			boost::algorithm::replace_all(m_zone.name,":","_");
			boost::algorithm::replace_all(m_zone.name,"/","_");
			m_zone.total_nbVertices = size[CGNS_VERT_IDX][0];
			m_zone.nbElements       = size[CGNS_CELL_IDX][0];
			m_zone.nbBdryVertices   = size[CGNS_BVRT_IDX][0];

			// get the number of grids
			CALL_CGNS(cg_ngrids(m_file.idx,m_base.idx,m_zone.idx,&m_zone.nbGrids));
			// nb coord dims
			CALL_CGNS(cg_ncoords(m_file.idx,m_base.idx,m_zone.idx, &m_zone.coord_dim));
			// find out number of solutions
			CALL_CGNS(cg_nsols(m_file.idx,m_base.idx,m_zone.idx,&m_zone.nbSols));
			// find out how many sections
			CALL_CGNS(cg_nsections(m_file.idx,m_base.idx,m_zone.idx,&m_zone.nbSections));
			m_section.unique = m_zone.nbSections == 1 ? true : false;
			// find out number of BCs that exist under this zone
			CALL_CGNS(cg_nbocos(m_file.idx,m_base.idx,m_zone.idx,&m_zone.nbBocos));
			m_boco.unique = m_zone.nbBocos == 1 ? true : false;
			// Add up all the nb elements from all sections
			//m_zone.total_nbElements = get_total_nbElements();			

			// read coordinates in this zone
			for (int i=1; i<=m_zone.nbGrids; ++i);
				read_coordinates();

			// read sections (or subregions) in this zone			
			for (m_section.idx=1; m_section.idx<=m_zone.nbSections; ++m_section.idx)
				read_section();
		
			// read boundaryconditions (or subregions) in this zone
			nbBCs = 0;
			bcNames.clear();
			bcElements.clear();
			for (m_boco.idx=1; m_boco.idx<=m_zone.nbBocos; ++m_boco.idx)
				read_boco();		

			// check that each face element has BC attached
			std::unordered_set<int> boundaryFaces;
			for (int i = 0; i<nbBCs; i++) {
				for (std::unordered_set<int>::iterator it = bcElements[i].begin(); it != bcElements[i].end(); it++) {
					boundaryFaces.insert(*it);
				};
			};

			bool isConsistent = true;			
			std::unordered_set<int>::iterator itbf = boundaryFaces.begin();
			if (faces.size() != boundaryFaces.size()) isConsistent = false;
			while (isConsistent) {
				if (itbf == boundaryFaces.end()) break;					
				if (faces.find(*itbf) == faces.end()) {
					isConsistent = false;
					break;
				};				
				itbf++;				
			};

			if (!isConsistent) {
				_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "CGNS: Not all face elements from grid have their index attached to boundary condition.");
				exit(0);
			};
		//	//m_mesh->geometry_fields().update_structures();
		//	//read_flowsolution();
		};
	};

	void read_base()
	{
		// get the name, dimension and physical dimension from the base
		char base_name_char[CGNS_CHAR_MAX];
		CALL_CGNS(cg_base_read(m_file.idx,m_base.idx,base_name_char,&m_base.cell_dim,&m_base.phys_dim));
		m_base.name=base_name_char;
		boost::algorithm::replace_all(m_base.name," ","_");
		boost::algorithm::replace_all(m_base.name,".","_");
		boost::algorithm::replace_all(m_base.name,":","_");
		boost::algorithm::replace_all(m_base.name,"/","_"); 	  
  	  
		// check how many zones we have
		CALL_CGNS(cg_nzones(m_file.idx,m_base.idx,&m_base.nbZones));
		m_zone.unique = m_base.nbZones == 1 ? true : false;

		//We assume that only one zone allowed
		if (m_base.nbZones != 1) {
			_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "CGNSReader error: we dont support multiple zones");
		};

		// Read every zone in this base
		for (m_zone.idx = 1; m_zone.idx<=m_base.nbZones; ++m_zone.idx)
			read_zone();

		//
		m_zone.idx = 1;
	};

public:

	//Constructor
	CGNSReader() {	
	};

	void Init(Logger& logger, ParallelHelper& helper)
	{
		//_logger = boost::shared_ptr<Logger>(&logger);
		//_parallelHelper = boost::shared_ptr<ParallelHelper>(&helper);
		_logger = &logger;
		_parallelHelper = &helper;
	};

	void Finalize() {		
	};

	//Read CGNS mesh from file
	Grid LoadGrid(std::string fname) {
		//Initialize
		cells.clear();
		faces.clear();
		nCells = 0;
		nNodes = 0;
		
		//Create new grid object
		Grid grid;

		//Open file
		grid.gridInfo.fileLocation = fname;			
		open_file(fname);

		//Read number of bases
		CALL_CGNS(cg_nbases(m_file.idx, &m_file.nbBases));		

		//We assume that only one base allowed
		if (m_file.nbBases != 1) {
			_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "CGNSReader error: we dont support multiple bases");
		};

		//Read base info and all other grid info recursivelly
		m_base.idx = 1;
		read_base();		

		//Close file
		CALL_CGNS(cg_close(m_file.idx));

		//Fill in grid info structure
		//TODO		
		grid.gridInfo.CellDimensions = m_base.cell_dim;

		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good.");	

		//Copy nodes information
		grid.localNodes = nodes;

		//Numerate cells
		grid.nProperCells = cells.size();
		_cellCToG.clear();
		_cellGToC.clear();
		int globalIndex = 0;		
		for (int cgnsIndex : cells) {
			_cellCToG[cgnsIndex] = globalIndex;
			_cellGToC[globalIndex++] = cgnsIndex;
		};

		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good.");	

		//Dummy cells\boundary faces
		grid.nDummyCells = faces.size();
		_faceCToG.clear();
		_faceGToC.clear();
		for (int cgnsIndex : faces) {
			_cellCToG[cgnsIndex] = globalIndex;
			_cellGToC[globalIndex] = cgnsIndex;
			_faceCToG[cgnsIndex] = globalIndex;
			_faceGToC[globalIndex++] = cgnsIndex;
		};

		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good.");	

		//Initially distribute vertices uniformly over processors		
		nCells = grid.nCells = globalIndex;
		int nProcessors = _parallelHelper->getProcessorNumber();
		int vProc = nCells / nProcessors;
		int vLeft = nCells % nProcessors;
		grid.vdist.resize(nProcessors + 1);
		grid.vdist[0] = 0;
		for (int i = 0; i<nProcessors; i++) {
			grid.vdist[i+1] = grid.vdist[i] + vProc;
			if (vLeft > 0) {
				grid.vdist[i+1]++;
				vLeft--;
			};
		};
		grid.vdist[nProcessors] = nCells;
				
		//Convert connectivity info to graph
		int rank = _parallelHelper->getRank();		
		int nLocalCells = grid.vdist[rank + 1] - grid.vdist[rank];
		std::vector<idx_t> eptr;
		std::vector<idx_t> eind;
		idx_t ne = nLocalCells;
		idx_t nn = 0;
		idx_t ncommon = grid.gridInfo.CellDimensions;
		idx_t numflag = 0;	

		//Restore grid connectivity
		//Fill in connectivity data structures		
		eptr.resize(ne+1);
		eptr[0] = 0;		
		int localI = 0;
		for (int i = grid.vdist[rank]; i < grid.vdist[rank + 1]; i++)
		{		
			int elementID = _cellGToC[i];
			int nCellNodes = elementNodes[elementID].size();			
			eptr[localI+1] = eptr[localI] + nCellNodes;			
			//For dummy cells 
			if (i >= grid.nProperCells) {
				eptr[localI+1] += 1;
			};
			localI++;
		};
		nn = eptr[ne];
		eind.resize(nn);
		localI = 0;
		int fakeVertexID = grid.localNodes.size();
		for (int i = grid.vdist[rank]; i < grid.vdist[rank + 1]; i++)
		{		
			int elementID = _cellGToC[i];
			int nCellNodes = elementNodes[elementID].size();
			for (int j = 0; j<nCellNodes; j++) {			
				eind[eptr[localI] + j] = elementNodes[elementID][j];
			};
			//For dummy cells 
			if (i >= grid.nProperCells) {
				eind[eptr[localI] + nCellNodes] = fakeVertexID++;
			};
			localI++;
		};

		//Call graph generation function
		idx_t* _xadj;
		idx_t* _adjncy;		
		MPI_Comm _comm = _parallelHelper->getComm();
		if (nProcessors != 1) {
			int result = ParMETIS_V3_Mesh2Dual(&grid.vdist[0], &eptr[0], &eind[0], &numflag, &ncommon, &_xadj, &_adjncy, &_comm);
			if (result != METIS_OK) {
				_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "ParMETIS_V3_Mesh2Dual failed.");
				exit(0);			
			};
		} else {
			//We use serial function (TO DO parallelize)
			int result = METIS_MeshToDual(&ne, &nn, &eptr[0], &eind[0], &ncommon, &numflag, &_xadj, &_adjncy);			
			if (result != METIS_OK) {				
				_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "METIS_MeshToDual failed.");
				exit(0);
			};
		};

		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good transform.");	
		
		//Make local copies of _xadj amd _adjncy
		std::vector<int> xadjLocal;
		xadjLocal.resize(grid.vdist[rank+1] - grid.vdist[rank]);
		for (int i = 0; i<xadjLocal.size(); i++) xadjLocal[i] = _xadj[i];
		int nEdgesLocal = _xadj[grid.vdist[rank+1] - grid.vdist[rank]];
		std::vector<int> adjncyLocal;
		adjncyLocal.resize(nEdgesLocal);
		for (int i = 0; i<adjncyLocal.size(); i++) adjncyLocal[i] = _adjncy[i];		
		
		
		//Free memory
		METIS_Free(_xadj);
		METIS_Free(_adjncy);

		//Gather connectivity information		
		std::vector<int> nEdges;
		_parallelHelper->AllgatherCounts(nEdgesLocal, nEdges);		

		std::ostringstream msg;
		msg<<"vdist[] = \n";
		for (int i = 0; i<=nProcessors; i++) msg<<grid.vdist[i]<<" ";
		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		msg.clear();
		msg<<"nEdges[] = \n";
		for (int i = 0; i<nProcessors; i++) msg<<nEdges[i]<<" ";
		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

		//Gather xadj structure
		std::vector<int> xadj;
		std::vector<int> adjncy;
		std::vector<int> counts;
		for (int i = 0; i<nProcessors; i++) counts.push_back(grid.vdist[i+1] - grid.vdist[i]);				
		_parallelHelper->Allgatherv(xadjLocal, counts, xadj);

		//Further renumbering for serial CSR
		int k = 1;	
		int totalEdges = 0;
		for (int i = 0; i<xadj.size(); i++) {
			if (i >= grid.vdist[k]) {
				totalEdges += nEdges[k-1];
				k++;				
			};
			xadj[i] += totalEdges;			
		};
		totalEdges += nEdges[nEdges.size() - 1];

		//Add total number of edges to xadj
		/*int totalEdges = 0;
		for (int i = 0; i<nEdges.size(); i++) totalEdges += nEdges[i];*/
		xadj.push_back(totalEdges);
		
		//Gather adjncy structure		
		_parallelHelper->Allgatherv(adjncyLocal, nEdges, adjncy);

		//Create basic cells structures (only type and nodes and connectivity)
		grid.nCells = nCells;
		grid.Cells.resize(nCells);
		for (int i = 0; i<grid.nCells; i++) {
			Cell newCell;			
			int elementID = _cellGToC[i];
			newCell.GlobalIndex = i;			
			newCell.CGNSType = elementType[elementID];
			//Fill in nodes
			newCell.Nodes.clear();
			int nCellNodes = elementNodes[elementID].size();
			for (int j = 0; j<nCellNodes; j++) {							
				newCell.Nodes.push_back(_nodeCToG[elementNodes[elementID][j]]);
			};
			//Set dumminess
			newCell.IsDummy = (i >= grid.nProperCells);						

			//Fill in neighbours
			newCell.NeigbourCells.clear();				
			for (int j = xadj[i]; j<xadj[i+1]; j++) {
				int neighbour = adjncy[j];
				int neighbourCGNS = _cellGToC[neighbour];
				std::vector<int> nnodes = elementNodes[neighbourCGNS];
				newCell.NeigbourCells.push_back(neighbour);
			};			
			grid.Cells[i] = newCell;
		};	
		
		//Modify dummy cells BCMarkers and create patches for grid
		for (int BCMarker = 0; BCMarker < bcNames.size(); BCMarker++) {
			grid.addPatch(bcNames[BCMarker], BCMarker);
			for (int faceCGNSIndex : bcElements[BCMarker]) {	
				int faceIndex = _faceCToG[faceIndex];
				grid.Cells[faceIndex].BCMarker = BCMarker;
				grid.patches[BCMarker].addFace(faceIndex);
				for (int nodeIndex : grid.Cells[faceIndex].Nodes) grid.patches[BCMarker].addNode(nodeIndex);
			};
		};

		//Now exclude dummy cells information from graph structure
		grid.xadj.clear();
		grid.adjncy.clear();

		//Modify vdist structure (exclude dummy)
		vProc = grid.nProperCells / nProcessors;
		vLeft = grid.nProperCells % nProcessors;
		grid.vdist.resize(nProcessors + 1);
		grid.vdist[0] = 0;
		for (int i = 0; i<nProcessors; i++) {
			grid.vdist[i+1] = grid.vdist[i] + vProc;
			if (vLeft > 0) {
				grid.vdist[i+1]++;
				vLeft--;
			};
		};
		grid.vdist[nProcessors] = grid.nProperCells;
		grid.nCellsLocal = grid.vdist[rank+1] - grid.vdist[rank];

		//Fill in xadj and adjncy (exclude dummy)
		int currentInd = 0;
		grid.xadj.push_back(currentInd);
		for (int i = grid.vdist[rank]; i < grid.vdist[rank + 1]; i++)
		{		
			int index = i;
			int elementID = _cellGToC[i];
			int start = xadj[index];
			int end = xadj[index+1];
			for (int j = start; j<end; j++) {
				int nIndex = adjncy[j];
				if (nIndex >= grid.nProperCells) continue; //Skip dummy
				grid.adjncy.push_back(nIndex);
				currentInd++;
			};
			grid.xadj.push_back(currentInd);
		};					
		
		msg.clear();
		msg<<"vdist[] = \n";
		for (int i = 0; i<=nProcessors; i++) msg<<grid.vdist[i]<<" ";
		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());			
		_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "Good read.");			

		//Synchronize
		_parallelHelper->Barrier();
		
		return grid;
	};

	//Read physical field from opened CGNS file
	void ReadField(Grid& grid, std::string solutionName, std::string fieldName, std::vector<double>& variable) {		
		//Open file		 
		open_file(grid.gridInfo.fileLocation, CG_MODE_READ);

		//Check if variable array lenght equals number of local cells
		if (variable.size() != grid.nCellsLocal) {
			_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::FATAL_ERROR, "Number of field variebles doesn't equal to number of cells.");
			exit(0);
		};

		//Check if such solution exist
		bool solutionFound = false;
		int nSols = 0;
		//Determine number of solutions
		CALL_CGNS(cg_nsols(m_file.idx, m_base.idx, m_zone.idx, &nSols));

		//Find solution with given name
		int S = 0;
		for (S = 1; S<=nSols; S++) {
			//Read solution name and grid location
			char solName[256];
			GridLocation_t gridLocation;
			CALL_CGNS(cg_sol_info(m_file.idx, m_base.idx, m_zone.idx, S, solName, &gridLocation));

			//Check if it's needed solution
			if (std::string(solName) == solutionName) {
				solutionFound = true;
				if (gridLocation != CellCenter) {
					_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Could not read solution with grid location not CellCenter.");
					_parallelHelper->Barrier();
					return;
				};
				break;
			};
		};

		//If we didnt find solution
		if (!solutionFound) {
			_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Could not find solution with such name.");
			_parallelHelper->Barrier();
			return;
		};

		//Find field with given name
		bool fieldFound = false;
		DataType_t datatype;
		int nFields;
		CALL_CGNS(cg_nfields(m_file.idx, m_base.idx, m_zone.idx, S, &nFields));
		for (int F = 1; F<=nFields; F++) {
			//Read solution name and grid location
			char fieldname[256];			
			CALL_CGNS(cg_field_info(m_file.idx, m_base.idx, m_zone.idx, S, F, &datatype, fieldname));

			//Check if it's needed field
			if (std::string(fieldname) == fieldName) {
				fieldFound = true;	
				break;
			};
		};

		//If we didnt find field
		if (!solutionFound) {
			_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Could not find field with such name.");
			_parallelHelper->Barrier();
			return;
		};

		//Read cell values for each cell
		for (int i = 0; i<grid.nCellsLocal; i++) {
			int cIndex = grid.localCells[i]->GlobalIndex;
			int cgnsIndex = _cellGToC[cIndex];
			cgsize_t range_min[1]; 
			range_min[0] = cgnsIndex;
			cgsize_t range_max[1];
			range_max[0] = cgnsIndex;
			cg_field_read(m_file.idx, m_base.idx, m_zone.idx, S, fieldName.c_str(), datatype, range_min, range_max, &variable[i]);
		};
	
		_parallelHelper->Barrier();		

		//Close file
		CALL_CGNS(cg_close(m_file.idx));
	};

}; //CGNSReader class

}; //GridLoading namespace

#endif

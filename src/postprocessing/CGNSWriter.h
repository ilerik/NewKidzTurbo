#ifndef TURBO_PostProcessing_CGNSWriter
#define TURBO_PostProcessing_CGNSWriter

#include "stdafx.h"

#include "grid.h"
#include "cgnslib.h"
#include "logger.h"
#include "parallelHelper.h"
#include <string>

namespace PostProcessing {

class CGNSWriter {
private:
	//Pointer to logger
	Logger* _logger;
	//Pointer to parrallel interface
	ParallelHelper* _parallelHelper;
protected:
	//Internal cgns tructure info
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
		cgsize_t size[3];
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
public:
	
	//Constructor
	CGNSWriter() {	
	};

	//Initialize required subsystems
	void Init(Logger& logger, ParallelHelper& helper)
	{		
		_logger = &logger;
		_parallelHelper = &helper;
	};

	//Finilize
	void Finalize() {		
	};

	//Create new cgns file or rewrite existing
	void CreateFile(std::string fname) {
		int mode = CG_MODE_WRITE;
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Opening cgns file " + fname + " for writing.");				
		if (_parallelHelper->IsMaster()) {
			CALL_CGNS(cg_open(fname.c_str(), mode, &m_file.idx));
		};
		_parallelHelper->Barrier();
	};

	//Write grid structure to file
	void WriteGridToFile(Grid& grid) {
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Started writing grid.");				
					
		//Gather number of cells on master node
		std::vector<int> nCellsPerProcessor;
		_parallelHelper->GatherCounts(grid.nCellsLocal, nCellsPerProcessor);

		//Cells elements connectivity for current processor
		std::vector<cgsize_t> LocalElements;
		for (int i = 0; i < grid.nCellsLocal; i++) {
			Cell* cell = grid.localCells[i];
			LocalElements.push_back(cell->CGNSType);
			for (int j = 0; j<cell->Nodes.size(); j++) LocalElements.push_back(cell->Nodes[j] + 1); //Adjust by 1 for cgns proper numbering
		};

		//Gather local elements arrays sizes
		std::vector<int> nLocalElementsSize;
		_parallelHelper->GatherCounts(LocalElements.size(), nLocalElementsSize);

		//Gather all elements on master node
		std::vector<cgsize_t> AllElements;
		_parallelHelper->GathervInt(LocalElements, nLocalElementsSize, AllElements);
		

		//Write to cgns database file
		if (_parallelHelper->IsMaster()) {
			//Create base
			m_base.cell_dim = grid.gridInfo.CellDimensions;
			m_base.phys_dim = grid.gridInfo.CellDimensions;
			CALL_CGNS(cg_base_write(m_file.idx, "Base", m_base.cell_dim, m_base.phys_dim, &m_base.idx));

			//Create zone
			m_zone.coord_dim = grid.gridInfo.CellDimensions;
			m_zone.name = "Zone";	
			m_zone.type = Unstructured;
			m_zone.size[0] = grid.localNodes.size(); //NVertex
			m_zone.size[1] = grid.nCells; //NCells
			m_zone.size[2] = 0; //NBoundaryVertex			
			CALL_CGNS(cg_zone_write(m_file.idx, m_base.idx, m_zone.name.c_str(), m_zone.size, m_zone.type, &m_zone.idx));

			//Write nodes coordinates		
			int G;
			//Note that the name "GridCoordinates" is reserved for the original grid and must be the first GridCoordinates_t node to be defined.
			CALL_CGNS(cg_grid_write(m_file.idx, m_base.idx, m_zone.idx, "GridCoordinates", &G));
			//Number of nodes
			int nNodes = grid.localNodes.size();
			std::vector<double> coords(nNodes);
			int C;
			if (m_base.phys_dim >= 1) {
				//Write x cartezian coordinates to buffer
				for (int i = 0; i<nNodes; i++) coords[i] = grid.localNodes[i].P.x;
				//Write coords to file
				CALL_CGNS(cg_coord_write(m_file.idx, m_base.idx, m_zone.idx, RealDouble, "CoordinateX", &coords[0], &C));
			};
			if (m_base.phys_dim >= 2) {
				//Write y cartezian coordinates to buffer
				for (int i = 0; i<nNodes; i++) coords[i] = grid.localNodes[i].P.y;
				//Write coords to file
				CALL_CGNS(cg_coord_write(m_file.idx, m_base.idx, m_zone.idx, RealDouble, "CoordinateY", &coords[0], &C));
			};
			if (m_base.phys_dim >= 3) {
				//Write y cartezian coordinates to buffer
				for (int i = 0; i<nNodes; i++) coords[i] = grid.localNodes[i].P.z;
				//Write coords to file
				CALL_CGNS(cg_coord_write(m_file.idx, m_base.idx, m_zone.idx, RealDouble, "CoordinateZ", &coords[0], &C));
			};			

			//Write cells elements info into main section			
			int start = 1;
			int end = start;
			for (int i = 0; i<nCellsPerProcessor.size(); i++) end += nCellsPerProcessor[i];
			end -= 1;
			int S;
			CALL_CGNS(cg_section_write(m_file.idx, m_base.idx, m_zone.idx, "Cells",  MIXED, start, end, 0, &AllElements[0], &S));
		};
		_parallelHelper->Barrier();
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished writing grid.");				
	};

	//Write physical field to file 
	void WriteField(Grid& grid, std::string fieldName, const std::vector<double>& variable) {		
		_parallelHelper->Barrier();
	};

	//Close file
	void CloseFile() {
		if (_parallelHelper->IsMaster()) {
			CALL_CGNS(cg_close(m_file.idx));
		};
		_parallelHelper->Barrier();
	};
};

}; //namespace

#endif
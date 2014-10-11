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

	std::string fileLocation; //File location

	//Dictionary of created solution nodes to their indexes
	std::map<std::string, int> solutionNameToS;

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
		fileLocation = fname;
		_parallelHelper->Barrier();
	};	

	//Open existing file
	void OpenFile(std::string fname) {		
		//Open file for modification
		int mode = CG_MODE_MODIFY;
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Opening cgns file " + fname + " for writing.");				
		if (_parallelHelper->IsMaster()) {
			CALL_CGNS(cg_open(fname.c_str(), mode, &m_file.idx));
		};
		fileLocation = fname;
		_parallelHelper->Barrier();

		//Obtain information about existing solution nodes
		int nSolutions = 0;		
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Getting information about existing solutions.");				
		if (_parallelHelper->IsMaster()) {
			//Trying to get solutions number
			int ierror = cg_nsols(m_file.idx, m_base.idx, m_zone.idx, &nSolutions );
			if (ierror) {
				//Error while reading asuuming no solutions
				_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Could't get any solutions info. Assuming no solutions exist.");	
			};
		};
		
		//Reading names of all exiting solutions
		char solName[256]; //TO DO remove hardcode (use cgns constant instead)
		solutionNameToS.clear();		
		for (int S = 0; S<nSolutions; S++) {
			if (_parallelHelper->IsMaster()) {				
				GridLocation_t gridLocation;
				CALL_CGNS(cg_sol_info(m_file.idx, m_base.idx, m_zone.idx, S, solName, &gridLocation));
				
				//We dont support other grid locations
				if (gridLocation == GridLocation_t::CellCenter) {					
					solutionNameToS[std::string(solName)] = S; //Add solution name
				};
			};			
		};

		_parallelHelper->Barrier();
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished reading solutions info.");	
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished opening file.");	
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
		

		//Compute total number of cells
		int totalCells = _parallelHelper->SumInt(grid.nCellsLocal);

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
			m_zone.size[1] = totalCells; //NCells
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

	//Write iterative data to base in file
	void WriteIterativeData(double time, int iteration) {
		if (_parallelHelper->IsMaster()) {	
			//Write iterative data node	with one iteration info				
			cgsize_t nSteps = 1;			
			char solutionName[32] = "Solution";
			CALL_CGNS(cg_biter_write(m_file.idx, m_base.idx, "IterativeData", nSteps));
			CALL_CGNS(cg_goto(m_file.idx, m_base.idx, "IterativeData", 0, "end"));
			CALL_CGNS(cg_array_write("TimeValues", DataType_t::RealDouble, 1, &nSteps, &time));			
			CALL_CGNS(cg_ziter_write(m_file.idx, m_base.idx, m_zone.idx, "IterativeData"));
			CALL_CGNS(cg_goto(m_file.idx, m_base.idx, "Zone", 0, "IterativeData", 0, "end"));
			int dimensionVector[2];
			dimensionVector[0] = 32;
			dimensionVector[1] = nSteps;
			CALL_CGNS(cg_array_write("FlowSolutionPointers", DataType_t::Character, 2, dimensionVector, solutionName));			
		};

		//Synchronize
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished writing iterative data node.");				
		_parallelHelper->Barrier();
	};

	//Write new solution node to file
	void WriteSolution(Grid& grid, std::string solutionName) {		
		if (_parallelHelper->IsMaster()) {	
			//Write solution node		
			int S;
			CALL_CGNS(cg_sol_write(m_file.idx, m_base.idx, m_zone.idx, solutionName.c_str(), CellCenter, &S));
			solutionNameToS[ solutionName ] = S;
		};

		//Synchronize
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished writing solution node.");				
		_parallelHelper->Barrier();
	};

	//Write physical field to file 
	void WriteField(Grid& grid, std::string solutionName, std::string fieldName, std::vector<double>& variable) {		
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Started writing field " + fieldName + " to solution" + solutionName + ".");				

		//Check if solution node was created
		int S;
		if (_parallelHelper->IsMaster()) {				
			bool isExistingSolution = (solutionNameToS.find(solutionName) != solutionNameToS.end());		
			if (!isExistingSolution) {
				//Error trying to write solution that wasnt created
				_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::FATAL_ERROR, "Solution with name " + solutionName + " doesn't exist. Cant write data");
				exit(0);
			} else {
				S = solutionNameToS[solutionName];
			};
		};

		//Check if variable array lenght equals number of local cells
		if (variable.size() != grid.nCellsLocal) {
			_logger->WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::FATAL_ERROR, "Number of field variebles doesn't equal to number of cells.");
			exit(0);
		};

		//Gather number of cells on master node
		std::vector<int> nCellsPerProcessor;
		_parallelHelper->GatherCounts(grid.nCellsLocal, nCellsPerProcessor);				

		//Gather all elements on master node
		std::vector<double> AllVariables;
		_parallelHelper->GathervDouble(variable, nCellsPerProcessor, AllVariables);

		if (_parallelHelper->IsMaster()) {			
			//Write flow field data
			int F;
			CALL_CGNS(cg_field_write(m_file.idx, m_base.idx, m_zone.idx, S, RealDouble, fieldName.c_str(), &AllVariables[0], &F));
		};

		//Synchronize
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "Finished writing field " + fieldName + " to solution" + solutionName + ".");				
		_parallelHelper->Barrier();		
	};

	//Close file
	void CloseFile() {
		if (_parallelHelper->IsMaster()) {
			CALL_CGNS(cg_close(m_file.idx));
		};
		_parallelHelper->Barrier();
		_logger->WriteMessage(LoggerMessageLevel::GLOBAL, LoggerMessageType::INFORMATION, "CGNS file closed.");				
	};
};

}; //namespace

#endif
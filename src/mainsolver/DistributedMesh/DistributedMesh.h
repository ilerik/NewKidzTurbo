#ifndef NewKidzTurbo_MainSolver_DistributedMesh_DistributedMeshMesh
#define NewKidzTurbo_MainSolver_DistributedMesh_DistributedMeshMesh

#include "basetypes.h"
#include "alglibmisc.h"
#include "cgnslib.h"
#include "patch.h"
#include "parmetis.h"
#include <map>
#include <set>
#include <algorithm>
#include <cassert>
#include <memory>

#include "ParallelManager.h"

//Grid types
struct Node {	
	int GlobalIndex; //Global index of node
	Vector P;		 //Position	
};

enum FaceType {
	LOCAL,
	BOUNDARY,
	INTERPROCESSOR
};

struct Face {
	int GlobalIndex;
	ElementType_t CGNSType;
	FaceType Type;		//Type of face (boundary, local, interprocessor)
	int FaceCell_1;		//Index of cell 1
	int FaceCell_2;		//Index of cell 2
	Vector FaceCenter;	//Center of face
	Vector FaceNormal;	//Surface normal * square
	double FaceSquare;	//Face square
	std::vector<int> FaceNodes;	//Indexes of all face nodes
	int isExternal;		//If its an external face (FaceCell_2 == dummy cell global index)	
	int BCMarker;		//Marker of boundary condition type	
	std::map<std::string,  void*> Data; //Arbitrary data	

	//Constructors
	Face() {
		BCMarker = 0;
	};

	Face(ElementType_t type, std::vector<int> nodes) {
		CGNSType = type;
		FaceNodes = nodes;
		BCMarker = 0;
	};

	bool operator<( const Face& second ) const
    {
		//Compare length
		if (FaceNodes.size() < second.FaceNodes.size()) return true;
		if (FaceNodes.size() > second.FaceNodes.size()) return false;
		//In case length are equal
		for (int i = 0; i<FaceNodes.size(); i++)
			if (FaceNodes[i] != second.FaceNodes[i]) return FaceNodes[i] < second.FaceNodes[i];
		//Last case they are equal
		return false;
    }

	//Utility function outputs face info
	std::string getInfo() {
		std::ostringstream msg;
		msg << "Face Info :\n";
		msg << "Global Index : " << GlobalIndex << "\n";
		msg << "CGNS Type : " << cg_ElementTypeName(CGNSType) << "\n";
		msg << "Nodes : ";
		for (int ind : FaceNodes) msg << ind << " ";
		msg << "\n";
		msg << "IsExternal = " << isExternal << "\n";
		msg << "FaceNormal = (" << FaceNormal.x << "," << FaceNormal.y << "," << FaceNormal.z << ")\n";
		return msg.str();
	};
};

struct FaceComparer
{
	//Is first less than second
    bool operator()( const Face& first , const Face& second) const
    {
		//Compare length
		if (first.FaceNodes.size() < second.FaceNodes.size()) return true;
		if (first.FaceNodes.size() > second.FaceNodes.size()) return false;
		//In case length are equal		
		for (int i = 0; i<first.FaceNodes.size(); i++)
			if (first.FaceNodes[i] != second.FaceNodes[i]) return first.FaceNodes[i] < second.FaceNodes[i];
		//Last case they are equal
		return false;
    }
};

struct Cell {
public:
	//Common properties
	int GlobalIndex;		//Global cell index
	bool IsDummy;			//Is dummy cell
	int BCMarker;			//Boundary condition marker
	ElementType_t CGNSType; //CGNS cell type	
	std::vector<int> Nodes; //Cell  node indexes
	std::vector<int> NeigbourCells; //Indexes of neighbour cells
	//Local properties
	int CGNSElementIndex;	//CGNS element index	
	std::vector<int> Faces; //Indexes of all cell faces	
	double CellVolume;		// Volume of cell
	Vector CellCenter;		// Center of cell
	double CellHSize;		// Size of cell		
	std::map<std::string,  void*> Data; //Arbitrary data	

	//Utility function outputs cell info
	std::string getInfo() {
		std::ostringstream msg;
		msg << "Cell Info :\n";
		msg << "Global Index : " << GlobalIndex << "\n";
		msg << "CGNS Type : " << cg_ElementTypeName(CGNSType) << "\n";
		msg << "Nodes : ";
		for (int ind : Nodes) msg << ind << " ";
		msg << "\n";	
		msg << "Neighbours : ";
		for (int ind : NeigbourCells) msg << ind << " ";
		msg << "\n";	
		msg << "CellCenter = (" << CellCenter.x << "," << CellCenter.y << "," << CellCenter.z << ")\n";
		return msg.str();
	};
};



//general info
class GridInfo {
public:
	//Basic grid properties and cgns file info
	std::string fileLocation;
	int nCoords;
	int CellDimensions;
	int GridDimensions;	
	int nNodes;
	int nCells;

	//CGNS info
	int MainBaseIndex; //base index (one base assumed)
	std::string MainBaseName; 
	int MainZoneIndex; //zone index (one zone assumed)
	std::string MainZoneName; 
	std::vector<int> CellsSections; // volume elements sections (possible many)
	std::vector<int> BoundarySections; // boundary elements sections (possible many)	

	//raw cgns connectivity info
	std::map<int, std::set<int> > cgnsIndexToElementNodes;	
	std::map<int, ElementType_t > cgnsIndexToElementType;
	std::vector<int> boundaryElements;  //Boundary faces
	std::vector<int> volumeElements;	//Cells

	//Indexes mapping and numeration
	std::map<idx_t, idx_t> GlobalIndexToNumber;	// cell global index to number from [0, nCells-1]
	std::map<idx_t, idx_t> NumberToGlobalIndex;	// number from [0, nCells-1] to cell global index	
};

class GridStructure {
public:
	//Flag set to true if grid was properly created
	bool Initialized;

	//Grid actual cell dimensions
	int nCellDimensions;

	//Total number of nodes and cells in grid
	int nNodesGrid;
	int nCellsGrid;

	//Cells and partitioning information
	int rank;	//Rank of current node
	int np;		//Number of computational nodes
	std::vector<Cell> cells;	//All cells stored for now
	std::vector<int> cellsPart; //Partitioning information for cells

	//Grid connectivity graph via METIS style datastructures
	std::vector<idx_t> vdist;	//cells distribution over processors
	std::vector<idx_t> xadj;	//adjacency list shift for local cells
	std::vector<idx_t> adjncy;	//concateneted adjacency lists for each local cell	

	//Local cells, faces and nodes
	int nProperCellsLocal;  //Number of local non dummy cells
	int nCellsLocal;		//Number of local cells (including dummy)
	int nFacesLocal;		//Number of local faces
	//std::vector<std::shared_ptr<Cell> > localNodes;
	//std::vector<std::shared_ptr<Cell> > localFaces;
	//std::vector<std::shared_ptr<Cell> > localCells;
};

class DistributedMesh
{		
	//Underlying MPI implementation
	std::shared_ptr<ParallelManager> _MPIManager;
public:	
	//basic grid structure
	GridInfo gridInfo;

	//Grid connectivity graph
	//METIS style datastructures
	std::vector<idx_t> vdist;	//cells distribution over processors
	std::vector<idx_t> xadj;	//adjacency list shift for local cells
	std::vector<idx_t> adjncy;	//concateneted adjacency lists for each local cell	

	//local part of grid

	//Nodes
	std::vector<Node> localNodes; // all grid nodes

	//Cells	
	int nCells; //total number of cells
	int nCellsLocal; //number of local cells (without dummy)
	int nProperCells; //number of non dummy cells
	int nDummyCells;  //number of dummy cells
	std::vector<Cell> Cells; // all cells
	std::vector<Cell*> localCells; // only local cells (including dummy cells)
	std::set<int> dummyLocalCells; // dummy cells created for each boundary face element	
	std::map<int, int> cellsGlobalToLocal; // cells global index map to local index

	//Topological boundary conditions (e.g. periodic) 
	std::map<int, std::set<int>> periodicNodesIdentityList; // for each node list of nodes identifyed with it via periodic boundary
	std::map<int, int> FaceToFaceIdentity; //For each face store it's paired face through topological mapping
	std::map<int, Vector> FaceToFaceTranslation; 
	std::map<int, Matrix> FaceToFaceTransformation; 

	//Faces
	int nFaces; //number of local faces
	std::vector<Face> localFaces; // only local faces

	//Utility functions
	inline bool IsDummyCell(int globalIndex) {
		return Cells[globalIndex].IsDummy;
	};

	inline bool IsBoundaryFace(Face& face) {
		return IsDummyCell(face.FaceCell_2);
	};

	inline int GetCellPart(int globalIndex) {
		return cellsPartitioning[globalIndex];
	};

	//partitioning info
	std::vector<int> cellsPartitioning; //map from cell number index to processor rank
	std::vector<int> localCellIndexes;  //indexes of local non dummy cells

	//boundaries information (BCMarker -> Patch)	
	std::map<int, Patch> patches;
	std::map<std::string, int> patchesNames;

	//refresh boundary nodes
	std::set<int> boundaryNodes;
	std::set<int> innerNodes;
	std::set<int> fixedBoundaryNodes;
	std::set<int> movableBoundaryNodes;

	//Refresh list of boundary nodes
	void RefreshBoundaryNodes()	{		
		for (int i = nCellsLocal; i<localCells.size(); i++) {
			Cell* dummyCell = localCells[i];			
			for (int j = 0; j < dummyCell->Nodes.size(); j++) {
				int index = dummyCell->Nodes[j];
				if (boundaryNodes.find(index) != boundaryNodes.end()) continue;
				boundaryNodes.insert(index);
			};			
		};					
	};

	//Geometric properties computation
	std::vector<Face> ObtainFaces(Cell* cell);
	void ComputeGeometricProperties(Cell* cell);
	void ComputeGeometricProperties(Face* face);
	
	//Patch operation section
	void addPatch(std::string name, int bcMarker) {		
		patches[bcMarker] = Patch(name, bcMarker);
		patchesNames[name] = bcMarker;
	};

	//Build alglib kdtree for all boundary nodes
	void BuildBoundaryKDTree(Patch& patch);

	//Partition grid
	void PartitionGrid(std::shared_ptr<ParallelManager>);

	//Generate local cells given partitioning
	void GenerateLocalCells(int rank, const std::vector<int>& cellsPart);

	//Generate local faces given partitioning
	void GenerateLocalFaces(int rank);	

	//Update geometric properties of cells and faces
	void UpdateGeometricProperties();

	//Construct patches (fill in nodes and faces)
	bool ConstructAndCheckPatches();

	//Refactored grid implementation
	GridStructure gridStructure; //Public information of grid structure

	//Constructors
	Grid() {
		gridStructure.Initialized = false;
	};
	Grid(std::shared_ptr<ParallelManager> MPIManager) : _MPIManager(MPIManager) {
		gridStructure.Initialized = false;
	};

	void GenerateLocalCells();	//Generate local cells given partitioning
	void GenerateLocalFaces();	//Generate local faces given partitioning
	void GenerateLocalGrid();	//Generate local part of grid
	//void PartitionGrid();		//Partition grid
	
};

/// Implementation part
//  Build alglib kdtree for all boundary nodes
void Grid::BuildBoundaryKDTree(Patch& patch) {			
	alglib::ae_int_t nx = 3;
	alglib::ae_int_t ny = 1;
	alglib::ae_int_t normtype = 2;		
	alglib::real_2d_array a;  
	a.setlength(nx + ny, patch.nodes_idx.size());
	int i = 0;
	for (std::set<int>::iterator it = patch.nodes_idx.begin(); it != patch.nodes_idx.end(); it++) {			
		a[0][i] = localNodes[*it].P.x;
		a[1][i] = localNodes[*it].P.y;
		a[2][i] = localNodes[*it].P.z;
		a[3][i] = localNodes[*it].GlobalIndex;
		i++;
	};
	
	// Build tree	
	alglib::kdtreebuild(a, nx, ny, normtype, patch.kdt);		
};

//Generate local cells given partitioning and processor rank
void Grid::GenerateLocalCells(int rank, const std::vector<int>& cellsPart) {
	cellsPartitioning.resize(nProperCells);	

	//Add dummy cells
	for (int i = nProperCells; i<nCells; i++) {
		Cell* cell = &Cells[i];
		int neigbour = cell->NeigbourCells[0];
		int p = cellsPartitioning[neigbour];
		cellsPartitioning.push_back(p);			
	};

	//Extract local cells indexes		
	int count = 0;
	localCells.clear();
	localCellIndexes.clear();
	for (int i = 0; i<nCells; i++) {
		if ((cellsPart[i] == rank) && (!Cells[i].IsDummy)) {
			localCells.push_back(&Cells[i]);
			localCellIndexes.push_back(i);				
		};
	};		
	nCellsLocal = localCellIndexes.size(); //Without dummy cells
	for (int i = 0; i<nCells; i++) {
		if ((cellsPart[i] == rank) && (Cells[i].IsDummy)) {
			localCells.push_back(&Cells[i]);
			localCellIndexes.push_back(i);
		};
	};		

	//Assign local indexes
	int localIndex = 0;
	for (int i = 0; i<cellsPart.size(); i++) {
		if (cellsPart[i] == rank) {
			Cell* newCell = &Cells[i];			
			localCells.push_back(newCell);
			cellsGlobalToLocal[i] = localIndex++;
		};
	};

	//Create local cells
	localCells.resize(localCellIndexes.size());		
	for (int i = 0; i < localCellIndexes.size(); i++) {	
		int cellGlobalIndex = localCellIndexes[i];
		int cellLocalIndex = i;
		localCells[i] = &Cells[cellGlobalIndex];	
		localCells[i]->Faces.clear();
		cellsGlobalToLocal[cellGlobalIndex] = cellLocalIndex;
	};	
};

//Generate local faces
void Grid::GenerateLocalFaces(int rank) {
	//Create faces		
	localFaces.clear();
	int faceIndex = 0;
	std::map<std::set<int>, int> faces;	//Face nodes to face index		
	for (int i = 0; i < localCellIndexes.size(); i++) {
		//Generate all faces for the cell
		Cell* cell = localCells[i];		
		std::vector<Face> newFaces;
		if (!cell->IsDummy) {
			newFaces = ObtainFaces(cell);
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
				localFaces.push_back(face);
				localFaces[faceIndex].isExternal = true;
				//localFaces[faceIndex].FaceCell_1 = i;
				localFaces[faceIndex].FaceCell_1 = cell->GlobalIndex;
				localCells[i]->Faces.push_back(face.GlobalIndex);
				faceIndex++;
			} else {					
				int fIndex = result.first->second;
				localFaces[fIndex].isExternal = false;
				if (cell->IsDummy) {
					localFaces[fIndex].FaceCell_2 = cell->GlobalIndex;
				} else {
					//localFaces[faceIndex].FaceCell_2 = i;
					localFaces[fIndex].FaceCell_2 = cell->GlobalIndex;
				};
				localCells[i]->Faces.push_back(fIndex);
			};
		};
	};

	//Finish generating external interprocessor faces
	std::stringstream msg;
	msg.clear();
	int nExternal = 0;
	for (Face& face : localFaces) if (face.isExternal) {
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External faceID ", face.GlobalIndex);				
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External face FaceCell_1 ", face.FaceCell_1);
		int faceNodes = face.FaceNodes.size();
		//msg<<"faceNodes = "<<faceNodes<<"\n";
		Cell& cell = Cells[face.FaceCell_1];
		for (int nCellID : cell.NeigbourCells) {
			//msg<<"nCellID = "<<faceNodes<<"\n";
			Cell& nCell = Cells[nCellID];
			int nCommonNodes = 0;
			std::set<int> cellNodes;
			for (int cellNodeID : nCell.Nodes) {
				cellNodes.insert(cellNodeID);
				//Handle periodic
				if (periodicNodesIdentityList.find(cellNodeID) != periodicNodesIdentityList.end()) {
					cellNodes.insert(periodicNodesIdentityList[cellNodeID].begin(), periodicNodesIdentityList[cellNodeID].end());
				};
				//msg<<"cellNodeID = "<<cellNodeID<<"\n";
			};
			for (int nodeID : face.FaceNodes) {
				//msg<<"nodeID = "<<nodeID<<"\n";
				if (cellNodes.find(nodeID) != cellNodes.end()) nCommonNodes++;
			};
			if (nCommonNodes == faceNodes) {
				//We found second cell
				if (GetCellPart(nCellID) == rank) {
					face.isExternal = false;
				} else {
					face.isExternal = true; nExternal++;
				};
				face.FaceCell_2 = nCellID;
				break;
			};
		};
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External face FaceCell_2 ", face.FaceCell_2);
	};

	nFaces = faceIndex;	
};

//Generate local part of grid
void Grid::GenerateLocalGrid() {
	if (!gridStructure.Initialized) throw new Exception("Grid wasn't initialized properly before calling GenerateLocalGrid()");
	GenerateLocalCells();
	GenerateLocalFaces();
	UpdateGeometricProperties();
};

//Given type of element and nodes fill the properties of cell
void Grid::ComputeGeometricProperties(Cell* cell) {
	//Check cell type and compute cell volume
	bool isTypeCheckPassed = false;

	if (cell->CGNSType == HEXA_8) {
		isTypeCheckPassed = true;

		cell->CellVolume = 0;
		Vector rOrigin = localNodes[cell->Nodes[0]].P;
		for (int i = 0; i<cell->Faces.size(); i++) {
			//Get face
			Face& face = localFaces[cell->Faces[i]];
			//Consider normal direction
			double direction = 1;
			if (cell->GlobalIndex != face.FaceCell_1) direction = -1;
			//Compute normal vector
			cell->CellVolume += direction * (face.FaceCenter - rOrigin) * face.FaceSquare * face.FaceNormal;
		};
		cell->CellVolume /= 3;
	};

	//
	if (cell->CGNSType == PENTA_6) {
		isTypeCheckPassed = true;

		cell->CellVolume = 0;
		Vector rOrigin = localNodes[cell->Nodes[0]].P;
		for (int i = 0; i<cell->Faces.size(); i++) {
			//Get face
			Face& face = localFaces[cell->Faces[i]];
			//Consider normal direction
			double direction = 1;
			if (cell->GlobalIndex != face.FaceCell_1) direction = -1;
			//Compute normal vector
			cell->CellVolume += direction * (face.FaceCenter - rOrigin) * face.FaceSquare * face.FaceNormal;
		};
		cell->CellVolume /= 3;
	};

	if (cell->CGNSType == TETRA_4) {
		isTypeCheckPassed = true;

		cell->CellVolume = 0;
		Vector rOrigin = localNodes[cell->Nodes[0]].P;
		Vector a = localNodes[cell->Nodes[1]].P - rOrigin;
		Vector b = localNodes[cell->Nodes[2]].P - rOrigin;
		Vector c = localNodes[cell->Nodes[3]].P - rOrigin;		
		cell->CellVolume = abs((a & b) * c);
		cell->CellVolume /= 3;
	};

	if (cell->CGNSType == QUAD_4) {
		isTypeCheckPassed = true;

		//Compute cell volume
		cell->CellVolume = 0;
		Vector a = localNodes[cell->Nodes[1]].P - localNodes[cell->Nodes[0]].P; 
		Vector b = localNodes[cell->Nodes[2]].P - localNodes[cell->Nodes[0]].P; 
		cell->CellVolume += (a & b).mod() / 2.0;
		a = localNodes[cell->Nodes[2]].P - localNodes[cell->Nodes[0]].P; 
		b = localNodes[cell->Nodes[3]].P - localNodes[cell->Nodes[0]].P; 
		cell->CellVolume += (a & b).mod() / 2.0;
	};

	if (cell->CGNSType == TRI_3) {
		isTypeCheckPassed = true;

		//Compute cell volume
		cell->CellVolume = 0;
		Vector a = localNodes[cell->Nodes[1]].P - localNodes[cell->Nodes[0]].P; 
		Vector b = localNodes[cell->Nodes[2]].P - localNodes[cell->Nodes[0]].P; 
		cell->CellVolume += (a & b).mod() / 2.0;		
	};

	if (cell->CGNSType == BAR_2) {
		isTypeCheckPassed = true;

		//Compute cell volume
		Vector a = localNodes[cell->Nodes[1]].P - localNodes[cell->Nodes[0]].P; 
		cell->CellVolume = a.mod();		
	};

	if (DEBUG_GRID) {
		//Check all volumes for non-negativity
		assert(cell->CellVolume >= 0.0);
	};
	if (!isTypeCheckPassed) throw Exception("Cell element type is not supported");

	//Compute cell center (independent from type)
	cell->CellCenter = Vector(0,0,0);
	for (int i = 0; i<cell->Nodes.size(); i++) cell->CellCenter += localNodes[cell->Nodes[i]].P;
	cell->CellCenter /= cell->Nodes.size();
};

void Grid::ComputeGeometricProperties(Face* face) {
	//Check cell type and compute face square and normal
	bool isTypeCheckPassed = false;

	if (face->CGNSType == QUAD_4) {
		isTypeCheckPassed = true;
		
		//Compute surface area and normal
		face->FaceNormal = Vector(0,0,0);
		Vector a = localNodes[face->FaceNodes[1]].P - localNodes[face->FaceNodes[0]].P; 
		Vector b = localNodes[face->FaceNodes[2]].P - localNodes[face->FaceNodes[0]].P; 
		face->FaceNormal += (a & b) / 2.0;
		a = localNodes[face->FaceNodes[2]].P - localNodes[face->FaceNodes[0]].P; 
		b = localNodes[face->FaceNodes[3]].P - localNodes[face->FaceNodes[0]].P; 
		face->FaceNormal += (a & b) / 2.0;

		//Normalize	
		face->FaceSquare = face->FaceNormal.mod();
		face->FaceNormal /= face->FaceSquare;		
	};

	if (face->CGNSType == TRI_3) {
		isTypeCheckPassed = true;

		//Compute surface area and normal
		face->FaceNormal = Vector(0,0,0);
		Vector a = localNodes[face->FaceNodes[1]].P - localNodes[face->FaceNodes[0]].P; 
		Vector b = localNodes[face->FaceNodes[2]].P - localNodes[face->FaceNodes[0]].P; 
		face->FaceNormal += (a & b) / 2.0;		

		//Normalize	
		face->FaceSquare = face->FaceNormal.mod();
		face->FaceNormal /= face->FaceSquare;		
	};

	if (face->CGNSType == BAR_2) {
		isTypeCheckPassed = true;

		//Compute surface area
		face->FaceSquare = (localNodes[face->FaceNodes[1]].P - localNodes[face->FaceNodes[0]].P).mod();		

		//Compute normal		
		face->FaceNormal.x = (localNodes[face->FaceNodes[1]].P - localNodes[face->FaceNodes[0]].P).y;
		face->FaceNormal.y = -(localNodes[face->FaceNodes[1]].P - localNodes[face->FaceNodes[0]].P).x;
		face->FaceNormal = face->FaceNormal / face->FaceNormal.mod();
	};

	if (face->CGNSType == NODE) {		
		isTypeCheckPassed = true;

		//Compute surface area
		face->FaceSquare = 1.0;

		//Compute normal		
		face->FaceNormal = Vector(1,0,0);		
		face->FaceNormal = face->FaceNormal / face->FaceNormal.mod();		
	};

	if (DEBUG_GRID) {
		//Check all face squares for non-negativity
		assert(face->FaceSquare >= 0.0);
		//And normals for being unit length
		assert(abs(face->FaceNormal.mod() - 1.0) <= 1e-12);
	};
	if (!isTypeCheckPassed) throw Exception("Face element type is not supported");
	
	//Compute face center (independent from type)
	face->FaceCenter = Vector(0,0,0);
	for (int i = 0; i<face->FaceNodes.size(); i++) face->FaceCenter += localNodes[face->FaceNodes[i]].P;
	face->FaceCenter /= face->FaceNodes.size();	
};

//Given type of cell element and nodes obtain face without geometric properties
std::vector<Face> Grid::ObtainFaces(Cell* cell) {
	std::vector<Face> res;

	//Check cell type
	bool isTypeCheckPassed = false;

	if (cell->CGNSType == HEXA_8) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;

		//Add faces one by one		
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[4]);
		nodes.push_back(cell->Nodes[7]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[5]);
		nodes.push_back(cell->Nodes[4]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[4]);
		nodes.push_back(cell->Nodes[5]);
		nodes.push_back(cell->Nodes[6]);
		nodes.push_back(cell->Nodes[7]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();
		
		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[5]);
		nodes.push_back(cell->Nodes[6]);
		nodes.push_back(cell->Nodes[2]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[6]);
		nodes.push_back(cell->Nodes[7]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();
	};

	if (cell->CGNSType == PENTA_6) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one		
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[4]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[5]);
		nodes.push_back(cell->Nodes[4]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[3]);
		nodes.push_back(cell->Nodes[5]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(TRI_3, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[3]);
		nodes.push_back(cell->Nodes[4]);
		nodes.push_back(cell->Nodes[5]);
		res.push_back(Face(TRI_3, nodes));
		nodes.clear();
	};

	if (cell->CGNSType == TETRA_4) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[2]);
		res.push_back(Face(TRI_3, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(TRI_3, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(TRI_3, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(TRI_3, nodes));
		nodes.clear();
	};

	if (cell->CGNSType == QUAD_4) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[2]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[3]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[3]);
		nodes.push_back(cell->Nodes[0]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();
	}

	if (cell->CGNSType == TRI_3) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one
		nodes.push_back(cell->Nodes[0]);
		nodes.push_back(cell->Nodes[1]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[1]);
		nodes.push_back(cell->Nodes[2]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[2]);
		nodes.push_back(cell->Nodes[0]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();		
	}

	if (cell->CGNSType == BAR_2) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one
		nodes.push_back(cell->Nodes[0]);		
		res.push_back(Face(NODE, nodes));
		nodes.clear();

		nodes.push_back(cell->Nodes[1]);		
		res.push_back(Face(NODE, nodes));
		nodes.clear();
	}

	if (!isTypeCheckPassed) throw Exception("Cell element type is not supported");
	return res;
};

void Grid::UpdateGeometricProperties() {
	//Generate face geometric properties		
	for (Face& face : localFaces) {
		int index = face.GlobalIndex;						
		ComputeGeometricProperties(&localFaces[index]);
	};

	//Update cells geometric properties
	for (int i = 0; i < nCellsLocal; i++) {
		ComputeGeometricProperties(localCells[i]);
		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, localCells[i]->getInfo());
	};	

	//Orient face normals			
	for (Face& face : localFaces) {
		int index = face.GlobalIndex;									

		//If boundary face make sure dummy if FaceCell_2
		if (face.FaceCell_1 >= nProperCells) {
			//Swap
			int tmp = face.FaceCell_1;
			face.FaceCell_1 = face.FaceCell_2;
			face.FaceCell_2 = tmp;
		};

		//Orient normal
		Vector cellCenter = Cells[localFaces[index].FaceCell_1].CellCenter;
		Vector faceCenter = localFaces[index].FaceCenter;
		if (((cellCenter - faceCenter) * localFaces[index].FaceNormal) > 0) {
			localFaces[index].FaceNormal *= -1;
		};

		//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, face.getInfo());
	};


	//Compute dummy cell geometric properties
	for (int i = nCellsLocal; i < localCells.size(); i++) {			
		int neighbour = localCells[i]->NeigbourCells[0];
		Cell& cell = Cells[neighbour];
		localCells[i]->CellVolume = cell.CellVolume;

		//Reflect cell center over boundary face plane
		Face& face = localFaces[ localCells[i]->Faces[0]];
		Vector dR = ((cell.CellCenter - face.FaceCenter) * face.FaceNormal) * face.FaceNormal / face.FaceNormal.mod();
		Vector dummyCenter = cell.CellCenter - 2 * dR;
		localCells[i]->CellCenter = dummyCenter;
		//s_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, localCells[i]->getInfo());
	};
};

/* Check that all boundary markers are set to some patchs 
And fill in information about nodes and faces */
bool Grid::ConstructAndCheckPatches() {
	bool checkResult = true;
	std::vector<std::string> reasons;	
	for (int i = nCellsLocal; i<localCells.size(); i++) {
		Cell* dummyCell = localCells[i];		
		int bcMarker = dummyCell->BCMarker;
		int globalIndex = dummyCell->GlobalIndex;
		if (patches.find(bcMarker) == patches.end()) {
			std::string msg;			
			std::cout<<"Not found patch for boundary marker ";
			std::cout<<bcMarker;
			std::cout<<" for faceID = ";
			std::cout<<globalIndex<<"\n";
			//reasons.push_back("Not found patch for boundary marker " + bcMarker + " for faceID = " + globalIndex);
			checkResult = false;
		} else {
			//Add nodes and faces
			patches[dummyCell->BCMarker].addFace(dummyCell->GlobalIndex);
			for (int j = 0; j<dummyCell->Nodes.size(); j++) patches[dummyCell->BCMarker].addNode(dummyCell->Nodes[j]);
		};
	};
	return checkResult;
};

//New grid implementation

//Partition grid
void Grid::PartitionGrid(std::shared_ptr<ParallelManager> pHelper) {
	_MPIManager = std::shared_ptr<ParallelManager>(pHelper);
	gridStructure.np = _MPIManager->np();
	gridStructure.rank = _MPIManager->rank();

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

	//Number of subdomains desired
	idx_t nparts = gridStructure.np; 

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
	std::vector<real_t> ubvec(ncon);	// imbalance tolerance for each vertex weight,
	ubvec[0] = 1.02;					// recomended default value 

	//Algoritm options for displaing information
	idx_t options[3];
	options[0] = 0; //No information and default values

	idx_t edgecut;	//Upon successful completion, the number of edges that are cut by the partitioning is written to this parameter.

	/* This is an array of size equal to the number of locally-stored vertices. Upon successful completion the
	partition vector of the locally-stored vertices is written to this array. (See discussion in Section 4.2.4). */		 
	std::vector<idx_t> part(nCellsLocal);

	//Call partitioning function		
	MPI_Comm _comm = pHelper->comm();
	pHelper->Barrier();
	if (gridStructure.np < nProperCells) { //Make sure that partitioning is needed
		int result = ParMETIS_V3_PartKway(&vdist[0], &xadj[0], adjncy._Myfirst, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, &tpwgts[0], &ubvec[0], options, &edgecut, &part[0], &_comm);			
		if (result != METIS_OK) {
			throw new Exception("ParMETIS_V3_PartKway failed.");
		};							
	};

	//Gather partitioning on every processor
	std::vector<int> recvcounts(gridStructure.np);
	for (int i = 0; i<gridStructure.np; i++) {
		recvcounts[i] = vdist[i+1] - vdist[i];
	};
		
	/*msg.clear();
	msg<<part.size()<<"\n";
	_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	

	msg.clear();
	msg<<"recvcounts[] = \n";
	for (int i = 0; i<_nProcessors; i++) msg<<recvcounts[i]<<" ";
	msg<<"\n";

	_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, msg.str());	*/
	_MPIManager->Allgatherv<int, MPI_INT>( part, recvcounts, cellsPartitioning);

	GenerateLocalCells(pHelper->rank(), cellsPartitioning);

	//Synchronize		
	pHelper->Barrier();
};

//Generate local cells given partitioning and processor rank
void Grid::GenerateLocalCells() {
	//Add dummy cells
	//for (int i = nProperCells; i<nCells; i++) {
	//	Cell* cell = &Cells[i];
	//	int neigbour = cell->NeigbourCells[0];
	//	int p = gridStructure.cellsPart[neigbour];
	//	gridStructure.cellsPart.push_back(p);			
	//};

	////Extract local cells indexes		
	//int count = 0;
	//localCells.clear();
	//localCellIndexes.clear();
	//for (int i = 0; i<nCells; i++) {
	//	if ((gridStructure.cellsPart[i] == gridStructure.rank) && (!gridStructure.cells[i].IsDummy)) {
	//		localCells.push_back(&Cells[i]);
	//		localCellIndexes.push_back(i);				
	//	};
	//};		
	//nCellsLocal = localCellIndexes.size(); //Without dummy cells
	//for (int i = 0; i<nCells; i++) {
	//	if ((gridStructure.cellsPart[i] == gridStructure.rank) && (gridStructure.cells[i].IsDummy)) {
	//		localCells.push_back(&Cells[i]);
	//		localCellIndexes.push_back(i);
	//	};
	//};		

	////Assign local indexes
	//int localIndex = 0;
	//for (int i = 0; i<cellsPart.size(); i++) {
	//	if (cellsPart[i] == rank) {
	//		Cell* newCell = &Cells[i];			
	//		localCells.push_back(newCell);
	//		cellsGlobalToLocal[i] = localIndex++;
	//	};
	//};

	////Create local cells
	//localCells.resize(localCellIndexes.size());		
	//for (int i = 0; i < localCellIndexes.size(); i++) {	
	//	int cellGlobalIndex = localCellIndexes[i];
	//	int cellLocalIndex = i;
	//	localCells[i] = &Cells[cellGlobalIndex];	
	//	localCells[i]->Faces.clear();
	//	cellsGlobalToLocal[cellGlobalIndex] = cellLocalIndex;
	//};	
};

//Generate local faces
void Grid::GenerateLocalFaces() {
	//Create faces		
	//localFaces.clear();
	//int faceIndex = 0;
	//std::map<std::set<int>, int> faces;	//Face nodes to face index		
	//for (int i = 0; i < localCellIndexes.size(); i++) {
	//	//Generate all faces for the cell
	//	Cell* cell = localCells[i];		
	//	std::vector<Face> newFaces;
	//	if (!cell->IsDummy) {
	//		newFaces = ObtainFaces(cell);
	//	} else {
	//		Face boundaryFace;
	//		boundaryFace.CGNSType = cell->CGNSType;
	//		boundaryFace.FaceNodes = cell->Nodes;
	//		boundaryFace.BCMarker = cell->BCMarker;
	//		newFaces.clear();
	//		newFaces.push_back(boundaryFace);
	//	};
	//	for (int j = 0; j<newFaces.size(); j++) {
	//		Face& face = newFaces[j];				
	//		std::set<int> nodes(face.FaceNodes.begin(), face.FaceNodes.end());
	//		std::pair<std::set<int>, int> newFaceInfo(nodes, 0);
	//		std::pair<std::map<std::set<int>, int>::iterator, bool> result = faces.insert( newFaceInfo );				
	//		//Add face and save connectivity info
	//		if ( result.second ) {	
	//			face.GlobalIndex = faceIndex;
	//			result.first->second = face.GlobalIndex;
	//			localFaces.push_back(face);
	//			localFaces[faceIndex].isExternal = true;
	//			//localFaces[faceIndex].FaceCell_1 = i;
	//			localFaces[faceIndex].FaceCell_1 = cell->GlobalIndex;
	//			localCells[i]->Faces.push_back(face.GlobalIndex);
	//			faceIndex++;
	//		} else {					
	//			int fIndex = result.first->second;
	//			localFaces[fIndex].isExternal = false;
	//			if (cell->IsDummy) {
	//				localFaces[fIndex].FaceCell_2 = cell->GlobalIndex;
	//			} else {
	//				//localFaces[faceIndex].FaceCell_2 = i;
	//				localFaces[fIndex].FaceCell_2 = cell->GlobalIndex;
	//			};
	//			localCells[i]->Faces.push_back(fIndex);
	//		};
	//	};
	//};

	////Finish generating external interprocessor faces
	//std::stringstream msg;
	//msg.clear();
	//int nExternal = 0;
	//for (Face& face : localFaces) if (face.isExternal) {
	//	//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External faceID ", face.GlobalIndex);				
	//	//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External face FaceCell_1 ", face.FaceCell_1);
	//	int faceNodes = face.FaceNodes.size();
	//	//msg<<"faceNodes = "<<faceNodes<<"\n";
	//	Cell& cell = Cells[face.FaceCell_1];
	//	for (int nCellID : cell.NeigbourCells) {
	//		//msg<<"nCellID = "<<faceNodes<<"\n";
	//		Cell& nCell = Cells[nCellID];
	//		int nCommonNodes = 0;
	//		std::set<int> cellNodes;
	//		for (int cellNodeID : nCell.Nodes) {
	//			cellNodes.insert(cellNodeID);
	//			//Handle periodic
	//			if (periodicNodesIdentityList.find(cellNodeID) != periodicNodesIdentityList.end()) {
	//				cellNodes.insert(periodicNodesIdentityList[cellNodeID].begin(), periodicNodesIdentityList[cellNodeID].end());
	//			};
	//			//msg<<"cellNodeID = "<<cellNodeID<<"\n";
	//		};
	//		for (int nodeID : face.FaceNodes) {
	//			//msg<<"nodeID = "<<nodeID<<"\n";
	//			if (cellNodes.find(nodeID) != cellNodes.end()) nCommonNodes++;
	//		};
	//		if (nCommonNodes == faceNodes) {
	//			//We found second cell
	//			if (GetCellPart(nCellID) == rank) {
	//				face.isExternal = false;
	//			} else {
	//				face.isExternal = true; nExternal++;
	//			};
	//			face.FaceCell_2 = nCellID;
	//			break;
	//		};
	//	};
	//	//_logger.WriteMessage(LoggerMessageLevel::LOCAL, LoggerMessageType::INFORMATION, "External face FaceCell_2 ", face.FaceCell_2);
	//};

	//nFaces = faceIndex;	
};

#endif


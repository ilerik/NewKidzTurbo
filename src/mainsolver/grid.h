#ifndef TURBO_GRID
#define TURBO_GRID

#include "basetypes.h"
#include "alglibmisc.h"
#include "cgnslib.h"
#include "patch.h"
#include "parmetis.h"
#include <map>
#include <set>
//#include <unordered_set>
//#include <unordered_map>
#include <algorithm>
#include <cassert>

#define DEBUG_GRID true

//Specify supported element types


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
  ElementType_t CGNSType;
  int GlobalIndex;
	
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
	};

	Face(ElementType_t type, std::vector<int> nodes) {
		CGNSType = type;
		FaceNodes = nodes;
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

class Grid
{		
public:	
	//basic grid properties	
	GridInfo gridInfo;

	//grid connectivity info
	//METIS style datastructures
	std::vector<idx_t> vdist;  //cells distribution over processors
	std::vector<idx_t> xadj;   //adjacency list shift for local cells
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
	std::map<int, std::set<int>> periodicNodesIdentityList; // for each node list of nodes identifyed with it via periodic boundary

	//Faces
	int nFaces; //number of local faces
	std::vector<Face> localFaces; // only local faces

	//Utility functions
	inline bool IsDummyCell(int globalIndex) {
		return Cells[globalIndex].IsDummy;
	};

	inline int GetCellPart(int globalIndex) {
		return cellsPartitioning[globalIndex];
	};

	inline bool IsBoundaryFace(Face& face) {
		return IsDummyCell(face.FaceCell_2);
	};

	//partitioning info
	std::vector<int> cellsPartitioning; // map from cell number index to processor rank
	std::vector<int> localCellIndexes;  // indexes of local non dummy cells

	//boundaries information (BCMarker -> Patch)	
	std::map<int, Patch> patches;
	std::map<std::string, int> patchesNames;

	//refresh boundary nodes
	std::set<int> boundaryNodes;
	std::set<int> innerNodes;
	std::set<int> fixedBoundaryNodes;
	std::set<int> movableBoundaryNodes;

	//Refresh list of boundary nodes
	void RefreshBoundaryNodes()
	{		
		for (int i = nCellsLocal; i<localCells.size(); i++) {
			Cell* dummyCell = localCells[i];			
			for (int j = 0; j < dummyCell->Nodes.size(); j++) {
				int index = dummyCell->Nodes[j];
				if (boundaryNodes.find(index) != boundaryNodes.end()) continue;
				boundaryNodes.insert(index);
			};			
		};					
	};
		
	Grid(void) {	
		gridInfo.GridDimensions = 2;		//Default value
	};

	~Grid(void){
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

	//Generate local cells given partitioning
	void GenerateLocalCells(int rank, std::vector<int>& cellsPart);

	//Update geometric properties of cells and faces
	void UpdateGeometricProperties();

	//Generate local faces given partitioning
	void GenerateLocalFaces();	

	//Construct patches (fill in nodes and faces)
	bool ConstructAndCheckPatches();

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
void Grid::GenerateLocalCells(int rank, std::vector<int>& cellsPart) {
	localCells.clear();

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
void Grid::GenerateLocalFaces() {

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

//Check that all boundary markers are set to some patchs
//And fill in information about nodes and faces
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

#endif


#ifndef TURBO_GRID
#define TURBO_GRID

#include "basetypes.h"
#include "alglibmisc.h"
#include "cgnslib.h"
#include "patch.h"
#include <map>
#include <set>
#include <unordered_set>

//Grid types
struct Node {
	int GlobalIndex; //Global index of node
	Vector P;	//Position	
};

struct Face {
	ElementType_t CGNSType;
	int GlobalIndex;
	int FaceCell_1;		//Index of cell 1
	int FaceCell_2;		//Index of cell 2
	Vector FaceCenter;	//Center of face
	Vector FaceNormal;	//Surface normal * square
	double FaceSquare;	//Face square
	std::vector<int> FaceNodes;	//Indexes of all face nodes
	int isExternal;		//If its an external face (FaceCell_2 == 0)
	int BCMarker;		//Marker of boundary condition type	
	std::map<std::string,  void*> Data; //Arbitrary data	

	//Constructors
	Face() {
	};

	Face(ElementType_t type, std::vector<int> nodes) {
		CGNSType = type;
		FaceNodes = nodes;
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

class Cell {
public:
	ElementType_t CGNSType; //CGNS cell type
	int GlobalIndex;		
	std::vector<int> Faces; //Indexes of all cell faces
	double CellVolume;		// Volume of cell
	Vector CellCenter;		// Center of cell
	double CellHSize;		// Size of cell	
	std::vector<int> Nodes; //Cell  node indexes
	std::map<std::string,  void*> Data; //Arbitrary data	
};


//Parallel data node manager
//Class that handles distributed collection of arbitrary data nodes
//DataNode must posess public property int GlobalIndex
template<class DataNode> 
class DistributedEntityManager {
	std::set<int> node_idx;
	std::map<int, DataNode> localNodes;

public:
	DataNode& operator[](int index) {
		//Check if index is ok
		if (node_idx.find(index) == node_idx.cend()) throw Exception("Invalid node global index");
		//Check if its local first
		std::map<int, DataNode>::iterator it = localNodes.find(index);
		if (it != localNodes.end()) return it->second;	//Its local and it exists		
	};

	std::vector<DataNode*> getLocalNodes() {
		std::vector<DataNode*> nodes;
		for (std::map<int, DataNode>::iterator it = localNodes.begin(); it!=localNodes.end(); it++) {
			nodes.push_back(&(it->second));
		};
		return nodes;
	};

	std::map<int, DataNode>& getLocalNodesWithIndex() {		
		return localNodes;
	};

	int size() {
		return node_idx.size();
	};	

	std::set<int> getAllIndexes() {
		return node_idx;
	};

	void add(DataNode n) {
		node_idx.insert(n.GlobalIndex);
		localNodes[n.GlobalIndex] = n;
	};

	void removeAt(int index) {
		if (node_idx.find(index) != node_idx.cend()) {
			node_idx.erase(index);
			localNodes.erase(index);
		};
	};

	bool exist(int index) {
		return node_idx.find(index) != node_idx.cend();
	};
};


//general info
class GridInfo {
public:
	int nCoords;
	int CellDimensions;
	int GridDimensions;	
};

class Grid
{		
public:	
	//basic grid properties
	int GridID;
	GridInfo gridInfo;

	//boundaries information (BCMarker -> Patch)	
	std::map<int, Patch> patches;
	std::map<std::string, int> patchesNames;

	//refresh boundary nodes
	std::set<int> boundaryNodes;
	std::set<int> innerNodes;
	std::set<int> fixedBoundaryNodes;
	std::set<int> movableBoundaryNodes;

	void RefreshBoundaryNodes()
	{
		std::vector<Face*> f = faces.getLocalNodes();
		for (int i = 0; i<f.size(); i++) {
			if (f[i]->isExternal) {
				for (int j = 0; j < f[i]->FaceNodes.size(); j++) {
					int index = f[i]->FaceNodes[j];
					if (boundaryNodes.find(index) != boundaryNodes.end()) continue;
					boundaryNodes.insert(index);
				};
			};			
		};					
	};
		

	Grid(void) {	
		gridInfo.GridDimensions = 2;		//Default value
	};

	~Grid(void){
	};

	DistributedEntityManager<Node> nodes;
	DistributedEntityManager<Cell> cells;
	DistributedEntityManager<Face> faces;

	//Geometric properties computation
	std::vector<Face> Grid::ObtainFaces(Cell& cell);
	void ComputeGeometricProperties(Cell& cell);
	void ComputeGeometricProperties(Face& face);
	
	//Patch operation section
	void addPatch(std::string name, int bcMarker) {		
		patches[bcMarker] = Patch(name, bcMarker);
		patchesNames[name] = bcMarker;
	};

	//Build alglib kdtree for all boundary nodes
	void BuildBoundaryKDTree(Patch& patch);

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
		a[0][i] = nodes[*it].P.x;
		a[1][i] = nodes[*it].P.y;
		a[2][i] = nodes[*it].P.z;
		a[3][i] = nodes[*it].GlobalIndex;
		i++;
	};
	
	// Build tree	
	alglib::kdtreebuild(a, nx, ny, normtype, patch.kdt);		
};

//Given type of element and nodes fill the properties of cell
void Grid::ComputeGeometricProperties(Cell& cell) {
	//Check cell type and compute cell volume
	bool isTypeCheckPassed = false;

	if (cell.CGNSType == HEXA_8) {
		isTypeCheckPassed = true;

		cell.CellVolume = 0;
		Vector rOrigin = nodes[cell.Nodes[0]].P;
		for (int i = 0; i<cell.Faces.size(); i++) {
			//Get face
			Face& face = faces[cell.Faces[i]];
			//Consider normal direction
			double direction = 1;
			if (cell.GlobalIndex != face.FaceCell_1) direction = -1;
			//Compute normal vector
			cell.CellVolume += direction * (face.FaceCenter - rOrigin) * face.FaceSquare * face.FaceNormal;
		};
		cell.CellVolume /= 3;
	};

	if (cell.CGNSType == QUAD_4) {
		isTypeCheckPassed = true;

		//Compute cell volume
		cell.CellVolume = 0;
		Vector a = nodes[cell.Nodes[1]].P - nodes[cell.Nodes[0]].P; 
		Vector b = nodes[cell.Nodes[2]].P - nodes[cell.Nodes[0]].P; 
		cell.CellVolume += (a & b).mod() / 2.0;
		a = nodes[cell.Nodes[2]].P - nodes[cell.Nodes[0]].P; 
		b = nodes[cell.Nodes[3]].P - nodes[cell.Nodes[0]].P; 
		cell.CellVolume += (a & b).mod() / 2.0;
	};

	if (cell.CGNSType == BAR_2) {
		isTypeCheckPassed = true;

		//Compute cell volume
		Vector a = nodes[cell.Nodes[1]].P - nodes[cell.Nodes[0]].P; 
		cell.CellVolume = a.mod();		
	};

	if (!isTypeCheckPassed) throw Exception("Cell element type is not supported");

	//Compute cell center (independent from type)
	cell.CellCenter = Vector(0,0,0);
	for (int i = 0; i<cell.Nodes.size(); i++) cell.CellCenter += nodes[cell.Nodes[i]].P;
	cell.CellCenter /= cell.Nodes.size();
};

void Grid::ComputeGeometricProperties(Face& face) {
	//Check cell type and compute face square and normal
	bool isTypeCheckPassed = false;

	if (face.CGNSType == QUAD_4) {
		isTypeCheckPassed = true;
		
		//Compute surface area and normal
		face.FaceNormal = Vector(0,0,0);
		Vector a = nodes[face.FaceNodes[1]].P - nodes[face.FaceNodes[0]].P; 
		Vector b = nodes[face.FaceNodes[2]].P - nodes[face.FaceNodes[0]].P; 
		face.FaceNormal += (a & b) / 2.0;
		a = nodes[face.FaceNodes[2]].P - nodes[face.FaceNodes[0]].P; 
		b = nodes[face.FaceNodes[3]].P - nodes[face.FaceNodes[0]].P; 
		face.FaceNormal += (a & b) / 2.0;

		//Normalize	
		face.FaceSquare = face.FaceNormal.mod();
		face.FaceNormal /= face.FaceSquare;		
	};

	if (face.CGNSType == BAR_2) {
		isTypeCheckPassed = true;

		//Compute surface area
		face.FaceSquare = (nodes[face.FaceNodes[1]].P - nodes[face.FaceNodes[0]].P).mod();		

		//Compute normal		
		face.FaceNormal.x = (nodes[face.FaceNodes[1]].P - nodes[face.FaceNodes[0]].P).y;
		face.FaceNormal.y = -(nodes[face.FaceNodes[1]].P - nodes[face.FaceNodes[0]].P).x;
		face.FaceNormal = face.FaceNormal / face.FaceNormal.mod();
	};

	if (face.CGNSType == NODE) {		
		isTypeCheckPassed = true;

		//Compute surface area
		face.FaceSquare = 1.0;

		//Compute normal		
		face.FaceNormal = (cells[face.FaceCell_1].CellCenter - face.FaceCenter);		
		face.FaceNormal = face.FaceNormal / face.FaceNormal.mod();		
	};

	if (!isTypeCheckPassed) throw Exception("Face element type is not supported");
	
	//Compute face center (independent from type)
	face.FaceCenter = Vector(0,0,0);
	for (int i = 0; i<face.FaceNodes.size(); i++) face.FaceCenter += nodes[face.FaceNodes[i]].P;
	face.FaceCenter /= face.FaceNodes.size();	
};

//Given type of cell element and nodes obtain face without geometric properties
std::vector<Face> Grid::ObtainFaces(Cell& cell) {
	std::vector<Face> res;

	//Check cell type
	bool isTypeCheckPassed = false;

	if (cell.CGNSType == HEXA_8) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;

		//Add faces one by one		
		nodes.push_back(cell.Nodes[0]);
		nodes.push_back(cell.Nodes[1]);
		nodes.push_back(cell.Nodes[2]);
		nodes.push_back(cell.Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[0]);
		nodes.push_back(cell.Nodes[4]);
		nodes.push_back(cell.Nodes[7]);
		nodes.push_back(cell.Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[0]);
		nodes.push_back(cell.Nodes[1]);
		nodes.push_back(cell.Nodes[5]);
		nodes.push_back(cell.Nodes[4]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[4]);
		nodes.push_back(cell.Nodes[5]);
		nodes.push_back(cell.Nodes[6]);
		nodes.push_back(cell.Nodes[7]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();
		
		nodes.push_back(cell.Nodes[1]);
		nodes.push_back(cell.Nodes[5]);
		nodes.push_back(cell.Nodes[6]);
		nodes.push_back(cell.Nodes[2]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[2]);
		nodes.push_back(cell.Nodes[6]);
		nodes.push_back(cell.Nodes[7]);
		nodes.push_back(cell.Nodes[3]);
		res.push_back(Face(QUAD_4, nodes));
		nodes.clear();
	};

	if (cell.CGNSType == QUAD_4) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one
		nodes.push_back(cell.Nodes[0]);
		nodes.push_back(cell.Nodes[1]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[1]);
		nodes.push_back(cell.Nodes[2]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[2]);
		nodes.push_back(cell.Nodes[3]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[3]);
		nodes.push_back(cell.Nodes[0]);
		res.push_back(Face(BAR_2, nodes));
		nodes.clear();
	}

	if (cell.CGNSType == BAR_2) {
		isTypeCheckPassed = true;
		std::vector<int> nodes;
		//Add faces one by one
		nodes.push_back(cell.Nodes[0]);		
		res.push_back(Face(NODE, nodes));
		nodes.clear();

		nodes.push_back(cell.Nodes[1]);		
		res.push_back(Face(NODE, nodes));
		nodes.clear();
	}

	if (!isTypeCheckPassed) throw Exception("Cell element type is not supported");
	return res;
};


//Check that all boundary markers are set to some patchs
//And fill in information about nodes and faces
bool Grid::ConstructAndCheckPatches() {
	bool checkResult = true;
	std::vector<std::string> reasons;
	std::vector<Face*> fcs = faces.getLocalNodes();
	for (int i = 0; i<fcs.size(); i++) {
		if (!fcs[i]->isExternal) continue;
		int bcMarker = fcs[i]->BCMarker;
		int globalIndex = fcs[i]->GlobalIndex;
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
			patches[fcs[i]->BCMarker].addFace(fcs[i]->GlobalIndex);
			for (int j = 0; j<fcs[i]->FaceNodes.size(); j++) patches[fcs[i]->BCMarker].addNode(fcs[i]->FaceNodes[j]);
		};
	};
	return checkResult;
};

#endif


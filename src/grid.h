#ifndef TURBO_GRID
#define TURBO_GRID

#include "basetypes.h"
#include "patch.h"
#include "alglibmisc.h"
#include "cgnslib.h"
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
	Vector FaceNormal;	//unit vector
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
	bool CarryOutNodeIdxControl;	//to check all nodes index exist in operator []

	DistributedEntityManager() {
		CarryOutNodeIdxControl = true;
	};

	DataNode& operator[](int index) {
		//Check if index is ok
		if(CarryOutNodeIdxControl == true)
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

	//fast working operator [] for DistributedEntityManager
	void SetFastGridOperations() {
		nodes.CarryOutNodeIdxControl = false;
		cells.CarryOutNodeIdxControl = false;
		faces.CarryOutNodeIdxControl = false;
	};

	//Geometric properties computation
	std::vector<Face> Grid::ObtainFaces(Cell& cell);
	void ComputeGeometricProperties(Cell& cell);
	void ComputeGeometricProperties(Face& face);

	//for axisymmetric 2.5D task with X symmetry axys for quadrolatiral grid
	void ComputeGeometryToAxisymmetric(double alpha, int BC_Marker);
	
	//Patch operation section
	void addPatch(std::string name, int bcMarker) {		
		patches[bcMarker] = Patch(name, bcMarker);
		patchesNames[name] = bcMarker;
	};

	//Build alglib kdtree for all boundary nodes
	void BuildBoundaryKDTree(Patch& patch);

	//Construct patches (fill in nodes and faces)
	bool ConstructAndCheckPatches();

	//write grid to TecPlot file
	void ToTecPlotGrid(char *fname);
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

	if (cell.CGNSType == TRI_3) {
		isTypeCheckPassed = true;

		//Compute cell volume
		cell.CellVolume = 0;
		Vector a = nodes[cell.Nodes[1]].P - nodes[cell.Nodes[0]].P; 
		Vector b = nodes[cell.Nodes[2]].P - nodes[cell.Nodes[0]].P; 
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

	//Compute face center (independent from type)
	face.FaceCenter = Vector(0,0,0);
	for (int i = 0; i<face.FaceNodes.size(); i++) face.FaceCenter += nodes[face.FaceNodes[i]].P;
	face.FaceCenter /= face.FaceNodes.size();	

	//compute normals
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
		face.FaceNormal = (face.FaceCenter - cells[face.FaceCell_1].CellCenter);		
		face.FaceNormal = face.FaceNormal / face.FaceNormal.mod();		
	};

	if (!isTypeCheckPassed) throw Exception("Face element type is not supported");
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

	if (cell.CGNSType == TRI_3) {
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

void Grid::ComputeGeometryToAxisymmetric(double alpha, int Periodic_BC_Marker)
{
	double ZERO = 1e-10;
	std::vector<Cell*> cells_array = cells.getLocalNodes();
	//compute max global index of faces
	int MaxFaceIndex = 0;
	std::vector<Face*> faces_array = faces.getLocalNodes();	//compute max global index of faces
	for(int i=0; i<faces.size(); i++)
	{
		int GI = faces_array[i]->GlobalIndex;
		if(MaxFaceIndex<GI) MaxFaceIndex = GI;
	};

	for(int i=0; i<cells.size(); i++)
	{
		Cell c = *cells_array[i];
		double R = c.CellCenter.y;
		int top, bottom, left, right;	//Global index of appropriate faces;
		int n_horiz = 0;		//horizontal faces number
		int n_vert = 0;			//vertical one
		//determine relative faces positions
		for(int j=0; j<4; j++)	//just four faces
		{
			Face f = faces[c.Faces[j]];
			if(abs(f.FaceNormal*Vector(1,0,0)) < ZERO)	//horizontal faces
			{
				if((f.FaceCenter.y - R)>0) top = f.GlobalIndex;
				else bottom = f.GlobalIndex;
				n_horiz++;
			}else  //vertical ones
			{
				if((f.FaceCenter.x - c.CellCenter.x)>0) right = f.GlobalIndex;
				else left = f.GlobalIndex;
				n_vert++;
			};
		};
		//check correctness
		if((n_vert!=2)||(n_horiz!=2))
		{
			std::cout << "can't recompute geometry properties for axisymmtric grid" << '\n';
			std::getchar();
			return;
		};

		double a = faces[top].FaceCenter.y - faces[bottom].FaceCenter.y;	//height of 2D cell
		double h = faces[right].FaceCenter.x - faces[left].FaceCenter.x;	//lenght of 2D cell

		//compute areas of all faces
		faces[top].FaceSquare = alpha*(R + 0.5*a)*h;
		faces[bottom].FaceSquare = alpha*(R - 0.5*a)*h;
		faces[left].FaceSquare = alpha*R*a;
		faces[right].FaceSquare = alpha*R*a;

		//add two new faces
		Face f1, f2;

		f1.GlobalIndex = ++MaxFaceIndex;
		f2.GlobalIndex = ++MaxFaceIndex;
		f1.isExternal = true;
		f2.isExternal = true;
		f1.FaceCell_1 = c.GlobalIndex;
		f2.FaceCell_1 = c.GlobalIndex;
		f1.FaceSquare = c.CellVolume;
		f2.FaceSquare = c.CellVolume;
		Vector n;
		n.x = 0;
		n.y = -sin(0.5*alpha);
		n.z = cos(0.5*alpha);
		f1.FaceNormal = n;
		n.z = -cos(0.5*alpha);
		f2.FaceNormal = n;
		Vector FaceCenter = c.CellCenter;
		FaceCenter.y = R*cos(0.5*alpha);
		FaceCenter.z = R*sin(0.5*alpha);
		f1.FaceCenter = FaceCenter;
		FaceCenter.z *= -1;
		f2.FaceCenter = FaceCenter;
		f1.BCMarker = Periodic_BC_Marker;
		f2.BCMarker = Periodic_BC_Marker;
		faces.add(f1);
		faces.add(f2);
		c.Faces.push_back(f1.GlobalIndex);
		c.Faces.push_back(f2.GlobalIndex);
		
		//compute cell volume
		c.CellVolume = alpha*R*a*h;
		if(c.CellVolume<=0) std::cout << "Negative Volume" << '\n';
		
	};
	return;
};

void Grid::ToTecPlotGrid(char *filename = "grid.dat")
{	
	std::set<std::pair<int,int>> l_faces;	//set of edges (linesegments)
	std::vector<Face*> fcs = faces.getLocalNodes();
	//fill set of all edges (linesegments)
	for(int i=0; i<fcs.size(); i++)
	{
		std::pair<int,int> line;
		Face* f = fcs[i];
		for(int j=1; j<f->FaceNodes.size(); j++)
		{
			line.first = f->FaceNodes[j];
			line.second = f->FaceNodes[j-1];
			if(l_faces.find(line)!=l_faces.cend()) continue;
			line.first = f->FaceNodes[j-1];
			line.second = f->FaceNodes[j];
			if(l_faces.find(line)!=l_faces.cend()) continue;
			l_faces.insert(line);
		};
		line.first = f->FaceNodes[0];
		line.second = f->FaceNodes[f->FaceNodes.size()-1];
		if(l_faces.find(line)!=l_faces.cend()) continue;
		line.first = f->FaceNodes[f->FaceNodes.size()-1];
		line.second = f->FaceNodes[0];
		if(l_faces.find(line)!=l_faces.cend()) continue;
		l_faces.insert(line);
	};
	//write grid in tecplot format
	std::ofstream ofs(filename);
	ofs<<"VARIABLES= \"X\", \"Y\", \"Z\"\n";
	ofs<<"ZONE T=\"D\"\n";
	ofs<<"N=" << nodes.size() << ", E=" << l_faces.size() <<", F=FEBLOCK," << " ET=LINESEG\n";
	ofs<<"VARLOCATION = (NODAL, NODAL, NODAL)\n";

	//Map all node global indexes to natural numbers
	std::map<int,int> toNaturalIndex;
	std::set<int> nodeIndexes = nodes.getAllIndexes();
	int counter = 1;
	for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;
	std::vector<Node*> nds = nodes.getLocalNodes();

	//X
	for (int i = 0; i<nodes.size(); i++) {
		ofs<<nds[i]->P.x<<"\n";
	};
	//Y
	for (int i = 0; i<nodes.size(); i++) {
		ofs<<nds[i]->P.y<<"\n";
	};
	//Z
	for (int i = 0; i<nodes.size(); i++) {
		ofs<<nds[i]->P.z<<"\n";
	};
	
	//Connectivity list
	for (std::set<std::pair<int, int>>::iterator it = l_faces.begin(); it!= l_faces.end(); it++) {
		ofs << toNaturalIndex[it->first] << " " << toNaturalIndex[it->second] << '\n';
	};
	ofs.close();

	return;
};

#endif


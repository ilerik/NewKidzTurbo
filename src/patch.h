#ifndef TURBO_PATCH
#define TURBO_PATCH
//
////FOr tesst commit
//
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include "alglibmisc.h"
#include "cgnslib.h"

//each patch describe one boundary
class Patch 
{	
public:	
	int boundaryMarker;
	std::string name;
	bool isFixed;

	//Patch constructor
	Patch() {};
	Patch(std::string _name, int _bcMarker) : name(_name), boundaryMarker(_bcMarker) {};

	//Patch faces and nodes
	std::set<int> faces_idx;
	std::set<int> nodes_idx;

	//connectivity list for periodic condition
	std::map<int, int> periodic_face;

	//KDTree structure and related operations	
	alglib::kdtree kdt;

	void addNode(int nodeIdx) {
		nodes_idx.insert(nodeIdx);
	};

	void addFace(int faceIdx) {
		faces_idx.insert(faceIdx);
	};
};

#endif


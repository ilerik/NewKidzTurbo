#ifndef TURBO_GRIDLOADING_TRIANGLELOAD
#define TURBO_GRIDLOADING_TRIANGLELOAD

#include "stdlib.h"
#include "grid.h"
#include "cgnslib.h"
#include "iostream"

Grid Load2DTriangleGrid(std::string fname) {
	Grid grid;	

	// grid info
	grid.gridInfo.CellDimensions = 2;
	grid.gridInfo.GridDimensions = 2;

	//open input
	std::ifstream fin(fname);

	// read number of nodes and cells
	int nNodes;
	int nCells;
	int nWhatever;
	fin>>nNodes>>nCells>>nWhatever;

	// read nodes coords	
	for (int i = 0; i<nNodes; i++) {
		Node newNode;
		newNode.GlobalIndex = i+1;
		fin>>newNode.P.x;
		fin>>newNode.P.y;
		newNode.P.z = 0;
		grid.nodes.add(newNode);
	};

	// read connectivity list and create faces and cells	
	std::map<std::set<int>, int> faces;
	int faceGlobalIndex = 1;	
	int N1, N2, N3;
	for (int ind = 0; ind<nCells; ind++) {
		Cell newCell;
		newCell.GlobalIndex = ind;
		fin>>N1>>N2>>N3;
		newCell.Nodes.push_back(N1);					
		newCell.Nodes.push_back(N2);					
		newCell.Nodes.push_back(N3);					
		newCell.CGNSType = TRI_3;	
		//Fill in geometric properties
		grid.ComputeGeometricProperties(newCell);

		//Obtain faces based on cell type
		std::vector<Face> newFaces = grid.ObtainFaces(newCell);
		for (int i = 0; i<newFaces.size(); i++) {				
			std::set<int> idx; for (int j = 0; j<newFaces[i].FaceNodes.size(); j++) idx.insert(newFaces[i].FaceNodes[j]);

			std::map<std::set<int>, int>::iterator it = faces.find(idx);
			if (it != faces.end()) {
				//Second time
				Face& face = grid.faces[it->second];
				face.FaceCell_2 = newCell.GlobalIndex;
				face.isExternal = false;				
				newCell.Faces.push_back(face.GlobalIndex);
			} else {
				//First time	
				Face& face = newFaces[i];
				face.GlobalIndex = faceGlobalIndex++;
				face.FaceCell_1 = newCell.GlobalIndex;
				face.isExternal = true;
				grid.ComputeGeometricProperties(face);	

				//Adjust face normal to point outside the cell
				double df = face.FaceNormal * (newCell.CellCenter - face.FaceCenter);
				if ( df > 0) face.FaceNormal *= -1;

				grid.faces.add(face);
				faces[idx] = face.GlobalIndex;
				newCell.Faces.push_back(face.GlobalIndex);
			};
		};

		grid.cells.add(newCell);
	};

	return grid;
};

#endif
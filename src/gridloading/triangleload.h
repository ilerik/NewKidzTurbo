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
				if (df > 0) face.FaceNormal *= -1;

				grid.faces.add(face);
				faces[idx] = face.GlobalIndex;
				newCell.Faces.push_back(face.GlobalIndex);
			};
		};

		grid.cells.add(newCell);
	};

	return grid;
};

//convert cgns grid to full format
void ConvertGrid(char* file_name)
{
	std::string gridFile = "D:\\Projects\\NewKidzTurbo\\Grids\\SimpleCircle.cgns";
	Grid grid = LoadCGNSGrid(gridFile);

	//write full info about grid
	std::ofstream ofs(file_name);
	std::vector<Node*> nodes = grid.nodes.getLocalNodes();
	std::vector<Cell*> cells = grid.cells.getLocalNodes();
	std::vector<Face*> faces = grid.faces.getLocalNodes();

	//write nodes info
	ofs << nodes.size() << '\n';
	for(int i=0; i<nodes.size(); i++)
	{
		ofs << nodes[i]->GlobalIndex << '\n';
		ofs << nodes[i]->P.x << ' ';
		ofs << nodes[i]->P.y << ' ';
		ofs << nodes[i]->P.z << '\n';
	};

	//write cells info
	ofs << cells.size() << '\n';
	for(int i=0; i<cells.size(); i++)
	{
		ofs << cells[i]->GlobalIndex << '\n';
		ofs << cells[i]->CellCenter.x << ' ';
		ofs << cells[i]->CellCenter.y << ' ';
		ofs << cells[i]->CellCenter.z << '\n';
		ofs << cells[i]->CellVolume << '\n';
		int faces_size = cells[i]->Faces.size();
		ofs << faces_size;
		for(int j=0; j<faces_size; j++) ofs << ' ' << cells[i]->Faces[j];
		ofs << '\n';
		int nodes_size = cells[i]->Nodes.size();
		ofs << nodes_size;
		for(int j=0; j<nodes_size; j++) ofs << ' ' << cells[i]->Nodes[j];
		ofs << '\n';
	};

	//write faces info
	ofs << faces.size() << '\n';
	for(int i=0; i<faces.size(); i++)
	{
		ofs << faces[i]->GlobalIndex << '\n';
		ofs << faces[i]->FaceCenter.x << ' ';
		ofs << faces[i]->FaceCenter.y << ' ';
		ofs << faces[i]->FaceCenter.z << '\n';
		ofs << faces[i]->isExternal << '\n';
		ofs << faces[i]->BCMarker << '\n';
		ofs << faces[i]->FaceSquare << '\n';
		ofs << faces[i]->FaceNormal.x << ' ';
		ofs << faces[i]->FaceNormal.y << ' ';
		ofs << faces[i]->FaceNormal.z << '\n';
		ofs << faces[i]->FaceCell_1 << ' ';
		ofs << faces[i]->FaceCell_2 << '\n';
		int nodes_size = faces[i]->FaceNodes.size();
		ofs << nodes_size;
		for(int j=0; j<nodes_size; j++) ofs << ' ' << faces[i]->FaceNodes[j];
		ofs << '\n';
	};
	ofs.close();

	return;
};

#endif
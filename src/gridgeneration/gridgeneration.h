#ifndef TURBO_GRIDGENERATION
#define TURBO_GRIDGENERATION

#include "gengrid2D.h"
#include "gengrid1D.h"


//check grid correctness
bool check_grid(Grid g)
{
	bool good_grid = true;
	std::vector<Face*> faces = g.faces.getLocalNodes();
	std::vector<Cell*> cells = g.cells.getLocalNodes();

	std::cout<< "cells number equals " << cells.size() << '\n';
	std::cout<< "nodes number equals " << g.nodes.size() << '\n';

	//Compute the sum of the directed areas for each cell
	Vector res = Vector(0,0,0);
	for(int i=0; i<cells.size(); i++)
	{
		Cell& c = *cells[i];
		for(int j=0; j<c.Faces.size(); j++)
		{
			Face& f = g.faces[c.Faces[j]];
			if(f.FaceCell_1==c.GlobalIndex) res += f.FaceSquare*f.FaceNormal;
			else res -= f.FaceSquare*f.FaceNormal;
		};
	};
	if(res.mod()!=0){
	std::cout << "Sum of the directed areas for each cell doesn't equal to zero " << res.mod() << '\n';
	good_grid = false;
	};

	//Global sum check for boundary face vector
	res = Vector(0,0,0);
	for(int i=0; i<faces.size(); i++)
	{
		Face& f = *faces[i];
		if(f.isExternal) res += f.FaceSquare*f.FaceNormal;
	};
	if(res.mod()!=0){
	std::cout << "Global sum of boundary directed areas doesn't equal to zero " << res.mod() << '\n';
	good_grid = false;
	};

	//Compute the sum of the directed areas for each inner cell
	res = Vector(0,0,0);
	for(int i=0; i<cells.size(); i++)
	{
		Cell& c = *cells[i];
		for(int j=0; j<c.Faces.size(); j++)
		{
			Face& f = g.faces[c.Faces[j]];
			if(f.isExternal) continue;
			if(f.FaceCell_1==c.GlobalIndex) res += f.FaceSquare*f.FaceNormal;
			else res -= f.FaceSquare*f.FaceNormal;
		};
	};
	if(res.mod()!=0){
	std::cout << "Sum of the directed areas for each cell doesn't equal to zero " << res.mod() << '\n';
	good_grid = false;
	};

	//Check elements volumes and global one
	double Volume = 0;
	for(int i=0; i<cells.size(); i++)
	{
		Cell& c = *cells[i];
		if(c.CellVolume<=0){
			std::cout << "Some cells have unpositive volume\n";
			good_grid = false;
			break;
		};
		Volume += c.CellVolume;
	};
	if(Volume<=0){
	std::cout << "Global Volume of the grid is negative number or zero" << res.mod() << '\n';
	good_grid = false;
	};
	std::cout << "Volume of the grid equals " << Volume << '\n';
	std::getchar();

	return good_grid;
};


//generate 3D uniform rectangular grid
Grid GenGrid3D(int N, int M, int K, double size_x, double size_y, double size_z)
{
	Grid g;
	g.gridInfo.GridDimensions = 3;
	N++;
	M++;
	K++;
	//Generate nodes
	std::vector<double> x_p(N);
	std::vector<double> y_p(M);
	std::vector<double> z_p(K);
	//double h_x = 1.0 * size_x * (1.0 -q_x) / (1.0 - pow(q_x, N-1));
	//double h_y = 1.0 * size_y * (1.0 -q_y) / (1.0 - pow(q_y, M-1));	
	for (int i = 0; i<N; i++) {		
		x_p[i] = (size_x*i) / (N-1);		
	};

	for (int i = 0; i<M; i++) {		
		y_p[i] = (size_y*i) / (M-1);		
	};

	for (int i = 0; i<K; i++) {		
		z_p[i] = (size_z*i) / (K-1);		
	};
	
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<M; j++) {
			for (int k = 0; k<K; k++) {
			Node new_node;
			new_node.GlobalIndex = i + j * N + k * M * N;
			new_node.P.x = x_p[i];
			new_node.P.y = y_p[j];
			new_node.P.z = z_p[k];
			g.nodes.add(new_node);
			};
		};
	};

	//Generate cells
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			for (int k = 0; k<(K-1); k++) {
			Cell c;
			c.GlobalIndex = i + j * (N-1) + k * (N-1) * (M-1);
			c.Nodes.resize(8);
			c.Nodes[0] = i + j * N + k * N * M;
			c.Nodes[1] = (i+1) + j * N + k * N * M;
			c.Nodes[2] = (i+1) + (j+1)*N + k * N * M;
			c.Nodes[3] = i + (j+1)*N + k * N * M;
			c.Nodes[4] = i + j * N + (k+1) * N * M;
			c.Nodes[5] = (i+1) + j * N + (k+1) * N * M;
			c.Nodes[6] = (i+1) + (j+1)*N + (k+1) * N * M;
			c.Nodes[7] = i + (j+1)*N + (k+1) * N * M;
			c.CellCenter = Vector(0,0,0);
			for (int k = 0; k<c.Nodes.size(); k++) {
				c.CellCenter =  c.CellCenter + g.nodes[c.Nodes[k]].P;
			};
			c.CellCenter = (1.0/c.Nodes.size()) * c.CellCenter;		
			c.CellVolume = ((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[0]].P).mod()) * ((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[2]].P).mod())*((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[5]].P).mod());
			c.Faces.resize(6);
			c.CGNSType = HEXA_8;
			c.CellHSize = 0;
			g.cells.add(c);
			};
		};
	};

	//Generate Faces
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			for (int k = 0; k<(K-1); k++) {
			Face nf_bottom, nf_front, nf_left;			// three new faces
			int CellIndex = i + j * (N-1) + k * (N-1) * (M-1);	// our cell

			//write Global Index
			nf_bottom.GlobalIndex = 3 * CellIndex;
			nf_front.GlobalIndex = 3 * CellIndex + 1;
			nf_left.GlobalIndex = 3 * CellIndex + 2;

			//write CGNS element type
			nf_bottom.CGNSType == QUAD_4;
			nf_front.CGNSType == QUAD_4;
			nf_left.CGNSType == QUAD_4;

			//write BCMarker
			nf_bottom.BCMarker = 0;
			nf_front.BCMarker = 0;
			nf_left.BCMarker = 0;
			if(i==0) nf_left.BCMarker = 1;	//left edge - 1, right edge - 2
			if(j==0) nf_front.BCMarker = 3;	//front edge - 3, back edge - 4
			if(k==0) nf_bottom.BCMarker = 5;	//bottom edge - 5, top edge - 6		

			//isExternal
			if(nf_bottom.BCMarker>0) nf_bottom.isExternal = 1;
			else nf_bottom.isExternal = 0;
			if(nf_front.BCMarker>0) nf_front.isExternal = 1;
			else nf_front.isExternal = 0;
			if(nf_left.BCMarker>0) nf_left.isExternal = 1;
			else nf_left.isExternal = 0;

			//set cell-face connection
			nf_bottom.FaceCell_1 = CellIndex;
			if(nf_bottom.isExternal==0) nf_bottom.FaceCell_2 = i + j * (N-1) + (k-1) * (N-1) * (M-1);
			else nf_bottom.FaceCell_2 = -1;

			nf_front.FaceCell_1 = CellIndex;
			if(nf_front.isExternal==0) nf_front.FaceCell_2 = i + (j-1) * (N-1) + k * (N-1) * (M-1);
			else nf_front.FaceCell_2 = -1;

			nf_left.FaceCell_1 = CellIndex;
			if(nf_left.isExternal==0) nf_left.FaceCell_2 = (i-1) + j * (N-1) + k * (N-1) * (M-1);
			else nf_left.FaceCell_2 = -1;

			g.cells[CellIndex].Faces[0] = nf_bottom.GlobalIndex;
			g.cells[CellIndex].Faces[1] = nf_front.GlobalIndex; 
			g.cells[CellIndex].Faces[2] = nf_left.GlobalIndex;
			if(nf_bottom.BCMarker==0) g.cells[i + j * (N-1) + (k-1) * (N-1) * (M-1)].Faces[3] = nf_bottom.GlobalIndex;
			if(nf_front.BCMarker==0) g.cells[i + (j-1) * (N-1) + k * (N-1) * (M-1)].Faces[4] = nf_front.GlobalIndex;
			if(nf_left.BCMarker==0) g.cells[(i-1) + j * (N-1) + k * (N-1) * (M-1)].Faces[5] = nf_left.GlobalIndex;

			//write normals
			nf_bottom.FaceNormal = Vector(0,0,-1);
			nf_front.FaceNormal = Vector(0,-1,0);
			nf_left.FaceNormal = Vector(-1,0,0);

			//write nodes
			nf_bottom.FaceNodes.resize(4);
			nf_front.FaceNodes.resize(4);
			nf_left.FaceNodes.resize(4);
			nf_bottom.FaceNodes[0] = g.cells[CellIndex].Nodes[0];
			nf_bottom.FaceNodes[1] = g.cells[CellIndex].Nodes[1];
			nf_bottom.FaceNodes[2] = g.cells[CellIndex].Nodes[2];
			nf_bottom.FaceNodes[3] = g.cells[CellIndex].Nodes[3];
			nf_front.FaceNodes[0] = g.cells[CellIndex].Nodes[0];
			nf_front.FaceNodes[1] = g.cells[CellIndex].Nodes[1];
			nf_front.FaceNodes[2] = g.cells[CellIndex].Nodes[5];
			nf_front.FaceNodes[3] = g.cells[CellIndex].Nodes[4];
			nf_left.FaceNodes[0] = g.cells[CellIndex].Nodes[0];
			nf_left.FaceNodes[1] = g.cells[CellIndex].Nodes[3];
			nf_left.FaceNodes[2] = g.cells[CellIndex].Nodes[7];
			nf_left.FaceNodes[3] = g.cells[CellIndex].Nodes[4];

			//calc FaceCenter
			nf_bottom.FaceCenter = 0.5*(g.nodes[nf_bottom.FaceNodes[0]].P + g.nodes[nf_bottom.FaceNodes[2]].P);
			nf_front.FaceCenter = 0.5*(g.nodes[nf_front.FaceNodes[0]].P + g.nodes[nf_front.FaceNodes[2]].P);
			nf_left.FaceCenter = 0.5*(g.nodes[nf_left.FaceNodes[0]].P + g.nodes[nf_left.FaceNodes[2]].P);

			//calc squares
			Vector a = g.nodes[nf_bottom.FaceNodes[1]].P - g.nodes[nf_bottom.FaceNodes[0]].P;
			Vector b = g.nodes[nf_bottom.FaceNodes[3]].P - g.nodes[nf_bottom.FaceNodes[0]].P;
			nf_bottom.FaceSquare = (a & b).mod();

			a = g.nodes[nf_front.FaceNodes[1]].P - g.nodes[nf_front.FaceNodes[0]].P;
			b = g.nodes[nf_front.FaceNodes[3]].P - g.nodes[nf_front.FaceNodes[0]].P;
			nf_front.FaceSquare = (a & b).mod();

			a = g.nodes[nf_left.FaceNodes[1]].P - g.nodes[nf_left.FaceNodes[0]].P;
			b = g.nodes[nf_left.FaceNodes[3]].P - g.nodes[nf_left.FaceNodes[0]].P;
			nf_left.FaceSquare = (a & b).mod();

			//write to the grid
			g.faces.add(nf_bottom);
			g.faces.add(nf_front);
			g.faces.add(nf_left);
			};			
		};
	};

	int FIndex = 3*(N-1)*(M-1)*(K-1);
	//write top faces
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			Face face_top;
			face_top.CGNSType = QUAD_4;
			int CellIndex = i + j*(N-1) + (K-2)*(N-1)*(M-1);
			face_top.GlobalIndex = FIndex + i + j*(N-1);
			face_top.BCMarker = 6;
			face_top.isExternal = 1;
			face_top.FaceCell_1 = CellIndex;
			face_top.FaceCell_2 = -1;
			g.cells[CellIndex].Faces[3] = face_top.GlobalIndex;
			face_top.FaceNodes.resize(4);
			face_top.FaceNodes[0] = g.cells[CellIndex].Nodes[4];
			face_top.FaceNodes[1] = g.cells[CellIndex].Nodes[5];
			face_top.FaceNodes[2] = g.cells[CellIndex].Nodes[6];
			face_top.FaceNodes[3] = g.cells[CellIndex].Nodes[7];
			face_top.FaceCenter = 0.5*(g.nodes[face_top.FaceNodes[0]].P + g.nodes[face_top.FaceNodes[2]].P);
			//square
			Vector a = g.nodes[face_top.FaceNodes[1]].P - g.nodes[face_top.FaceNodes[0]].P;
			Vector b = g.nodes[face_top.FaceNodes[3]].P - g.nodes[face_top.FaceNodes[0]].P;
			face_top.FaceSquare = (a & b).mod();

			face_top.FaceNormal = Vector(0, 0, 1);
			g.faces.add(face_top);
		};
	};

	FIndex += (M-1)*(N-1);
	//write back faces
	for (int i = 0; i<(N-1); i++) {
		for (int k = 0; k<(K-1); k++) {
			Face face_back;
			face_back.CGNSType = QUAD_4;
			int CellIndex = i + (M-2)*(N-1) + k*(N-1)*(M-1);
			face_back.GlobalIndex = FIndex + i + k*(N-1);
			face_back.BCMarker = 4;
			face_back.isExternal = 1;
			face_back.FaceCell_1 = CellIndex;
			face_back.FaceCell_2 = -1;
			g.cells[CellIndex].Faces[4] = face_back.GlobalIndex;
			face_back.FaceNodes.resize(4);
			face_back.FaceNodes[0] = g.cells[CellIndex].Nodes[3];
			face_back.FaceNodes[1] = g.cells[CellIndex].Nodes[2];
			face_back.FaceNodes[2] = g.cells[CellIndex].Nodes[6];
			face_back.FaceNodes[3] = g.cells[CellIndex].Nodes[7];
			face_back.FaceCenter = 0.5*(g.nodes[face_back.FaceNodes[0]].P + g.nodes[face_back.FaceNodes[2]].P);
			//square
			Vector a = g.nodes[face_back.FaceNodes[1]].P - g.nodes[face_back.FaceNodes[0]].P;
			Vector b = g.nodes[face_back.FaceNodes[3]].P - g.nodes[face_back.FaceNodes[0]].P;
			face_back.FaceSquare = (a & b).mod();

			face_back.FaceNormal = Vector(0, 1, 0);
			g.faces.add(face_back);
		};
	};

	FIndex += (K-1)*(N-1);
	//write right faces
	for (int j = 0; j<(M-1); j++) {
		for (int k = 0; k<(K-1); k++) {
			Face face_right;
			face_right.CGNSType = QUAD_4;
			int CellIndex = (N-2) + j*(N-1) + k*(N-1)*(M-1);
			face_right.GlobalIndex = FIndex + j + k*(M-1);
			face_right.BCMarker = 2;
			face_right.isExternal = 1;
			face_right.FaceCell_1 = CellIndex;
			face_right.FaceCell_2 = -1;
			g.cells[CellIndex].Faces[5] = face_right.GlobalIndex;
			face_right.FaceNodes.resize(4);
			face_right.FaceNodes[0] = g.cells[CellIndex].Nodes[1];
			face_right.FaceNodes[1] = g.cells[CellIndex].Nodes[2];
			face_right.FaceNodes[2] = g.cells[CellIndex].Nodes[6];
			face_right.FaceNodes[3] = g.cells[CellIndex].Nodes[5];
			face_right.FaceCenter = 0.5*(g.nodes[face_right.FaceNodes[0]].P + g.nodes[face_right.FaceNodes[2]].P);
			
			//square
			Vector a = g.nodes[face_right.FaceNodes[1]].P - g.nodes[face_right.FaceNodes[0]].P;
			Vector b = g.nodes[face_right.FaceNodes[3]].P - g.nodes[face_right.FaceNodes[0]].P;
			face_right.FaceSquare = (a & b).mod();

			face_right.FaceNormal = Vector(1, 0, 0);
			g.faces.add(face_right);
		};
	};

	//Add boundary names
	g.addPatch("left", 1);
	g.addPatch("right", 2);
	g.addPatch("front", 3);
	g.addPatch("back", 4);
	g.addPatch("bottom", 5);
	g.addPatch("top", 6);
	g.ConstructAndCheckPatches();	

	//periodic boundary condition part
	//top and bottom periodic
	FIndex = 3*(N-1)*(M-1)*(K-1);
	//write top faces
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			int face_top_idx = FIndex + i + j*(N-1);	//index of top cell
			int CellIndex = i + j * (N-1);				// bottom cell index
			int face_bot_idx = 3*CellIndex;				//bottom face index
			g.patches[6].periodic_face[face_top_idx] = face_bot_idx;
			g.patches[5].periodic_face[face_bot_idx] = face_top_idx;
		};
	};
	//front and back periodic
	FIndex += (M-1)*(N-1);
	//write back faces
	for (int i = 0; i<(N-1); i++) {
		for (int k = 0; k<(K-1); k++) {
			int face_back_idx = FIndex + i + k*(N-1);			//index of back face
			int CellIndex = i + k * (N-1) * (M-1);				//index of corresponding front cell
			int face_front_idx = 3 * CellIndex + 1;				//index of front face
			g.patches[4].periodic_face[face_back_idx] = face_front_idx;
			g.patches[3].periodic_face[face_front_idx] = face_back_idx;
		};
	};
	//left and right periodic
	FIndex += (K-1)*(N-1);
	//write right faces
	for (int j = 0; j<(M-1); j++) {
		for (int k = 0; k<(K-1); k++) {
			int face_right_idx = FIndex + j + k*(M-1);			//right face index
			int CellIndex = j*(N-1) + k*(N-1)*(M-1);			//left cell index
			int face_left_idx = 3*CellIndex + 2;				//left face index
			g.patches[2].periodic_face[face_right_idx] = face_left_idx;
			g.patches[1].periodic_face[face_left_idx] = face_right_idx;
		};
	};


	return g;
};
Grid GenGrid3D(int N, int M, int K, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
	Grid g;
	g.gridInfo.GridDimensions = 3;
	N++;
	M++;
	K++;
	//Generate nodes
	std::vector<double> x_p(N);
	std::vector<double> y_p(M);
	std::vector<double> z_p(K);

	double size_x = x_max - x_min;
	double size_y = y_max - y_min;
	double size_z = z_max - z_min;
	//double h_x = 1.0 * size_x * (1.0 -q_x) / (1.0 - pow(q_x, N-1));
	//double h_y = 1.0 * size_y * (1.0 -q_y) / (1.0 - pow(q_y, M-1));	
	for (int i = 0; i<N; i++) {		
		x_p[i] = x_min + (size_x*i) / (N-1);		
	};

	for (int i = 0; i<M; i++) {		
		y_p[i] = y_min + (size_y*i) / (M-1);		
	};

	for (int i = 0; i<K; i++) {		
		z_p[i] = z_min + (size_z*i) / (K-1);		
	};
	
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<M; j++) {
			for (int k = 0; k<K; k++) {
			Node new_node;
			new_node.GlobalIndex = i + j * N + k * M * N;
			new_node.P.x = x_p[i];
			new_node.P.y = y_p[j];
			new_node.P.z = z_p[k];
			g.nodes.add(new_node);
			};
		};
	};

	//Generate cells
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			for (int k = 0; k<(K-1); k++) {
			Cell c;
			c.GlobalIndex = i + j * (N-1) + k * (N-1) * (M-1);
			c.Nodes.resize(8);
			c.Nodes[0] = i + j * N + k * N * M;
			c.Nodes[1] = (i+1) + j * N + k * N * M;
			c.Nodes[2] = (i+1) + (j+1)*N + k * N * M;
			c.Nodes[3] = i + (j+1)*N + k * N * M;
			c.Nodes[4] = i + j * N + (k+1) * N * M;
			c.Nodes[5] = (i+1) + j * N + (k+1) * N * M;
			c.Nodes[6] = (i+1) + (j+1)*N + (k+1) * N * M;
			c.Nodes[7] = i + (j+1)*N + (k+1) * N * M;
			c.CellCenter = Vector(0,0,0);
			for (int k = 0; k<c.Nodes.size(); k++) {
				c.CellCenter =  c.CellCenter + g.nodes[c.Nodes[k]].P;
			};
			c.CellCenter = (1.0/c.Nodes.size()) * c.CellCenter;		
			c.CellVolume = ((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[0]].P).mod()) * ((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[2]].P).mod())*((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[5]].P).mod());
			c.Faces.resize(6);
			c.CGNSType = HEXA_8;
			c.CellHSize = 0;
			g.cells.add(c);
			};
		};
	};

	//Generate Faces
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			for (int k = 0; k<(K-1); k++) {
			Face nf_bottom, nf_front, nf_left;			// three new faces
			int CellIndex = i + j * (N-1) + k * (N-1) * (M-1);	// our cell

			//write Global Index
			nf_bottom.GlobalIndex = 3 * CellIndex;
			nf_front.GlobalIndex = 3 * CellIndex + 1;
			nf_left.GlobalIndex = 3 * CellIndex + 2;

			//write CGNS element type
			nf_bottom.CGNSType == QUAD_4;
			nf_front.CGNSType == QUAD_4;
			nf_left.CGNSType == QUAD_4;

			//write BCMarker
			nf_bottom.BCMarker = 0;
			nf_front.BCMarker = 0;
			nf_left.BCMarker = 0;
			if(i==0) nf_left.BCMarker = 1;	//left edge - 1, right edge - 2
			if(j==0) nf_front.BCMarker = 3;	//front edge - 3, back edge - 4
			if(k==0) nf_bottom.BCMarker = 5;	//bottom edge - 5, top edge - 6		

			//isExternal
			if(nf_bottom.BCMarker>0) nf_bottom.isExternal = 1;
			else nf_bottom.isExternal = 0;
			if(nf_front.BCMarker>0) nf_front.isExternal = 1;
			else nf_front.isExternal = 0;
			if(nf_left.BCMarker>0) nf_left.isExternal = 1;
			else nf_left.isExternal = 0;

			//set cell-face connection
			nf_bottom.FaceCell_1 = CellIndex;
			if(nf_bottom.isExternal==0) nf_bottom.FaceCell_2 = i + j * (N-1) + (k-1) * (N-1) * (M-1);
			else nf_bottom.FaceCell_2 = -1;

			nf_front.FaceCell_1 = CellIndex;
			if(nf_front.isExternal==0) nf_front.FaceCell_2 = i + (j-1) * (N-1) + k * (N-1) * (M-1);
			else nf_front.FaceCell_2 = -1;

			nf_left.FaceCell_1 = CellIndex;
			if(nf_left.isExternal==0) nf_left.FaceCell_2 = (i-1) + j * (N-1) + k * (N-1) * (M-1);
			else nf_left.FaceCell_2 = -1;

			g.cells[CellIndex].Faces[0] = nf_bottom.GlobalIndex;
			g.cells[CellIndex].Faces[1] = nf_front.GlobalIndex; 
			g.cells[CellIndex].Faces[2] = nf_left.GlobalIndex;
			if(nf_bottom.BCMarker==0) g.cells[i + j * (N-1) + (k-1) * (N-1) * (M-1)].Faces[3] = nf_bottom.GlobalIndex;
			if(nf_front.BCMarker==0) g.cells[i + (j-1) * (N-1) + k * (N-1) * (M-1)].Faces[4] = nf_front.GlobalIndex;
			if(nf_left.BCMarker==0) g.cells[(i-1) + j * (N-1) + k * (N-1) * (M-1)].Faces[5] = nf_left.GlobalIndex;

			//write normals
			nf_bottom.FaceNormal = Vector(0,0,-1);
			nf_front.FaceNormal = Vector(0,-1,0);
			nf_left.FaceNormal = Vector(-1,0,0);

			//write nodes
			nf_bottom.FaceNodes.resize(4);
			nf_front.FaceNodes.resize(4);
			nf_left.FaceNodes.resize(4);
			nf_bottom.FaceNodes[0] = g.cells[CellIndex].Nodes[0];
			nf_bottom.FaceNodes[1] = g.cells[CellIndex].Nodes[1];
			nf_bottom.FaceNodes[2] = g.cells[CellIndex].Nodes[2];
			nf_bottom.FaceNodes[3] = g.cells[CellIndex].Nodes[3];
			nf_front.FaceNodes[0] = g.cells[CellIndex].Nodes[0];
			nf_front.FaceNodes[1] = g.cells[CellIndex].Nodes[1];
			nf_front.FaceNodes[2] = g.cells[CellIndex].Nodes[5];
			nf_front.FaceNodes[3] = g.cells[CellIndex].Nodes[4];
			nf_left.FaceNodes[0] = g.cells[CellIndex].Nodes[0];
			nf_left.FaceNodes[1] = g.cells[CellIndex].Nodes[3];
			nf_left.FaceNodes[2] = g.cells[CellIndex].Nodes[7];
			nf_left.FaceNodes[3] = g.cells[CellIndex].Nodes[4];

			//calc FaceCenter
			nf_bottom.FaceCenter = 0.5*(g.nodes[nf_bottom.FaceNodes[0]].P + g.nodes[nf_bottom.FaceNodes[2]].P);
			nf_front.FaceCenter = 0.5*(g.nodes[nf_front.FaceNodes[0]].P + g.nodes[nf_front.FaceNodes[2]].P);
			nf_left.FaceCenter = 0.5*(g.nodes[nf_left.FaceNodes[0]].P + g.nodes[nf_left.FaceNodes[2]].P);

			//calc squares
			Vector a = g.nodes[nf_bottom.FaceNodes[1]].P - g.nodes[nf_bottom.FaceNodes[0]].P;
			Vector b = g.nodes[nf_bottom.FaceNodes[3]].P - g.nodes[nf_bottom.FaceNodes[0]].P;
			nf_bottom.FaceSquare = (a & b).mod();

			a = g.nodes[nf_front.FaceNodes[1]].P - g.nodes[nf_front.FaceNodes[0]].P;
			b = g.nodes[nf_front.FaceNodes[3]].P - g.nodes[nf_front.FaceNodes[0]].P;
			nf_front.FaceSquare = (a & b).mod();

			a = g.nodes[nf_left.FaceNodes[1]].P - g.nodes[nf_left.FaceNodes[0]].P;
			b = g.nodes[nf_left.FaceNodes[3]].P - g.nodes[nf_left.FaceNodes[0]].P;
			nf_left.FaceSquare = (a & b).mod();

			//write to the grid
			g.faces.add(nf_bottom);
			g.faces.add(nf_front);
			g.faces.add(nf_left);
			};			
		};
	};

	int FIndex = 3*(N-1)*(M-1)*(K-1);
	//write top faces
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			Face face_top;
			face_top.CGNSType = QUAD_4;
			int CellIndex = i + j*(N-1) + (K-2)*(N-1)*(M-1);
			face_top.GlobalIndex = FIndex + i + j*(N-1);
			face_top.BCMarker = 6;
			face_top.isExternal = 1;
			face_top.FaceCell_1 = CellIndex;
			face_top.FaceCell_2 = -1;
			g.cells[CellIndex].Faces[3] = face_top.GlobalIndex;
			face_top.FaceNodes.resize(4);
			face_top.FaceNodes[0] = g.cells[CellIndex].Nodes[4];
			face_top.FaceNodes[1] = g.cells[CellIndex].Nodes[5];
			face_top.FaceNodes[2] = g.cells[CellIndex].Nodes[6];
			face_top.FaceNodes[3] = g.cells[CellIndex].Nodes[7];
			face_top.FaceCenter = 0.5*(g.nodes[face_top.FaceNodes[0]].P + g.nodes[face_top.FaceNodes[2]].P);
			//square
			Vector a = g.nodes[face_top.FaceNodes[1]].P - g.nodes[face_top.FaceNodes[0]].P;
			Vector b = g.nodes[face_top.FaceNodes[3]].P - g.nodes[face_top.FaceNodes[0]].P;
			face_top.FaceSquare = (a & b).mod();

			face_top.FaceNormal = Vector(0, 0, 1);
			g.faces.add(face_top);
		};
	};

	FIndex += (M-1)*(N-1);
	//write back faces
	for (int i = 0; i<(N-1); i++) {
		for (int k = 0; k<(K-1); k++) {
			Face face_back;
			face_back.CGNSType = QUAD_4;
			int CellIndex = i + (M-2)*(N-1) + k*(N-1)*(M-1);
			face_back.GlobalIndex = FIndex + i + k*(N-1);
			face_back.BCMarker = 4;
			face_back.isExternal = 1;
			face_back.FaceCell_1 = CellIndex;
			face_back.FaceCell_2 = -1;
			g.cells[CellIndex].Faces[4] = face_back.GlobalIndex;
			face_back.FaceNodes.resize(4);
			face_back.FaceNodes[0] = g.cells[CellIndex].Nodes[3];
			face_back.FaceNodes[1] = g.cells[CellIndex].Nodes[2];
			face_back.FaceNodes[2] = g.cells[CellIndex].Nodes[6];
			face_back.FaceNodes[3] = g.cells[CellIndex].Nodes[7];
			face_back.FaceCenter = 0.5*(g.nodes[face_back.FaceNodes[0]].P + g.nodes[face_back.FaceNodes[2]].P);
			//square
			Vector a = g.nodes[face_back.FaceNodes[1]].P - g.nodes[face_back.FaceNodes[0]].P;
			Vector b = g.nodes[face_back.FaceNodes[3]].P - g.nodes[face_back.FaceNodes[0]].P;
			face_back.FaceSquare = (a & b).mod();

			face_back.FaceNormal = Vector(0, 1, 0);
			g.faces.add(face_back);
		};
	};

	FIndex += (K-1)*(N-1);
	//write right faces
	for (int j = 0; j<(M-1); j++) {
		for (int k = 0; k<(K-1); k++) {
			Face face_right;
			face_right.CGNSType = QUAD_4;
			int CellIndex = (N-2) + j*(N-1) + k*(N-1)*(M-1);
			face_right.GlobalIndex = FIndex + j + k*(M-1);
			face_right.BCMarker = 2;
			face_right.isExternal = 1;
			face_right.FaceCell_1 = CellIndex;
			face_right.FaceCell_2 = -1;
			g.cells[CellIndex].Faces[5] = face_right.GlobalIndex;
			face_right.FaceNodes.resize(4);
			face_right.FaceNodes[0] = g.cells[CellIndex].Nodes[1];
			face_right.FaceNodes[1] = g.cells[CellIndex].Nodes[2];
			face_right.FaceNodes[2] = g.cells[CellIndex].Nodes[6];
			face_right.FaceNodes[3] = g.cells[CellIndex].Nodes[5];
			face_right.FaceCenter = 0.5*(g.nodes[face_right.FaceNodes[0]].P + g.nodes[face_right.FaceNodes[2]].P);
			
			//square
			Vector a = g.nodes[face_right.FaceNodes[1]].P - g.nodes[face_right.FaceNodes[0]].P;
			Vector b = g.nodes[face_right.FaceNodes[3]].P - g.nodes[face_right.FaceNodes[0]].P;
			face_right.FaceSquare = (a & b).mod();

			face_right.FaceNormal = Vector(1, 0, 0);
			g.faces.add(face_right);
		};
	};

	//Add boundary names
	g.addPatch("left", 1);
	g.addPatch("right", 2);
	g.addPatch("front", 3);
	g.addPatch("back", 4);
	g.addPatch("bottom", 5);
	g.addPatch("top", 6);
	g.ConstructAndCheckPatches();	

	//periodic boundary condition part
	//top and bottom periodic
	FIndex = 3*(N-1)*(M-1)*(K-1);
	//write top faces
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			int face_top_idx = FIndex + i + j*(N-1);	//index of top cell
			int CellIndex = i + j * (N-1);				// bottom cell index
			int face_bot_idx = 3*CellIndex;				//bottom face index
			g.patches[6].periodic_face[face_top_idx] = face_bot_idx;
			g.patches[5].periodic_face[face_bot_idx] = face_top_idx;
		};
	};
	//front and back periodic
	FIndex += (M-1)*(N-1);
	//write back faces
	for (int i = 0; i<(N-1); i++) {
		for (int k = 0; k<(K-1); k++) {
			int face_back_idx = FIndex + i + k*(N-1);			//index of back face
			int CellIndex = i + k * (N-1) * (M-1);				//index of corresponding front cell
			int face_front_idx = 3 * CellIndex + 1;				//index of front face
			g.patches[4].periodic_face[face_back_idx] = face_front_idx;
			g.patches[3].periodic_face[face_front_idx] = face_back_idx;
		};
	};
	//left and right periodic
	FIndex += (K-1)*(N-1);
	//write right faces
	for (int j = 0; j<(M-1); j++) {
		for (int k = 0; k<(K-1); k++) {
			int face_right_idx = FIndex + j + k*(M-1);			//right face index
			int CellIndex = j*(N-1) + k*(N-1)*(M-1);			//left cell index
			int face_left_idx = 3*CellIndex + 2;				//left face index
			g.patches[2].periodic_face[face_right_idx] = face_left_idx;
			g.patches[1].periodic_face[face_left_idx] = face_right_idx;
		};
	};


	return g;
};

#endif
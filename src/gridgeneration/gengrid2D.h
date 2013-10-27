#ifndef TURBO_GRIDGENERATION_GENGRID2D
#define TURBO_GRIDGENERATION_GENGRID2D

#include "grid.h"

Grid GenGrid2D(int N, int M, double size_x, double size_y, double q_x, double q_y)
{
	//Generate grid topology
	Grid g;
	N++;
	M++;
	std::vector<double> x_p(N);
	std::vector<double> y_p(M);
	double h_x = 1.0 * size_x * (1.0 -q_x) / (1.0 - pow(q_x, N-1));
	double h_y = 1.0 * size_y * (1.0 -q_y) / (1.0 - pow(q_y, M-1));
	for (int i = 0; i<N; i++) {
		if (q_x == 1.0) {
			x_p[i] = (1.0*size_x*i) / (N-1);
		} else {
			x_p[i] = h_x * (1.0 - pow(q_x, i)) / (1.0-q_x);
		};
	};
	x_p[0] = 0;
	x_p[N-1] = size_x;
	for (int i = 0; i<M; i++) {
		if (q_y == 1.0) {
			y_p[i] = (1.0 *size_y*i) / (M-1);
		} else {
			y_p[i] = 1.0 * h_y * (1.0 - pow(q_y, i)) / (1.0-q_y);
		};
	};
	y_p[0] = 0;
	y_p[M-1] = size_y;	
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<M; j++) {
			Node new_node;
			new_node.GlobalIndex = i + j * N;
			new_node.P.x = x_p[i];//(size_x*i) / (N-1);
			new_node.P.y = y_p[j];//(size_y*j) / (M-1);
			new_node.P.z = 0;
			g.nodes.add(new_node); 
		};
	};	
	for (int i = 0; i<(N-1); i++) {
		for (int j = 0; j<(M-1); j++) {
			Cell c;
			c.GlobalIndex = i + j*(N-1);
			c.Nodes.resize(4);
			c.Nodes[0] = i + j*N;
			c.Nodes[1] = i + (j+1)*N;
			c.Nodes[2] = (i+1) + (j+1)*N;
			c.Nodes[3] = (i+1) + j*N;			
			c.CellCenter = Vector(0,0,0);
			for (int k = 0; k<c.Nodes.size(); k++) {
				c.CellCenter =  c.CellCenter + g.nodes[c.Nodes[k]].P;
			};
			c.CellCenter = (1.0/c.Nodes.size()) * c.CellCenter;		
			c.CellVolume = ((g.nodes[c.Nodes[1]].P-g.nodes[c.Nodes[0]].P).mod()) * ((g.nodes[c.Nodes[2]].P-g.nodes[c.Nodes[1]].P).mod());
			c.Faces.resize(4);
			c.Faces[0] = i*(M-1) + j;
			c.Faces[1] = (i+1)*(M-1) + j;
			c.Faces[2] = (i + j*(N-1)) + N*(M-1);
			c.Faces[3] = (i + (j+1)*(N-1)) + N*(M-1);			
			c.CellHSize = 0;
			g.cells.add(c);
		};
	};

	//Vertical first
	for (int i = 0; i<N; i++) {
		for (int j = 0; j<M-1; j++) {
			Face new_face;		
			new_face.GlobalIndex = i*(M-1) + j;
			new_face.FaceNodes.clear();
			new_face.FaceNodes.push_back(i + j * N);
			new_face.FaceNodes.push_back(i + (j+1) * N);
			new_face.FaceSquare = (g.nodes[new_face.FaceNodes[1]].P - g.nodes[new_face.FaceNodes[0]].P).mod();
			new_face.FaceCenter = Vector(0,0,0);
			for (int k = 0; k<new_face.FaceNodes.size(); k++) {
				new_face.FaceCenter = new_face.FaceCenter + g.nodes[new_face.FaceNodes[k]].P;
			};
			new_face.FaceCenter = (1.0/new_face.FaceNodes.size()) * new_face.FaceCenter;
			if (i == 0) {
				//Left border
				new_face.FaceNormal = Vector(-1, 0, 0);
				new_face.isExternal = true;
				new_face.BCMarker = 1;
				new_face.FaceCell_1 = i + j*(N-1);
				new_face.FaceCell_2 = -1;
			};
			if (i == N-1) {
				//Right border
				new_face.FaceNormal = Vector(1, 0, 0);
				new_face.isExternal = true;
				new_face.BCMarker = 2;
				new_face.FaceCell_1 = (i-1) + j*(N-1);
				new_face.FaceCell_2 = -1;
			};
			if ((i!=0) && (i!=N-1)) {
				//Inner face
				new_face.FaceNormal = Vector(1, 0, 0);
				new_face.isExternal = false;
				new_face.BCMarker = -1;
				new_face.FaceCell_1 = (i-1) + j*(N-1);
				new_face.FaceCell_2 = i + j*(N-1);
			};
			//new_face.FaceNormal = new_face.FaceNormal * new_face.FaceSquare;
			g.faces.add(new_face);
		};
	};
	//Horisontal now
	for (int i = 0; i<N-1; i++) {
		for (int j = 0; j<M; j++) {
			Face new_face;		
			new_face.GlobalIndex = N*(M-1) + j * (N-1) + i;
			new_face.FaceNodes.clear();
			new_face.FaceNodes.push_back(i + j * N);
			new_face.FaceNodes.push_back(i+1 + j * N);
			new_face.FaceSquare = (g.nodes[new_face.FaceNodes[1]].P - g.nodes[new_face.FaceNodes[0]].P).mod();
			new_face.FaceCenter = Vector(0,0,0);
			for (int k = 0; k<new_face.FaceNodes.size(); k++) {
				new_face.FaceCenter = new_face.FaceCenter + g.nodes[new_face.FaceNodes[k]].P;
			};
			new_face.FaceCenter = (1.0/new_face.FaceNodes.size()) * new_face.FaceCenter;
			if (j == 0) {
				//Bottom border
				new_face.FaceNormal = Vector(0, -1, 0);
				new_face.isExternal = true;
				new_face.BCMarker = 3;
				new_face.FaceCell_1 = i + j*(N-1);
				new_face.FaceCell_2 = -1;
			};
			if (j == M-1) {
				//Top border
				new_face.FaceNormal = Vector(0, 1, 0);
				new_face.isExternal = true;
				new_face.BCMarker = 4;
				new_face.FaceCell_1 = i + (j-1)*(N-1);
				new_face.FaceCell_2 = -1;
			};
			if ((j!=0) && (j!=M-1)) {
				//Inner face
				new_face.FaceNormal = Vector(0, 1, 0);
				new_face.isExternal = false;
				new_face.BCMarker = -1;
				new_face.FaceCell_1 = i + (j-1)*(N-1);
				new_face.FaceCell_2 = i + j*(N-1);
			};
			//new_face.FaceNormal = new_face.FaceNormal * new_face.FaceSquare;
			g.faces.add(new_face);
		};
	};		

	//Fill in grid info

	//Add boundary names
	g.addPatch("left", 1);
	g.addPatch("right", 2);
	g.addPatch("bottom", 3);
	g.addPatch("top", 4);
	g.ConstructAndCheckPatches();
	
	return g;
};

#endif
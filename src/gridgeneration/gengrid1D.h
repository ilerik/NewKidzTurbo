#ifndef TURBO_GRIDGENERATION_GENGRID1D
#define TURBO_GRIDGENERATION_GENGRID1D

#include "grid.h"

Grid GenGrid1D(int N, double lBegin, double lEnd, Vector direction)
{
	Grid g;

	//Normalize direction
	direction = direction / direction.mod();

	//Grid parameters
	g.GridID = 1;
	g.gridInfo.GridDimensions = 1;
	g.gridInfo.CellDimensions = 1;
	g.gridInfo.nCoords = 3;

	//Generate grid topology	
	N++;
	std::vector<double> x_p(N);
	std::vector<double> y_p(N);
	std::vector<double> z_p(N);
	double L = lEnd - lBegin;	
	for (int i = 0; i<N; i++) {
		double curL = lBegin + L * i / (N-1);
		Vector r = curL * direction;
		x_p[i] = r.x;
		y_p[i] = r.y;
		z_p[i] = r.z;
	};
		
	//Fill in nodes
	for (int i = 0; i<N; i++) {		
		Node new_node;
		new_node.GlobalIndex = i;
		new_node.P.x = x_p[i];
		new_node.P.y = y_p[i];
		new_node.P.z = z_p[i];		
		//Add node
		g.nodes.add(new_node); 		
	};	

	//Fill in cells
	for (int i = 0; i<N-1; i++) {		
		Cell c;
		c.GlobalIndex = i;
		c.Nodes.resize(2);
		c.Nodes[0] = i;
		c.Nodes[1] = i + 1;	
		c.Faces.resize(2);
		c.Faces[0] = i;
		c.Faces[1] = (i+1);		
		c.CellHSize = 0;
		c.CGNSType = BAR_2;
		g.ComputeGeometricProperties(c);				
		//Add cell
		g.cells.add(c);		
	};

	//Generate faces
	for (int i = 0; i<N; i++) {		
		Face new_face;		
		new_face.GlobalIndex = i;
		new_face.FaceNodes.clear();
		new_face.FaceNodes.push_back(i);			
		new_face.FaceCell_1 = i-1;
		new_face.FaceCell_2 = i;
		new_face.isExternal = false;
		if (i == 0) {
			//Left border			
			new_face.isExternal = true;
			new_face.BCMarker = 1;			
			new_face.FaceCell_1 = 0;
			new_face.FaceCell_2 = -1;
		};
		if (i == (N-1)) {
			//Right border			
			new_face.isExternal = true;
			new_face.BCMarker = 2;
			new_face.FaceCell_1 = N-2;
			new_face.FaceCell_2 = -1;
		};		
		new_face.CGNSType = NODE;
		g.ComputeGeometricProperties(new_face);
		g.faces.add(new_face);		
	};

	//Add boundary names
	g.addPatch("left", 1);
	g.addPatch("right", 2);	
	if (!g.ConstructAndCheckPatches()) {
		std::cout<<"Grid generation failed";
		exit(0);
	};
	
	return g;
};

#endif
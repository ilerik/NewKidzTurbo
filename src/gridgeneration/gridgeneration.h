#ifndef TURBO_GRIDGENERATION
#define TURBO_GRIDGENERATION

#include "gengrid2D.h"
#include "gengrid1D.h"



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


#endif
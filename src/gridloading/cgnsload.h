#ifndef TURBO_GRIDLOADING_CGNSLOAD
#define TURBO_GRIDLOADING_CGNSLOAD

#include "stdlib.h"
#include "grid.h"
#include "cgnslib.h"

int GetElementDimensions(ElementType_t type) {
	int dim = -1;
	switch (type) {
	//0D Elements
	case NODE:
		return 0;
	//1D Elements
	case BAR_2:
		return 1;
	//2D Elements
	case QUAD_4:
		return 2;
	//3D Elements
	case HEXA_8:
		return 3;
	}

	if (dim == -1) std::cout<<"Unknown element type\n";
	return dim;
};


Grid LoadCGNSGrid(std::string fname) {
	Grid grid;	

	BCType_t bcType;
	int cgfile, mode = CG_MODE_READ;

	//Open CGNS file
	printf ("opening cgns file <%s> ...\n", fname.c_str());
    fflush (stdout);
	if (cg_open (fname.c_str(), mode, &cgfile)) cg_error_exit();

	/* Determine the of bases in the grid. This example assumes */
	/* one base. However it is allowed to have multiple bases. */
   
	int nBases;
	if(cg_nbases(cgfile, &nBases)!= CG_OK) cg_error_exit();
	if(nBases != 1) {
		std::cout<< "This CGNS loader assumes one base\n";
		exit(1);
	};
	int base = 1;

	/* Check the cell and physical dimensions of the base. */
	
	int physDim;
	int cellDim;
	char cgnsName[255];
	if(cg_base_read(cgfile, base, cgnsName, &cellDim, &physDim) != CG_OK) cg_error_exit();

	// grid info
	grid.gridInfo.CellDimensions = cellDim;
	grid.gridInfo.GridDimensions = physDim;


	/* Read the number of zones in the grid. */
	/* This example assumes one zone. */

	int nZones;
	if(cg_nzones(cgfile, base, &nZones) != CG_OK) cg_error_exit();     
	if(nZones != 1) {
		std::cout<< "This CGNS loader assumes one zone\n";
		exit(1);
	}
	int zone = 1;

	/* Check the zone type. This should be Unstructured. */

	ZoneType_t zoneType;
	if(cg_zone_type(cgfile, base, zone, &zoneType) != CG_OK) cg_error_exit(); 		
	if(zoneType != Unstructured) {
		std::cout<< "Unstructured zone expected\n";
		exit(1);
	};

	/* Determine the number of vertices and volume elements in this */
	/* zone (and thus in the grid, because one zone is assumed). */

	char zoneName[255];
	cgsize_t sizes[3];
	if(cg_zone_read(cgfile, base, zone, zoneName, sizes) != CG_OK) cg_error_exit();		
	int nVertices    = sizes[0];
	int nVolElements = sizes[1];

	/* Determine the number and names of the coordinates. */

	int nCoords;
	if(cg_ncoords(cgfile, base, zone, &nCoords) != CG_OK)
		cg_error_exit();	

	char name[255];
	DataType_t dataType;
	if(cg_coord_info(cgfile, base, zone, 1, &dataType, name) != CG_OK)
		cg_error_exit();	

	/* Read the x-coordinates. The y and z-coordinates can be read */
	/* similarly. Just replace CoordinateX by CoordinateY and */
	/* CoordinateZ respectively. This assumes Cartesian coordinates */
	/* in double precision. Note that CGNS starts the numbering at */
	/* 1 even if C is used. */
	
	int one = 1;	
	std::vector<double> coorX(nVertices);
	if(cg_coord_read(cgfile, base, zone, "CoordinateX", RealDouble, &one,
					&nVertices, &coorX[0]) != CG_OK) cg_error_exit();	
	
	//Check if we have Y coordinate
	std::vector<double> coorY(nVertices);
	for (int i = 0; i<nVertices; i++) coorY[i] = 0;
	if (nCoords > 1) {
		if(cg_coord_read(cgfile, base, zone, "CoordinateY", RealDouble, &one,
					&nVertices, &coorY[0]) != CG_OK) cg_error_exit();	
	};
	
	//Check if its not 3D case and make all Z equal zero
	std::vector<double> coorZ(nVertices);
	for (int i = 0; i<nVertices; i++) coorZ[i] = 0;
	if (nCoords > 2) {
		if(cg_coord_read(cgfile, base, zone, "CoordinateZ", RealDouble, &one,
					&nVertices, &coorZ[0]) != CG_OK) cg_error_exit();	
	};

	//Create grid nodes
	for (int i = 0; i<nVertices; i++) {
		Node newNode;
		newNode.GlobalIndex = i+1;
		newNode.P.x = coorX[i];
		newNode.P.y = coorY[i];
		newNode.P.z = coorZ[i];
		grid.nodes.add(newNode);
	};

	/* Determine the number of sections for this zone. Note that */
	/* surface elements can be stored in a volume zone, but they */
	/* are NOT taken into account in the number obtained from */
	/* cg_zone_read. */

	int nSections;
	if(cg_nsections(cgfile, base, zone, &nSections) != CG_OK)
		cg_error_exit();	


	/* Loop over the number of sections and read the element */
	/* connectivities. As CGNS starts the numbering at 1 the */
	/* for-loop starts at 1 as well. */


	std::map<std::set<int>, int> faces;
	std::map<int, std::set<int> > indexToElementNodes;
	std::map<int, ElementType_t > indexToElementType;
	for(int sec=1; sec<=nSections; sec++)
	{
		/* Determine the element type and set the pointer for the */
		/* connectivity accordingly. */

		char secName[255];
		ElementType_t type;
		cgsize_t eBeg;
		cgsize_t eEnd;
		int nBdry;
		int parentFlag;
		if(cg_section_read(cgfile, base, zone, sec, secName, &type,
						&eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
			cg_error_exit();
		int nElements = (eEnd-eBeg+1);

		int nNodes;
		int* conn;
		switch (type)
		{
		case QUAD_4: 
			nNodes = 4;			
			break;	
		case BAR_2:
			nNodes = 2;			
			break;	
		case HEXA_8:
			nNodes = 8;
			break;

		default:
			std::cout<< "Unsupported element encountered.\n";
			break;
		}

		conn = new int[nNodes*nElements] ; 			
   
		/* Read the connectivity. Again, the node numbering of the */
		/* connectivities start at 1. If internally a starting index */
		/* of 0 is used (typical for C-codes) 1 must be substracted */
		/* from the connectivities read. */
   		
		if(cg_elements_read(cgfile, base, zone, sec, conn, NULL) != CG_OK)
			cg_error_exit();	

		//Add all elements read to approriate structures
		int counter = 0;
		for (int ind = eBeg; ind<=eEnd; ind++) {
			indexToElementType[ind] = type;
			for (int j = 0; j<nNodes; j++) indexToElementNodes[ind].insert(conn[counter*nNodes + j]);
			counter++;
		};

		//Create cells if its main section (contains elements with dim = physDim)
		//For now we assume that all volume elements are in one section so following condition must hold
		bool isVolumeElementsSection = false;
		if ((nVolElements == nElements) && (GetElementDimensions(type) == grid.gridInfo.CellDimensions)) {
			isVolumeElementsSection = true;
		};

		if (!isVolumeElementsSection) continue;
	
		int faceGlobalIndex = 1;
		counter = 0;
		for (int ind = eBeg; ind<=eEnd; ind++) {
			Cell newCell;
			newCell.GlobalIndex = ind;
			for (int j = 0; j<nNodes; j++) newCell.Nodes.push_back(conn[counter*nNodes + j]);			
			counter++;
			newCell.CGNSType = type;	
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

		delete[] conn;
	};

	/* Determine the number of boundary conditions for this zone. */

	int nBocos;
	if(cg_nbocos(cgfile, base, zone, &nBocos) != CG_OK)
		cg_error_exit();

	/* Loop over the number of boundary conditions. */

	for(int boco=1; boco<=nBocos; boco++)
	{
		/* Read the info for this boundary condition. */
		char bocoName[255];
		BCType_t bocoType;
		PointSetType_t ptsetType;
		int normalIndex;
		cgsize_t nBCElem;
		cgsize_t normListFlag;
		DataType_t normDataType;
		int nDataSet;
		if(cg_boco_info(cgfile, base, zone, boco, bocoName, &bocoType, &ptsetType,
						&nBCElem, &normalIndex, &normListFlag, &normDataType,
						&nDataSet) != CG_OK) cg_error_exit();

		/* Read the element ID's. */

		cgsize_t* BCElemRead = new cgsize_t[nBCElem];
		if(cg_boco_read(cgfile, base, zone, boco, BCElemRead, NULL) != CG_OK) cg_error_exit();

		/* And much more to make it fit into the */
		/* internal data structures. */

		//Assign to all faces boundary marker
		bool notBoundary = false;
		for (int i = 0; i<nBCElem; i++) {
			std::set<int> idx = indexToElementNodes[BCElemRead[i]];
			if (faces.find(idx) != faces.end()) {
				//Element is a face
				int globalInd = faces[idx];
				grid.faces[globalInd].BCMarker = boco;
				if (!grid.faces[globalInd].isExternal) throw Exception("Boundary face must be marked isExternal");
			} else {
				//Its not a face so this is not valid Boundary condition
				//??
				notBoundary = true;
			};
		};

		//Add new patch to grid
		if (!notBoundary) grid.addPatch(bocoName, boco);		
	}

	//Make all boundary data consistent
	grid.ConstructAndCheckPatches();

	//Close CGNS file
	printf ("closing cgns file ...\n");
    fflush (stdout);
	cg_close (cgfile);

	return grid;
};

void ReadGridInfo(Grid& grid, std::string fname) {
	BCType_t bcType;
	int cgfile, mode = CG_MODE_READ;

	//Open CGNS file
	//printf ("opening cgns file <%s> ...\n", fname.c_str());
    //fflush (stdout);
	if (cg_open (fname.c_str(), mode, &cgfile)) cg_error_exit();

	/* Determine the of bases in the grid. This example assumes */
	/* one base. However it is allowed to have multiple bases. */
   
	int nBases;
	if(cg_nbases(cgfile, &nBases)!= CG_OK) cg_error_exit();
	if(nBases != 1) {
		std::cout<< "This CGNS loader assumes one base\n";
		exit(1);
	};
	int base = 1;

	/* Check the cell and physical dimensions of the base. */
	
	int physDim;
	int cellDim;
	char cgnsName[255];
	if(cg_base_read(cgfile, base, cgnsName, &cellDim, &physDim) != CG_OK) cg_error_exit();

	// grid info
	grid.gridInfo.CellDimensions = cellDim;
	grid.gridInfo.GridDimensions = physDim;


	/* Read the number of zones in the grid. */
	/* This example assumes one zone. */

	int nZones;
	if(cg_nzones(cgfile, base, &nZones) != CG_OK) cg_error_exit();     
	if(nZones != 1) {
		std::cout<< "This CGNS loader assumes one zone\n";
		exit(1);
	}
	int zone = 1;

	/* Check the zone type. This should be Unstructured. */

	ZoneType_t zoneType;
	if(cg_zone_type(cgfile, base, zone, &zoneType) != CG_OK) cg_error_exit(); 		
	if(zoneType != Unstructured) {
		std::cout<< "Unstructured zone expected\n";
		exit(1);
	};

	/* Determine the number of vertices and volume elements in this */
	/* zone (and thus in the grid, because one zone is assumed). */

	char zoneName[255];
	cgsize_t sizes[3];
	if(cg_zone_read(cgfile, base, zone, zoneName, sizes) != CG_OK) cg_error_exit();		
	int nVertices    = sizes[0];
	int nVolElements = sizes[1];

	/* Determine the number and names of the coordinates. */

	int nCoords;
	if(cg_ncoords(cgfile, base, zone, &nCoords) != CG_OK)
		cg_error_exit();	

	char name[255];
	DataType_t dataType;
	if(cg_coord_info(cgfile, base, zone, 1, &dataType, name) != CG_OK)
		cg_error_exit();	

};

#endif


#include "grid.h"
#include "mpi.h"
#include "logger.h"

#define err_t int

////Calculation kernel
class Kernel {
public:
	//Logger
	std::string _logfilename;
	Logger _logger;

	//Grid data
	Grid _grid;

	//Task parameters
	int _rungeKuttaOrder;	
	int _nProcessors, _rank; 

	//Internal storage	
	double* Uconservative;	//Conservative variables
	double* Residual;		//Residual

	//Initialize computational kernel
	err_t Initilize(int *argc, char **argv[]) {
		MPI_Init(argc, argv);		
		MPI_Comm_size(MPI_COMM_WORLD, &_nProcessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

		_logfilename = "kernel.log";
		_logger.InitLogging(_logfilename);
	};

	//Finilize kernel work
	err_t Finalize() {
		MPI_Finalize();							
	};

	//Load grid topology info
	err_t LoadGridTopology(std::string filename) {		
		//Open cgns file
		_grid.gridInfo.fileLocation = filename;
		int cgfile, mode = CG_MODE_READ;		
		_logger.WriteMessage(0, "Opening cgns file"+filename+"\n");		
		if (cg_open (filename.c_str(), mode, &cgfile)) {
			_logger.WriteMessage(-1, "Cannot open grid file\n");		
			cg_error_exit();
		};
				
		//Determine the number of bases in the grid. This example assumes 
		//one base. However it is allowed to have multiple bases. 
   
		int nBases;
		if(cg_nbases(cgfile, &nBases)!= CG_OK) {			
			_logger.WriteMessage(-1, "Cannot read number of bases.\n");		
			cg_error_exit();
		};
		if (nBases != 1) {
			_logger.WriteMessage(-1, "Cannot process grid. Must be one base.\n");
			return 1;
		};
		int base = 1;		

		//Check the cell and physical dimensions of the base.
		int physDim;
		int cellDim;
		char cgnsName[255];
		if(cg_base_read(cgfile, base, cgnsName, &cellDim, &physDim) != CG_OK) {
			_logger.WriteMessage(-1, "Cannot read base info.\n");
			cg_error_exit();
			return 1;
		};		
		_grid.gridInfo.CellDimensions = cellDim;
		_grid.gridInfo.GridDimensions = physDim;


		//Read the number of zones in the grid.
		//This example assumes one zone.
		int nZones;
		if(cg_nzones(cgfile, base, &nZones) != CG_OK) {
			_logger.WriteMessage(-1, "Cannot read number of zones.\n");
			cg_error_exit();     
		};
		if(nZones != 1) {
			_logger.WriteMessage(-1, "This CGNS loader assumes one zone\n");
			return 1;
		}
		int zone = 1;

		//Check the zone type. This should be Unstructured.
		ZoneType_t zoneType;
		if(cg_zone_type(cgfile, base, zone, &zoneType) != CG_OK) {
			_logger.WriteMessage(-1, "Cannot read zone type.\n");
			cg_error_exit(); 		
		};
		if(zoneType != Unstructured) {			
			_logger.WriteMessage(-1, "Unstructured zone expected\n");
			return 1;
		};

		//Determine the number of vertices and volume elements in this */
		//zone (and thus in the grid, because one zone is assumed). */

		char zoneName[255];
		cgsize_t sizes[3];
		if(cg_zone_read(cgfile, base, zone, zoneName, sizes) != CG_OK) cg_error_exit();		
		int nVertices    = sizes[0];
		int nVolElements = sizes[1];

		//Determine the number and names of the coordinates.
		int nCoords;
		if(cg_ncoords(cgfile, base, zone, &nCoords) != CG_OK)
			cg_error_exit();	

		char name[255];
		DataType_t dataType;
		if(cg_coord_info(cgfile, base, zone, 1, &dataType, name) != CG_OK)
		cg_error_exit();	


		//Load connectivity info

		return 0;
	};

	//Partion loaded grid
	err_t PartitionGrid() {
	};

	//Load required geometric grid data
	err_t LoadGridGeometry() {
	};

	//Compute residual
	void ComputeResidual(double *R, double *U){
		
	};

	//Compute convective flux through face
	void ComputeConvectiveFlux() {
	};	

	//Explicit time step
	void ExplicitTimeStep() {
	};

	//Implicit time step
	void ImplicitTimeStep() {
	};

};